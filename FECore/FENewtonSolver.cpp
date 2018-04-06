#include "stdafx.h"
#include "FENewtonSolver.h"
#include "NumCore/NumCore.h"
#include "FENodeReorder.h"
#include "FEModel.h"
#include "FEGlobalMatrix.h"
#include "BFGSSolver.h"
#include "BFGSSolver2.h"
#include "FEBroydenStrategy.h"
#include "log.h"

//-----------------------------------------------------------------------------
// define the parameter list
BEGIN_PARAMETER_LIST(FENewtonSolver, FESolver)
	ADD_PARAMETER2(m_LStol             , FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "lstol"   );
	ADD_PARAMETER2(m_LSmin             , FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "lsmin"   );
	ADD_PARAMETER2(m_LSiter            , FE_PARAM_INT   , FE_RANGE_GREATER_OR_EQUAL(0), "lsiter"  );
	ADD_PARAMETER2(m_maxref            , FE_PARAM_INT   , FE_RANGE_GREATER_OR_EQUAL(0.0), "max_refs");
	ADD_PARAMETER2(m_maxups            , FE_PARAM_INT   , FE_RANGE_GREATER_OR_EQUAL(0.0), "max_ups" );
	ADD_PARAMETER2(m_max_buf_size      , FE_PARAM_INT   , FE_RANGE_GREATER_OR_EQUAL(0), "qn_max_buffer_size");
	ADD_PARAMETER(m_cycle_buffer       , FE_PARAM_BOOL  , "qn_cycle_buffer");
	ADD_PARAMETER2(m_cmax              , FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "cmax"    );
	ADD_PARAMETER(m_nqnmethod          , FE_PARAM_INT   , "qnmethod");
	ADD_PARAMETER(m_bzero_diagonal     , FE_PARAM_BOOL  , "check_zero_diagonal");
	ADD_PARAMETER(m_zero_tol           , FE_PARAM_DOUBLE, "zero_diagonal_tol"  );
	ADD_PARAMETER(m_profileUpdateMethod, FE_PARAM_INT   , "profile_update_method");
	ADD_PARAMETER(m_eq_scheme          , FE_PARAM_INT   , "equation_scheme");
	ADD_PARAMETER(m_force_partition    , FE_PARAM_INT   , "force_partition");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FENewtonSolver::FENewtonSolver(FEModel* pfem) : FESolver(pfem)
{
	// default parameters
	m_LSmin = 0.01;
	m_LStol = 0.9;
	m_LSiter = 5;
	m_maxref = 15;

    m_neq = 0;
    m_plinsolve = 0;
	m_pK = 0;

	m_cmax   = 1e5;
	m_maxups = 10;
	m_max_buf_size = 0;
	m_cycle_buffer = true;
	m_nqnmethod = QN_BFGS;
	m_pbfgs = 0;

	m_profileUpdateMethod = 0;

	m_bzero_diagonal = true;
	m_zero_tol = 0.0;

	m_force_partition = 0;

	m_eq_scheme = EQUATION_SCHEME::STAGGERED;
}

//-----------------------------------------------------------------------------
//! Set the default solution strategy
void FENewtonSolver::SetDefaultStrategy(QN_STRATEGY qn)
{
	m_nqnmethod = qn;
}

//-----------------------------------------------------------------------------
void FENewtonSolver::SetSolutionStrategy(FENewtonStrategy* pstrategy)
{
	if (m_pbfgs) delete m_pbfgs;
	m_pbfgs = pstrategy;
}

//-----------------------------------------------------------------------------
FENewtonSolver::~FENewtonSolver()
{
	delete m_plinsolve;	// clean up linear solver data
	delete m_pK;		// clean up stiffnes matrix data
	if (m_pbfgs) delete m_pbfgs;
}

//-----------------------------------------------------------------------------
FEGlobalMatrix& FENewtonSolver::GetStiffnessMatrix()
{
	return *m_pK;
}

//-----------------------------------------------------------------------------
//! Check the zero diagonal
void FENewtonSolver::CheckZeroDiagonal(bool bcheck, double ztol)
{
	m_bzero_diagonal = bcheck;
	m_zero_tol = fabs(ztol);
}

//-----------------------------------------------------------------------------
//! Reforms a stiffness matrix and factorizes it
bool FENewtonSolver::ReformStiffness(const FETimeInfo& tp)
{
	felog.printf("Reforming stiffness matrix: reformation #%d\n\n", m_nref);

    // first, let's make sure we have not reached the max nr of reformations allowed
    if (m_nref >= m_maxref) throw MaxStiffnessReformations();
    
    // recalculate the shape of the stiffness matrix if necessary
    if (m_breshape)
    {
        // TODO: I don't think I need to update here
        //		if (m_fem.m_bcontact) UpdateContact();
        
        // reshape the stiffness matrix
        if (!CreateStiffness(m_niter == 0)) return false;
        
        // reset reshape flag, except for contact
		m_breshape = (((m_fem.SurfacePairConstraints() > 0) || (m_fem.NonlinearConstraints() > 0)) ? true : false);
    }
    
    // calculate the global stiffness matrix
	bool bret = false;
	{
		TRACK_TIME("stiffness");
	    bret = StiffnessMatrix(tp);

		// check for zero diagonals
		if (m_bzero_diagonal)
		{
			// get the stiffness matrix
			SparseMatrix& K = *m_pK;
			vector<int> zd;
			int neq = K.Size();
			for (int i=0; i<neq; ++i)
			{
				double di = fabs(K.diag(i));
				if (di <= m_zero_tol) zd.push_back(i);
			}

			if (zd.empty() == false) throw ZeroDiagonal(-1, -1);
		}
	}

	// if the stiffness matrix was evaluated successfully,
	// we factor it.
    if (bret)
    {
        {
			TRACK_TIME("solve");
			// factorize the stiffness matrix
            m_plinsolve->Factor();
        }
        
        // increase total nr of reformations
        m_nref++;
        m_ntotref++;
        
        // reset bfgs update counter
        m_pbfgs->m_nups = 0;
    }
    
    return bret;
}

//-----------------------------------------------------------------------------
//!  Creates the global stiffness matrix
//! \todo Can we move this to the FEGlobalMatrix::Create function?
bool FENewtonSolver::CreateStiffness(bool breset)
{
	{
		TRACK_TIME("reform");
		// clean up the solver
		if (m_pK->NonZeroes()) m_plinsolve->Destroy();

		// clean up the stiffness matrix
		m_pK->Clear();

		// create the stiffness matrix
		felog.printf("===== reforming stiffness matrix:\n");
		SparseMatrixProfile::UpdateMethod updateMethod = (m_profileUpdateMethod == 0 ? SparseMatrixProfile::Method1 : SparseMatrixProfile::Method2);
		if (m_pK->Create(&GetFEModel(), m_neq, breset, updateMethod) == false) 
		{
			felog.printf("FATAL ERROR: An error occured while building the stiffness matrix\n\n");
			return false;
		}
		else
		{
			// output some information about the direct linear solver
			int neq = m_pK->Rows();
			int nnz = m_pK->NonZeroes();
			felog.printf("\tNr of equations ........................... : %d\n", neq);
			felog.printf("\tNr of nonzeroes in stiffness matrix ....... : %d\n", nnz);
			felog.printf("\n");
		}
		// let's flush the logfile to make sure the last output will not get lost
		felog.flush();
	}

	// Do the preprocessing of the solver
	{
		TRACK_TIME("solve");
		if (!m_plinsolve->PreProcess())
		{
			// TODO: get rid of throwing this exception. We should just return false.
			throw FatalError();
		}
	}

	// done!
	return true;
}

//-----------------------------------------------------------------------------
bool FENewtonSolver::Init()
{
	// Base class initialization and validation
	if (FESolver::Init() == false) return false;

	// choose a solution strategy
	switch (m_nqnmethod)
	{
	case QN_BFGS   : SetSolutionStrategy(new BFGSSolver ); break;
	case QN_BROYDEN: SetSolutionStrategy(new FEBroydenStrategy); break;
	// NOTE: Temporary hack for backward compatibility since the BFGSSolver2 was deprecated
	// This solver used to have the value 1 and Broyden 2, but 1 is now used for Broyden.
	case 2: SetSolutionStrategy(new FEBroydenStrategy); break;
	default:
		return false;
	}

	// set the solution parameters
	m_pbfgs->m_maxups = m_maxups;
	m_pbfgs->m_max_buf_size = m_max_buf_size;
	m_pbfgs->m_cycle_buffer = m_cycle_buffer;
	m_pbfgs->m_cmax   = m_cmax;

    // Now that we have determined the equation numbers we can continue
    // with creating the stiffness matrix. First we select the linear solver
    // The stiffness matrix is created in CreateStiffness
    // Note that if a particular solver was requested in the input file
    // then the solver might already be allocated. That's way we need to check it.
    if (m_plinsolve == 0)
    {
		FECoreKernel& fecore = FECoreKernel::GetInstance();
		m_plinsolve = fecore.CreateLinearSolver(m_fem.GetLinearSolverType());
        if (m_plinsolve == 0)
        {
            felog.printbox("FATAL ERROR","Unknown solver type selected\n");
            return false;
        }
    }

	// allocate data vectors
	m_R0.assign(m_neq, 0);
	m_R1.assign(m_neq, 0);
	m_ui.assign(m_neq, 0);

	// initialize BFGS data
	// Must be done after initialization of linear solver
	m_pbfgs->Init(m_neq, m_plinsolve);

    // set the create stiffness matrix flag
    m_breshape = true;

	// allocate storage for the sparse matrix that will hold the stiffness matrix data
	// we let the linear solver allocate the correct type of matrix format
	SparseMatrix* pS = m_plinsolve->CreateSparseMatrix(m_bsymm? REAL_SYMMETRIC : REAL_UNSYMMETRIC);
	if (pS == 0)
	{
		felog.printbox("FATAL ERROR", "The selected linear solver does not support the requested\n matrix format.\nPlease select a different linear solver.\n");
		return false;
	}

	// clean up the stiffness matrix if we have one
	if (m_pK) delete m_pK; m_pK = 0;

	// Create the stiffness matrix.
	// Note that this does not construct the stiffness matrix. This
	// is done later in the CreateStiffness routine.
	m_pK = new FEGlobalMatrix(pS);
	if (m_pK == 0)
	{
		felog.printbox("FATAL ERROR", "Failed allocating stiffness matrix\n\n");
		return false;
	}

	// Set the partitioning of the global matrix
	// This is only used for debugging block solvers for problems that
	// usually don't generate a block structure
	if (m_force_partition > 0) m_plinsolve->SetPartition(m_force_partition);

	return true;
}

//-----------------------------------------------------------------------------
//!	This function initializes the equation system.
//! It is assumed that all free dofs up until now have been given an ID >= 0
//! and the fixed or rigid dofs an ID < 0.
//! After this operation the nodal ID array will contain the equation
//! number assigned to the corresponding degree of freedom. To distinguish
//! between free or unconstrained dofs and constrained ones the following rules
//! apply to the ID array:
//!
//!           /
//!          |  >=  0 --> dof j of node i is a free dof
//! ID[i][j] <  == -1 --> dof j of node i is a fixed (no equation assigned too)
//!          |  <  -1 --> dof j of node i is constrained and has equation nr = -ID[i][j]-2
//!           \
//!
bool FENewtonSolver::InitEquations()
{
    // get the mesh
    FEMesh& mesh = m_fem.GetMesh();
    
    // initialize nr of equations
    int neq = 0;
    
    // see if we need to optimize the bandwidth
    if (m_fem.OptimizeBandwidth())
    {
		assert(m_eq_scheme == EQUATION_SCHEME::STAGGERED);
        // reorder the node numbers
        vector<int> P(mesh.Nodes());
        FENodeReorder mod;
        mod.Apply(mesh, P);
        
        // set the equation numbers
        for (int i=0; i<mesh.Nodes(); ++i)
        {
            FENode& node = mesh.Node(P[i]);
            for (int j=0; j<(int)node.m_ID.size(); ++j)
            {
                if      (node.m_ID[j] == DOF_FIXED     ) { node.m_ID[j] = -1; }
                else if (node.m_ID[j] == DOF_OPEN      ) { node.m_ID[j] =  neq++; }
                else if (node.m_ID[j] == DOF_PRESCRIBED) { node.m_ID[j] = -neq-2; neq++; }
                else { assert(false); return false; }
            }
        }
    }
    else
    {
		if (m_eq_scheme == EQUATION_SCHEME::STAGGERED)
		{
			// give all free dofs an equation number
			for (int i=0; i<mesh.Nodes(); ++i)
			{
				FENode& node = mesh.Node(i);
				for (int j=0; j<(int)node.m_ID.size(); ++j)
				{
					if      (node.m_ID[j] == DOF_FIXED     ) { node.m_ID[j] = -1; }
					else if (node.m_ID[j] == DOF_OPEN      ) { node.m_ID[j] =  neq++; }
					else if (node.m_ID[j] == DOF_PRESCRIBED) { node.m_ID[j] = -neq-2; neq++; }
					else { assert(false); return false; }
				}
			}
		}
		else
		{
			assert(m_eq_scheme == EQUATION_SCHEME::BLOCK);

			// Assign equations numbers in blocks
			DOFS& dofs = m_fem.GetDOFS();
			for (int nv=0; nv<dofs.Variables(); ++nv)
			{
				int n = dofs.GetVariableSize(nv);
				for (int l=0; l<n; ++l)
				{
					int nl = dofs.GetDOF(nv, l);

					for (int i = 0; i<mesh.Nodes(); ++i)
					{
						FENode& node = mesh.Node(i);
						if      (node.m_ID[nl] == DOF_FIXED     ) { node.m_ID[nl] = -1; }
						else if (node.m_ID[nl] == DOF_OPEN      ) { node.m_ID[nl] = neq++; }
						else if (node.m_ID[nl] == DOF_PRESCRIBED) { node.m_ID[nl] = -neq - 2; neq++; }
						else { assert(false); return false; }
					}
				}
			}
		}
    }
    
    // store the number of equations
    m_neq = neq;
    
    // All initialization is done
    return true;
}

//-----------------------------------------------------------------------------
//! Clean
//! \todo Why can this not be done in destructor?
void FENewtonSolver::Clean()
{
	if (m_plinsolve) m_plinsolve->Destroy();
}

//-----------------------------------------------------------------------------
void FENewtonSolver::Serialize(DumpStream& ar)
{
	FESolver::Serialize(ar);

	if (ar.IsSaving())
	{
		ar << m_neq;
		ar << m_LStol << m_LSiter << m_LSmin;
		ar << m_maxref;

		if (m_pbfgs == 0) ar << 0; else ar << 1;
		if (m_pbfgs)
		{
			ar << m_nqnmethod;
			ar << m_pbfgs->m_maxups;
			ar << m_pbfgs->m_max_buf_size;
			ar << m_pbfgs->m_cycle_buffer;
			ar << m_pbfgs->m_cmax;
			ar << m_pbfgs->m_nups;
		}
	}
	else
	{
		ar >> m_neq;
		ar >> m_LStol >> m_LSiter >> m_LSmin;
		ar >> m_maxref;

		int n = -1;
		ar >> n;
		if (n)
		{
			ar >> m_nqnmethod;
			switch (m_nqnmethod)
			{
			case QN_BFGS: SetSolutionStrategy(new BFGSSolver); break;
			case QN_BROYDEN: SetSolutionStrategy(new FEBroydenStrategy); break;
			// NOTE: Temporary hack for backward compatibility since the BFGSSolver2 was deprecated
			// This solver used to have the value 1 and Broyden 2, but 1 is now used for Broyden.
			case 2: SetSolutionStrategy(new FEBroydenStrategy); break;
			default:
				return;
			}

			ar >> m_pbfgs->m_maxups;
			ar >> m_pbfgs->m_max_buf_size;
			ar >> m_pbfgs->m_cycle_buffer;
			ar >> m_pbfgs->m_cmax;
			ar >> m_pbfgs->m_nups;
		}

		// realloc data
		m_R0.assign(m_neq, 0);
		m_R1.assign(m_neq, 0);
		m_ui.assign(m_neq, 0);
	}
}

//-----------------------------------------------------------------------------
//!  This function mainly calls the Quasin routine 
//!  and deals with exceptions that require the immediate termination of
//!	quasi-Newton iterations.
bool FENewtonSolver::SolveStep(double time)
{
	bool bret;

	try
	{
		// let's try to call Quasin
		bret = Quasin(time);
	}
	catch (NegativeJacobian e)
	{
		// A negative jacobian was detected
		felog.printbox("ERROR","Negative jacobian was detected at element %d at gauss point %d\njacobian = %lg\n", e.m_iel, e.m_ng+1, e.m_vol);
		return false;
	}
	catch (MaxStiffnessReformations)
	{
		// max nr of reformations is reached
		felog.printbox("ERROR", "Max nr of reformations reached.");
		return false;
	}
	catch (ForceConversion)
	{
		// user forced conversion of problem
		felog.printbox("WARNING", "User forced conversion.\nSolution might not be stable.");
		return true;
	}
	catch (IterationFailure)
	{
		// user caused a forced iteration failure
		felog.printbox("WARNING", "User forced iteration failure.");
		return false;
	}
	catch (ZeroLinestepSize)
	{
		// a zero line step size was detected
		felog.printbox("ERROR", "Zero line step size.");
		return false;
	}
	catch (EnergyDiverging)
	{
		// problem was diverging after stiffness reformation
		felog.printbox("ERROR", "Problem diverging uncontrollably.");
		return false;
	}
	catch (FEMultiScaleException e)
	{
		// the RVE problem didn't solve
		// logging was turned off during multi-scale runs
		// so we need to turn it back on
		felog.SetMode(Logfile::LOG_SCREEN);
		felog.printbox("ERROR", "The RVE problem has failed at element %d, gauss point %d.\nAborting macro run.", e.elemId, e.gptIndex+1);

		return false;
	}
	catch (DoRunningRestart)
	{
		// a request to fail the iteration and restart the time step
		return false;
	}

	if (bret)
	{
		// print a convergence summary to the felog file
		Logfile::MODE mode = felog.GetMode();
		if (mode != Logfile::LOG_NEVER)
		{
			felog.SetMode(Logfile::LOG_FILE);
			felog.printf("\nconvergence summary\n");
			felog.printf("    number of iterations   : %d\n", m_niter);
			felog.printf("    number of reformations : %d\n", m_nref);
			felog.SetMode(mode);
		}
	}

	return bret;
}

//-----------------------------------------------------------------------------
//! Solve the linear system of equations.
//! x is the solution vector
//! R is the right-hand-side vector
void FENewtonSolver::SolveLinearSystem(vector<double>& x, vector<double>& R)
{
	// solve the equations
	if (m_plinsolve->BackSolve(x, R) == false)
		throw LinearSolverFailed();
}

//-----------------------------------------------------------------------------
void FENewtonSolver::QNSolve(vector<double>& ui, vector<double>& R)
{
	TRACK_TIME("solve");
	m_pbfgs->SolveEquations(ui, R);
}

//-----------------------------------------------------------------------------
//! Do a QN update
bool FENewtonSolver::QNUpdate(double ls, vector<double>& ui, vector<double>& R0, vector<double>& R1)
{
	TRACK_TIME("qn_update");

	// make sure we didn't reach max updates
	if (m_pbfgs->m_nups >= m_pbfgs->m_maxups - 1)
	{
		// print a warning only if the user did not intent full-Newton
		if (m_pbfgs->m_maxups > 0)
			felog.printbox("WARNING", "Max nr of iterations reached.\nStiffness matrix will now be reformed.");
		return false;
	}

	// try to do an update
	bool bret = m_pbfgs->Update(ls, m_ui, m_R0, m_R1);

	if (bret == false)
	{
		// Stiffness update has failed.
		// this might be due a too large condition number
		// or the update was no longer positive definite.
		felog.printbox("WARNING", "The QN update has failed.\nStiffness matrix will now be reformed.");
	}

	return bret;
}

//-----------------------------------------------------------------------------
//! Performs a linesearch on a NR iteration
//! The description of this method can be found in:
//!    "Nonlinear Continuum Mechanics for Finite Element Analysis", 	Bonet & Wood.
//
//! \todo Find a different way to update the deformation based on the ls.
//! For instance, define a di so that ui = s*di. Also, define the 
//! position of the nodes at the previous iteration.

double FENewtonSolver::LineSearch(double s)
{
	double smin = s;

	double a, A, B, D;
	double r0, r1, r;

	// max nr of line search iterations
	int nmax = m_LSiter;
	int n = 0;

	// initial energy
	r0 = m_ui*m_R0;

	double rmin = fabs(r0);

	// ul = ls*ui
	vector<double> ul(m_ui.size());
	do
	{
		// Update geometry
		vcopys(ul, m_ui, s);
		Update(ul);

		// calculate residual at this point
		Residual(m_R1);

		// make sure we are still in a valid range
		if (s < m_LSmin)
		{
			// it appears that we are not converging
			// I found in the NIKE3D code that when this happens,
			// the line search step is simply set to 0.5.
			// so let's try it here too
			s = 0.5;

			// reupdate  
			vcopys(ul, m_ui, s);
			Update(ul);

			// recalculate residual at this point
			Residual(m_R1);

			// return and hope for the best
			break;
		}

		// calculate energies
		r1 = m_ui*m_R1;

		if ((n == 0) || (fabs(r1) < rmin))
		{
			smin = s;
			rmin = fabs(r1);
		}

		// make sure that r1 does not happen to be really close to zero,
		// since in that case we won't find any better solution.
		if (fabs(r1) < 1.e-17) r = 0;
		else r = fabs(r1 / r0);

		if (r > m_LStol)
		{
			// calculate the line search step
			a = r0 / r1;

			A = 1 + a*(s - 1);
			B = a*(s*s);
			D = B*B - 4 * A*B;

			// update the line step
			if (D >= 0)
			{
				s = (B + sqrt(D)) / (2 * A);
				if (s < 0) s = (B - sqrt(D)) / (2 * A);
				if (s < 0) s = 0;
			}
			else
			{
				s = 0.5*B / A;
			}

			++n;
		}
	} while ((r > m_LStol) && (n < nmax));

	if (n >= nmax)
	{
		// max nr of iterations reached.
		// we choose the line step that reached the smallest energy
		s = smin;
		vcopys(ul, m_ui, s);
		Update(ul);
		Residual(m_R1);
	}
	return s;
}
