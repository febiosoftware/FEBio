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
	ADD_PARAMETER2(m_LStol    , FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "lstol"   );
	ADD_PARAMETER2(m_LSmin    , FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "lsmin"   );
	ADD_PARAMETER2(m_LSiter   , FE_PARAM_INT   , FE_RANGE_GREATER_OR_EQUAL(0), "lsiter"  );
	ADD_PARAMETER2(m_maxref   , FE_PARAM_INT   , FE_RANGE_GREATER_OR_EQUAL(0.0), "max_refs");
	ADD_PARAMETER2(m_maxups   , FE_PARAM_INT   , FE_RANGE_GREATER_OR_EQUAL(0.0), "max_ups" );
	ADD_PARAMETER2(m_cmax     , FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "cmax"    );
	ADD_PARAMETER(m_nqnsolver, FE_PARAM_INT   , "qnmethod");
	ADD_PARAMETER(m_bzero_diagonal, FE_PARAM_BOOL  , "check_zero_diagonal");
	ADD_PARAMETER(m_zero_tol      , FE_PARAM_DOUBLE, "zero_diagonal_tol"  );
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
	m_nqnsolver = QN_BFGS;
	m_pbfgs = 0;

	m_bzero_diagonal = true;
	m_zero_tol = 0.0;
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
        m_breshape = (((m_fem.SurfacePairInteractions() > 0) || (m_fem.NonlinearConstraints() > 0)) ? true : false);
    }
    
    // calculate the global stiffness matrix
	m_StiffnessTime.start();
    bool bret = StiffnessMatrix(tp);

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
	m_StiffnessTime.stop();

	// if the stiffness matrix was evaluated successfully,
	// we factor it.
    if (bret)
    {
        m_SolverTime.start();
        {
            // factorize the stiffness matrix
            m_plinsolve->Factor();
        }
        m_SolverTime.stop();
        
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
	m_ReformTime.start();
	{
		// clean up the solver
		if (m_pK->NonZeroes()) m_plinsolve->Destroy();

		// clean up the stiffness matrix
		m_pK->Clear();

		// create the stiffness matrix
		felog.printf("===== reforming stiffness matrix:\n");
		if (m_pK->Create(&GetFEModel(), m_neq, breset) == false) 
		{
			felog.printf("FATAL ERROR: An error occured while building the stiffness matrix\n\n");
			m_ReformTime.stop();
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
	m_ReformTime.stop();

	// Do the preprocessing of the solver
	m_SolverTime.start();
	{
		if (!m_plinsolve->PreProcess()) 
		{
			// TODO: get rid of throwing this exception. We should just return false.
			m_SolverTime.stop();
			throw FatalError();
		}
	}
	m_SolverTime.stop();

	// done!
	return true;
}

//-----------------------------------------------------------------------------
bool FENewtonSolver::Init()
{
	// Base class initialization and validation
	if (FESolver::Init() == false) return false;

	// choose a solution strategy
	switch (m_nqnsolver)
	{
	case QN_BFGS   : SetSolutionStrategy(new BFGSSolver ); break;
	case QN_BFGS2  : SetSolutionStrategy(new BFGSSolver2); break;
	case QN_BROYDEN: SetSolutionStrategy(new FEBroydenStrategy); break;
	default:
		return false;
	}

	// set the solution parameters
	m_pbfgs->m_maxups = m_maxups;
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
			ar << m_nqnsolver;
			ar << m_pbfgs->m_maxups;
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
			ar >> m_nqnsolver;
			switch (m_nqnsolver)
			{
			case QN_BFGS: SetSolutionStrategy(new BFGSSolver); break;
			case QN_BFGS2: SetSolutionStrategy(new BFGSSolver2); break;
			case QN_BROYDEN: SetSolutionStrategy(new FEBroydenStrategy); break;
			default:
				return;
			}

			ar >> m_pbfgs->m_maxups;
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
		if (fabs(r1) < 1.e-20) r = 0;
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
