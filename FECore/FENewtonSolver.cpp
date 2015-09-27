#include "stdafx.h"
#include "FENewtonSolver.h"
#include "NumCore/NumCore.h"
#include "FECore/FENodeReorder.h"
#include "FEModel.h"
#include "log.h"

//-----------------------------------------------------------------------------
// define the parameter list
BEGIN_PARAMETER_LIST(FENewtonSolver, FESolver)
	ADD_PARAMETER(m_LStol        , FE_PARAM_DOUBLE, "lstol"   );
	ADD_PARAMETER(m_LSmin        , FE_PARAM_DOUBLE, "lsmin"   );
	ADD_PARAMETER(m_LSiter       , FE_PARAM_INT   , "lsiter"  );
	ADD_PARAMETER(m_bfgs.m_maxref, FE_PARAM_INT   , "max_refs");
	ADD_PARAMETER(m_bfgs.m_maxups, FE_PARAM_INT   , "max_ups" );
	ADD_PARAMETER(m_bfgs.m_cmax  , FE_PARAM_DOUBLE, "cmax"    );
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FENewtonSolver::FENewtonSolver(FEModel* pfem) : FESolver(pfem)
{
	// default line search parameters
	m_LSmin = 0.01;
	m_LStol = 0.9;
	m_LSiter = 5;

    m_neq = 0;
    m_plinsolve = 0;
}

//-----------------------------------------------------------------------------
FENewtonSolver::~FENewtonSolver()
{
	delete m_plinsolve;	// clean up linear solver data
}

//-----------------------------------------------------------------------------
bool FENewtonSolver::Init()
{
	// check parameters
	if (m_LStol  < 0.0) { felog.printf("Error: lstol must be nonnegative.\n" ); return false; }
	if (m_LSmin  < 0.0) { felog.printf("Error: lsmin must be nonnegative.\n" ); return false; }
	if (m_LSiter < 0  ) { felog.printf("Error: lsiter must be nonnegative.\n"  ); return false; }
	if (m_bfgs.m_maxref < 0) { felog.printf("Error: max_refs must be nonnegative.\n"); return false; }
	if (m_bfgs.m_maxups < 0) { felog.printf("Error: max_ups must be nonnegative.\n" ); return false; }
	if (m_bfgs.m_cmax   < 0) { felog.printf("Error: cmax must be nonnegative.\n"    ); return false; }

    // Now that we have determined the equation numbers we can continue
    // with creating the stiffness matrix. First we select the linear solver
    // The stiffness matrix is created in CreateStiffness
    // Note that if a particular solver was requested in the input file
    // then the solver might already be allocated. That's way we need to check it.
    if (m_plinsolve == 0)
    {
        m_plinsolve = NumCore::CreateLinearSolver(m_fem.m_nsolver);
        if (m_plinsolve == 0)
        {
            felog.printbox("FATAL ERROR","Unknown solver type selected\n");
            return false;
        }
    }

	// initialize BFGS data
	// Must be done after initialization of linear solver
	m_bfgs.Init(m_neq, this, m_plinsolve);

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
    if (m_fem.m_bwopt)
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
void FENewtonSolver::Serialize(DumpFile& ar)
{
	FESolver::Serialize(ar);

	if (ar.IsSaving())
	{
		ar << m_neq;
		ar << m_LStol << m_LSiter << m_LSmin;
		ar << m_bfgs.m_maxups;
		ar << m_bfgs.m_maxref;
		ar << m_bfgs.m_cmax;
		ar << m_bfgs.m_nups;
	}
	else
	{
		ar >> m_neq;
		ar >> m_LStol >> m_LSiter >> m_LSmin;
		ar >> m_bfgs.m_maxups;
		ar >> m_bfgs.m_maxref;
		ar >> m_bfgs.m_cmax;
		ar >> m_bfgs.m_nups;
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
	catch (FEMultiScaleException)
	{
		// the RVE problem didn't solve
		felog.printbox("ERROR", "The RVE problem has failed. Aborting macro run.");
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
	r0 = m_bfgs.m_ui*m_bfgs.m_R0;

	double rmin = fabs(r0);

	// ul = ls*ui
	vector<double> ul(m_bfgs.m_ui.size());
	do
	{
		// Update geometry
		vcopys(ul, m_bfgs.m_ui, s);
		Update(ul);

		// calculate residual at this point
		Evaluate(m_bfgs.m_R1);

		// make sure we are still in a valid range
		if (s < m_LSmin)
		{
			// it appears that we are not converging
			// I found in the NIKE3D code that when this happens,
			// the line search step is simply set to 0.5.
			// so let's try it here too
			s = 0.5;

			// reupdate  
			vcopys(ul, m_bfgs.m_ui, s);
			Update(ul);

			// recalculate residual at this point
			Evaluate(m_bfgs.m_R1);

			// return and hope for the best
			break;
		}

		// calculate energies
		r1 = m_bfgs.m_ui*m_bfgs.m_R1;

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
		vcopys(ul, m_bfgs.m_ui, s);
		Update(ul);
		Evaluate(m_bfgs.m_R1);
	}
	return s;
}
