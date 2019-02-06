#include "stdafx.h"
#include "FESolver.h"
#include "FEModel.h"
#include "FENodeReorder.h"

REGISTER_SUPER_CLASS(FESolver, FESOLVER_ID);

BEGIN_FECORE_CLASS(FESolver, FECoreBase)
	ADD_PARAMETER(m_bsymm    , "symmetric_stiffness");
	ADD_PARAMETER(m_eq_scheme, "equation_scheme");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FESolver::FESolver(FEModel* fem) : FECoreBase(fem, FESOLVER_ID)
{ 
	m_bsymm = true; // assume symmetric stiffness matrix
	m_niter = 0;

	m_nref = 0;
	
	m_baugment = false;
	m_naug = 0;

	m_eq_scheme = EQUATION_SCHEME::STAGGERED;
}

//-----------------------------------------------------------------------------
FESolver::~FESolver()
{
}

//-----------------------------------------------------------------------------
void FESolver::SetEquationScheme(EQUATION_SCHEME scheme)
{
	m_eq_scheme = scheme;
}

//-----------------------------------------------------------------------------
//! set the linear system partitions
void FESolver::SetPartitions(const vector<int>& part)
{
	m_part = part;
}

//-----------------------------------------------------------------------------
//! Get the size of a partition
int FESolver::GetPartitionSize(int partition)
{
	assert((partition >= 0) && (partition < (int)m_part.size()));
	if ((partition >= 0) && (partition < (int)m_part.size())) return m_part[partition];
	else return 0;
}

//-----------------------------------------------------------------------------
//! get the current stiffness matrix
FEGlobalMatrix* FESolver::GetStiffnessMatrix()
{
	return nullptr;
}

//-----------------------------------------------------------------------------
//! get the current load vector
std::vector<double> FESolver::GetLoadVector()
{
	return std::vector<double>();
}

//-----------------------------------------------------------------------------
void FESolver::Clean()
{
}

//-----------------------------------------------------------------------------
// get the linear solver
LinearSolver* FESolver::GetLinearSolver()
{
	return nullptr;
}

//-----------------------------------------------------------------------------
//! This function is called right before SolveStep and should be used to initialize
//! time dependent information and other settings.
bool FESolver::InitStep(double time)
{
	FEModel& fem = *GetFEModel();

	// evaluate load controllers values at current time
	fem.EvaluateLoadControllers(time);

	// evaluate load parameters
	fem.EvaluateLoadParameters();

	// re-validate materials
	// This is necessary since the material parameters can have changed (e.g. via load curves) and thus 
	// a new validation needs to be done to see if the material parameters are still valid. 
	if (fem.ValidateMaterials() == false) return false;

	return true;
}

//-----------------------------------------------------------------------------
bool FESolver::Init()
{
	// parameter checking
	return Validate();
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
bool FESolver::InitEquations()
{
   // get the mesh
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
    
    // initialize nr of equations
    int neq = 0;
	m_part.clear();
    
    // see if we need to optimize the bandwidth
    if (fem.OptimizeBandwidth())
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

		// assign partition
		m_part.push_back(neq);
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

			// assign partition
			m_part.push_back(neq);
		}
		else
		{
			assert(m_eq_scheme == EQUATION_SCHEME::BLOCK);

			// Assign equations numbers in blocks
			DOFS& dofs = fem.GetDOFS();
			for (int nv=0; nv<dofs.Variables(); ++nv)
			{
				for (int i = 0; i<mesh.Nodes(); ++i)
				{
					FENode& node = mesh.Node(i);
					int n = dofs.GetVariableSize(nv);
					for (int l=0; l<n; ++l)
					{
						int nl = dofs.GetDOF(nv, l);

						if      (node.m_ID[nl] == DOF_FIXED     ) { node.m_ID[nl] = -1; }
						else if (node.m_ID[nl] == DOF_OPEN      ) { node.m_ID[nl] = neq++; }
						else if (node.m_ID[nl] == DOF_PRESCRIBED) { node.m_ID[nl] = -neq - 2; neq++; }
						else { assert(false); return false; }
					}
				}

				// assign partitions
				if (nv == 0) m_part.push_back(neq);
				else m_part.push_back(neq - m_part[nv - 1]);
			}
		}
    }
    
    // store the number of equations
    m_neq = neq;
    
    // All initialization is done
    return true;
}

//-----------------------------------------------------------------------------
//! add equations
void FESolver::AddEquations(int neq, int partition)
{
	m_neq += neq;
	m_part[partition] += neq;
}

//-----------------------------------------------------------------------------
void FESolver::Serialize(DumpStream& ar)
{
	FECoreBase::Serialize(ar);
	if (ar.IsShallow())
	{
		if (ar.IsSaving())
		{
			ar << m_bsymm;
			ar << m_nrhs << m_niter << m_nref << m_ntotref << m_naug;
		}
		else
		{
			ar >> m_bsymm;
			ar >> m_nrhs >> m_niter >> m_nref >> m_ntotref >> m_naug;
		}
	}
}
