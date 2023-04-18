/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/



#include "stdafx.h"
#include "FESolver.h"
#include "FEModel.h"
#include "FENodeReorder.h"
#include "DumpStream.h"
#include "FEDomain.h"
#include "FESurfacePairConstraint.h"
#include "FENLConstraint.h"
#include "FELinearConstraintManager.h"
#include "FENodalLoad.h"
#include "LinearSolver.h"

BEGIN_FECORE_CLASS(FESolver, FECoreBase)
	BEGIN_PARAM_GROUP("linear system");
		ADD_PARAMETER(m_msymm    , "symmetric_stiffness", 0, "non-symmetric\0symmetric\0symmetric structure\0");
		ADD_PARAMETER(m_eq_scheme, "equation_scheme", 0, "staggered\0block\0");
		ADD_PARAMETER(m_eq_order , "equation_order", 0, "default\0reverse\0febio2\0");
		ADD_PARAMETER(m_bwopt    , "optimize_bw");
	END_PARAM_GROUP();
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FESolver::FESolver(FEModel* fem) : FECoreBase(fem)
{ 
	m_msymm = REAL_SYMMETRIC; // assume symmetric stiffness matrix
	m_niter = 0;

	m_nref = 0;
	
	m_baugment = false;
	m_naug = 0;

	m_neq = 0;

	m_bwopt = false;

	m_eq_scheme = EQUATION_SCHEME::STAGGERED;
	m_eq_order = EQUATION_ORDER::NORMAL_ORDER;
}

//-----------------------------------------------------------------------------
FESolver::~FESolver()
{
}

//-----------------------------------------------------------------------------
void FESolver::SetEquationScheme(int scheme)
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
void FESolver::Reset()
{
	m_niter = 0;
	m_nref = 0;
	m_naug = 0;
}

//-----------------------------------------------------------------------------
// get the linear solver
LinearSolver* FESolver::GetLinearSolver()
{
	return nullptr;
}

//-----------------------------------------------------------------------------
//! Matrix symmetry flag
int FESolver::MatrixSymmetryFlag() const
{ 
	return m_msymm; 
}

//-----------------------------------------------------------------------------
//! get matrix type
Matrix_Type FESolver::MatrixType() const
{
	Matrix_Type mtype;
	switch (m_msymm)
	{
	case REAL_UNSYMMETRIC   : mtype = REAL_UNSYMMETRIC; break;
	case REAL_SYMMETRIC     : mtype = REAL_SYMMETRIC; break;
	case REAL_SYMM_STRUCTURE: mtype = REAL_SYMM_STRUCTURE; break;
	}
	return mtype;
}

//-----------------------------------------------------------------------------
// extract the (square) norm of a solution vector
double FESolver::ExtractSolutionNorm(const vector<double>& v, const FEDofList& dofs) const
{
	assert(v.size() == m_dofMap.size());
	double norm = 0;
	for (int n = 0; n < dofs.Size(); ++n)
	{
		for (int i = 0; i < v.size(); ++i)
		{
			if (m_dofMap[i] == dofs[n]) norm += v[i] * v[i];
		}
	}
	return norm;
}

//-----------------------------------------------------------------------------
// return the solution vector
std::vector<double> FESolver::GetSolutionVector() const
{
	return std::vector<double>();
}

//-----------------------------------------------------------------------------
// see if the dofs in the dof list are active in this solver
bool FESolver::HasActiveDofs(const FEDofList& dof)
{
	assert(dof.IsEmpty() == false);
	if (dof.IsEmpty()) return true;

	assert(m_Var.size());
	for (int i = 0; i < dof.Size(); ++i)
	{
		int dof_i = dof[i];

		for (int i = 0; i < m_Var.size(); ++i)
		{
			FESolutionVariable& vi = m_Var[i];
			if (vi.m_dofs->Contains(dof_i))
			{
				return true;
			}
		}
	}
	return false;
}

//-----------------------------------------------------------------------------
// get the active dof map (returns nr of functions)
int FESolver::GetActiveDofMap(vector<int>& activeDofMap)
{
	// get the dof map
	int neq = (int)m_dofMap.size();
	if (m_dofMap.empty() || (m_dofMap.size() < neq)) return -1;

	// We need the partitions here, but for now we assume that
	// it is the first partition

	// The dof map indices point to the dofs as defined by the variables.
	// Since there could be more dofs than actually used in the linear system
	// we need to reindex this map. 
	// First, find the min and max
	int imin = m_dofMap[0], imax = m_dofMap[0];
	for (size_t i = 0; i < neq; ++i)
	{
		if (m_dofMap[i] > imax) imax = m_dofMap[i];
		if (m_dofMap[i] < imin) imin = m_dofMap[i];
	}

	// create the conversion table
	int nsize = imax - imin + 1;
	vector<int> LUT(nsize, -1);
	for (size_t i = 0; i < neq; ++i)
	{
		LUT[m_dofMap[i] - imin] = 1;
	}

	// count how many dofs are actually used
	int nfunc = 0;
	for (size_t i = 0; i < nsize; ++i)
	{
		if (LUT[i] != -1) LUT[i] = nfunc++;
	}

	// now, reindex the dof map
	// allocate dof map
	activeDofMap.resize(neq);
	for (size_t i = 0; i < neq; ++i)
	{
		activeDofMap[i] = LUT[m_dofMap[i] - imin];
	}

	return nfunc;
}

//-----------------------------------------------------------------------------
//! build the matrix profile
void FESolver::BuildMatrixProfile(FEGlobalMatrix& G, bool breset)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	DOFS& fedofs = fem.GetDOFS();
	int MAX_NDOFS = fedofs.GetTotalDOFS();

	// when reset is true we build the entire matrix profile
	// (otherwise we only build the "dynamic" profile)
	if (breset)
	{
		vector<int> elm;

		// Add all elements to the profile
		// Loop over all active domains
		for (int nd = 0; nd<mesh.Domains(); ++nd)
		{
			FEDomain& d = mesh.Domain(nd);
			d.BuildMatrixProfile(G);
		}

		// linear constraints
		FELinearConstraintManager& LCM = fem.GetLinearConstraintManager();
		LCM.BuildMatrixProfile(G);
	}
	else
	{
		// Do the "dynamic" profile. That is the part of the profile that always changes
		// This is mostly contact
		// do the nonlinear constraints
		int M = fem.NonlinearConstraints();
		for (int m = 0; m<M; ++m)
		{
			FENLConstraint* pnlc = fem.NonlinearConstraint(m);
			if (pnlc->IsActive()) pnlc->BuildMatrixProfile(G);
		}

		// All following "elements" are nonstatic. That is, they can change
		// connectivity between calls to this function. All of these elements
		// are related to contact analysis (at this point).
		if (fem.SurfacePairConstraints() > 0)
		{
			// Add all contact interface elements
			for (int i = 0; i<fem.SurfacePairConstraints(); ++i)
			{
				FESurfacePairConstraint* pci = fem.SurfacePairConstraint(i);
				if (pci->IsActive()) pci->BuildMatrixProfile(G);
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! This function is called right before SolveStep and should be used to initialize
//! time dependent information and other settings.
bool FESolver::InitStep(double time)
{
	FEModel& fem = *GetFEModel();

	// evaluate load controllers values at current time
	fem.EvaluateLoadControllers(time);

	// evaluate data generators at current time
	fem.EvaluateDataGenerators(time);

	// evaluate load parameters
	fem.EvaluateLoadParameters();

	// re-validate materials
	// This is necessary since the material parameters can have changed (e.g. via load curves) and thus 
	// a new validation needs to be done to see if the material parameters are still valid. 
	if (fem.ValidateMaterials() == false) return false;

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
bool FESolver::InitEquations()
{
   // get the mesh
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
    
    // clear partitions
	m_part.clear();

	// reorder the node numbers
	int NN = mesh.Nodes();
	vector<int> P(NN);
    
    // see if we need to optimize the bandwidth
	if (m_bwopt)
	{
		FENodeReorder mod;
		mod.Apply(mesh, P);
	}
	else for (int i = 0; i < NN; ++i) P[i] = i;

	for (int i = 0; i < mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(P[i]);
		if (node.HasFlags(FENode::EXCLUDE))
			for (int j = 0; j < (int)node.m_ID.size(); ++j) node.m_ID[j] = -1;
	}
	m_dofMap.clear();

	// assign equations based on allocation scheme
	int neq = 0;
	if (m_eq_scheme == EQUATION_SCHEME::STAGGERED)
	{
		DOFS& dofs = fem.GetDOFS();
		if (m_eq_order == EQUATION_ORDER::NORMAL_ORDER)
		{
			for (int i = 0; i < mesh.Nodes(); ++i)
			{
				FENode& node = mesh.Node(P[i]);
				if (node.HasFlags(FENode::EXCLUDE) == false) {
					for (int nv = 0; nv < dofs.Variables(); ++nv)
					{
						int n = dofs.GetVariableSize(nv);
						for (int l = 0; l < n; ++l)
						{
							int nl = dofs.GetDOF(nv, l);
							if (node.is_active(nl))
							{
								int bcj = node.get_bc(nl);
								if      (bcj == DOF_OPEN      ) { node.m_ID[nl] = neq++; m_dofMap.push_back(nl); }
								else if (bcj == DOF_FIXED     ) { node.m_ID[nl] = -1; }
								else if (bcj == DOF_PRESCRIBED) { node.m_ID[nl] = -neq - 2; neq++; m_dofMap.push_back(nl); }
								else { assert(false); return false; }
							}
							else node.m_ID[nl] = -1;
						}
					}
				}
			}
		}
		else
		{
			int NN = mesh.Nodes();
			for (int i = NN-1; i >= 0; --i)
			{
				FENode& node = mesh.Node(P[i]);
				if (node.HasFlags(FENode::EXCLUDE) == false) {
					int dofs = (int)node.m_ID.size();
					for (int j = dofs - 1; j >= 0; --j)
					{
						if (node.is_active(j))
						{
							int bcj = node.get_bc(j);
							if      (bcj == DOF_OPEN      ) { node.m_ID[j] = neq++; m_dofMap.push_back(j); }
							else if (bcj == DOF_FIXED     ) { node.m_ID[j] = -1; }
							else if (bcj == DOF_PRESCRIBED) { node.m_ID[j] = -neq - 2; neq++; m_dofMap.push_back(j); }
							else { assert(false); return false; }
						}
						else node.m_ID[j] = -1;
					}
				}
			}
		}

		// assign partition
		m_part.push_back(neq);
	}
	else
	{
		// Assign equations numbers in blocks
		assert(m_eq_scheme == EQUATION_SCHEME::BLOCK);
		DOFS& dofs = fem.GetDOFS();

		if (m_eq_order == EQUATION_ORDER::NORMAL_ORDER)
		{
			for (int nv = 0; nv < dofs.Variables(); ++nv)
			{
				int neq0 = neq;
				for (int i = 0; i < NN; ++i)
				{
					FENode& node = mesh.Node(P[i]);
					if (node.HasFlags(FENode::EXCLUDE) == false) {
						int n = dofs.GetVariableSize(nv);
						for (int l = 0; l < n; ++l)
						{
							int nl = dofs.GetDOF(nv, l);

							if (node.is_active(nl))
							{
								int bcl = node.get_bc(nl);
								if      (bcl == DOF_FIXED) { node.m_ID[nl] = -1; }
								else if (bcl == DOF_OPEN) { node.m_ID[nl] = neq++; m_dofMap.push_back(nl); }
								else if (bcl == DOF_PRESCRIBED) { node.m_ID[nl] = -neq - 2; neq++; m_dofMap.push_back(nl); }
								else { assert(false); return false; }
							}
							else node.m_ID[nl] = -1;
						}
					}
				}

				// assign partitions
				if (neq - neq0 > 0)
					m_part.push_back(neq - neq0);
			}
		}
		else if (m_eq_order == EQUATION_ORDER::REVERSE_ORDER)
		{
			int vars = dofs.Variables();
			for (int nv = vars-1; nv >= 0; --nv)
			{
				for (int i = 0; i <NN; ++i)
				{
					FENode& node = mesh.Node(P[i]);
					if (node.HasFlags(FENode::EXCLUDE) == false) {
						int n = dofs.GetVariableSize(nv);
						for (int l = 0; l < n; ++l)
						{
							int nl = dofs.GetDOF(nv, l);
							if (node.is_active(nl))
							{
								int bcl = node.get_bc(nl);
								if      (bcl == DOF_FIXED     ) { node.m_ID[nl] = -1; }
								else if (bcl == DOF_OPEN      ) { node.m_ID[nl] = neq++; m_dofMap.push_back(nl); }
								else if (bcl == DOF_PRESCRIBED) { node.m_ID[nl] = -neq - 2; neq++; m_dofMap.push_back(nl); }
								else { assert(false); return false; }
							}
							else node.m_ID[nl] = -1;
						}
					}
				}

				// assign partitions
				if (nv == vars-1) m_part.push_back(neq);
				else m_part.push_back(neq - m_part[(vars-1) - nv - 1]);
			}
		}
		else
		{
			assert(m_eq_order == FEBIO2_ORDER);

			// Assign equations numbers in blocks
			for (int nv = 0; nv < dofs.Variables(); ++nv)
			{
				int n = dofs.GetVariableSize(nv);
				int neq0 = neq;
				for (int l = 0; l < n; ++l)
				{
					int nl = dofs.GetDOF(nv, l);

					for (int i = 0; i < mesh.Nodes(); ++i)
					{
						FENode& node = mesh.Node(i);
						if (node.is_active(nl))
						{
							int bcl = node.get_bc(nl);
							if      (bcl == DOF_FIXED     ) { node.m_ID[nl] = -1; }
							else if (bcl == DOF_OPEN      ) { node.m_ID[nl] = neq++; m_dofMap.push_back(nl); }
							else if (bcl == DOF_PRESCRIBED) { node.m_ID[nl] = -neq - 2; neq++; m_dofMap.push_back(nl); }
							else { assert(false); return false; }
						}
						else node.m_ID[nl] = -1;
					}
				}

				// assign partitions
				if (neq - neq0 > 0)
					m_part.push_back(neq - neq0);
			}
		}
	}
    
    // store the number of equations
    m_neq = neq;

	assert(m_dofMap.size() == m_neq);
    
    // All initialization is done
    return true;
}

//-----------------------------------------------------------------------------
bool FESolver::InitEquations2()
{
	// get the mesh
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// clear partitions
	m_part.clear();

	// reorder the node numbers
	int NN = mesh.Nodes();
	vector<int> P(NN);

	// see if we need to optimize the bandwidth
	if (m_bwopt)
	{
		FENodeReorder mod;
		mod.Apply(mesh, P);
	}
	else for (int i = 0; i < NN; ++i) P[i] = i;

	// reset all equation numbers
	// first, on all nodes
	for (int i = 0; i < mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(P[i]);
		if (node.HasFlags(FENode::EXCLUDE))
			for (int j = 0; j < (int)node.m_ID.size(); ++j) node.m_ID[j] = -1;
	}
	// then, on all elements
	for (int i = 0; i < mesh.Domains(); ++i)
	{
		FEDomain& dom = mesh.Domain(i);
		for (int j = 0; j < dom.Elements(); ++j)
		{
			FEElement& el = dom.ElementRef(j);
			el.m_lm = -1;
		}
	}
	m_dofMap.clear();

	// see if we need to deactivate some nodal dofs based on requested interpolation order
	for (int i = 0; i < mesh.Domains(); ++i)
	{
		FEDomain& dom = mesh.Domain(i);

		// get the interpolation orders for the different variables.
		for (int n = 0; n < m_Var.size(); ++n)
		{
			FESolutionVariable& var = m_Var[n];
			FEDofList& dofs = *var.m_dofs;
			int P = var.m_order;

			// set unused dofs to -1
			for (int j = 0; j < dom.Elements(); ++j)
			{
				FEElement& el = dom.ElementRef(j);
				int ne = el.Nodes();
				int ne_p = el.ShapeFunctions(P);

				for (int n = ne_p; n < ne; ++n)
				{
					FENode& node = mesh.Node(el.m_node[n]);
					for (int k=0; k<dofs.Size(); ++k)
						node.set_bc(dofs[k], DOF_FIXED);
				}
			}
		}
	}

	// assign equations based on allocation scheme
	int neq = 0;
	if ((m_eq_scheme == EQUATION_SCHEME::STAGGERED) && (m_eq_order == EQUATION_ORDER::NORMAL_ORDER))
	{
		for (int i = 0; i < mesh.Nodes(); ++i)
		{
			FENode& node = mesh.Node(P[i]);
			if (node.HasFlags(FENode::EXCLUDE) == false) 
			{
				int nvar = (int)m_Var.size();
				for (int j = 0; j < nvar; ++j)
				{
					FESolutionVariable& var = m_Var[j];
					FEDofList& dofs = *var.m_dofs;
					for (int k = 0; k < dofs.Size(); ++k)
					{
						int nk = dofs[k];
						if (node.is_active(nk))
						{
							int bck = node.get_bc(nk);
							if      (bck == DOF_OPEN      ) { node.m_ID[nk] = neq++; m_dofMap.push_back(nk); }
							else if (bck == DOF_FIXED     ) { node.m_ID[nk] = -1; }
							else if (bck == DOF_PRESCRIBED) { node.m_ID[nk] = -neq - 2; neq++; m_dofMap.push_back(nk); }
						}
					}
				}
			}
		}

		// assign element dofs
		for (int i = 0; i < mesh.Domains(); ++i)
		{
			FEDomain& dom = mesh.Domain(i);
			for (int j = 0; j < dom.Elements(); ++j)
			{
				FEElement& el = dom.ElementRef(j);
				for (int n = 0; n < m_Var.size(); ++n)
				{
					FESolutionVariable& var = m_Var[n];
					if (var.m_order == 0)
					{
						FEDofList& dofs = *var.m_dofs;
						assert(dofs.Size() == 1);
						assert(el.m_lm == -1);
						el.m_lm = neq++;
						m_dofMap.push_back(dofs[0]);
					}
				}
			}
		}

		// only one partition for this allocation scheme
		m_part.push_back(neq);
	}
	else if ((m_eq_scheme == EQUATION_SCHEME::BLOCK) && (m_eq_order == EQUATION_ORDER::NORMAL_ORDER))
	{
		int neq0 = 0;
		int nvar = (int)m_Var.size();
		for (int j = 0; j < nvar; ++j)
		{
			neq0 = neq;
			FESolutionVariable& var = m_Var[j];
			FEDofList& dofs = *var.m_dofs;
			if (var.m_order != 0)
			{
				for (int i = 0; i < mesh.Nodes(); ++i)
				{
					FENode& node = mesh.Node(P[i]);
					if (node.HasFlags(FENode::EXCLUDE) == false)
					{
						for (int k = 0; k < dofs.Size(); ++k)
						{
							int nk = dofs[k];
							if (node.is_active(nk))
							{
								int bck = node.get_bc(nk);
								if      (bck == DOF_OPEN      ) { node.m_ID[nk] = neq++; m_dofMap.push_back(nk); }
								else if (bck == DOF_FIXED     ) { node.m_ID[nk] = -1; }
								else if (bck == DOF_PRESCRIBED) { node.m_ID[nk] = -neq - 2; neq++; m_dofMap.push_back(nk); }
							}
						}
					}
				}
			}
			else
			{
				assert(dofs.Size() == 1);
				// assign element dofs
				for (int i = 0; i < mesh.Domains(); ++i)
				{
					FEDomain& dom = mesh.Domain(i);
					for (int j = 0; j < dom.Elements(); ++j)
					{
						FEElement& el = dom.ElementRef(j);
						FEDofList& dofs = *var.m_dofs;
						assert(dofs.Size() == 1);
						assert(el.m_lm == -1);
						el.m_lm = neq++;
						m_dofMap.push_back(dofs[0]);
					}
				}
			}

			// add a partition
			int partitionSize = neq - neq0;
			if (partitionSize > 0) m_part.push_back(partitionSize);
		}
	}
	else assert(false);

	// store the number of equations
	m_neq = neq;

	assert(m_dofMap.size() == m_neq);

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
	ar & m_nrhs & m_niter & m_nref & m_ntotref & m_naug;
}

//-----------------------------------------------------------------------------
//! Update the state of the model
void FESolver::Update(std::vector<double>& u)
{ 
	assert(false); 
};

//-----------------------------------------------------------------------------
// The augmentation is done after a time step converges and gives model components
// an opportunity to modify the model's state. This will usually require that the time 
// step is solved again.
bool FESolver::Augment()
{ 
	FEModel& fem = *GetFEModel();

	const FETimeInfo& tp = fem.GetTime();

	// Assume we will pass (can't hurt to be optimistic)
	bool bconv = true;

	// Do contact augmentations
	for (int i = 0; i<fem.SurfacePairConstraints(); ++i)
	{
		FESurfacePairConstraint* pci = fem.SurfacePairConstraint(i);
		if (pci->IsActive()) bconv = (pci->Augment(m_naug, tp) && bconv);
	}

	// do nonlinear constraint augmentations
	for (int i = 0; i<fem.NonlinearConstraints(); ++i)
	{
		FENLConstraint* plc = fem.NonlinearConstraint(i);
		if (plc->IsActive()) bconv = plc->Augment(m_naug, tp) && bconv;
	}

	// do domain augmentations
	FEMesh& mesh = fem.GetMesh();
	for (int i = 0; i<mesh.Domains(); ++i)
	{
		FEDomain& dom = mesh.Domain(i);
		bconv = dom.Augment(m_naug) && bconv;
	}

	fem.GetTime().augmentation++;

	return bconv;
}

//-----------------------------------------------------------------------------
// return the node (mesh index) from an equation number
FENodalDofInfo FESolver::GetDOFInfoFromEquation(int ieq)
{
	FENodalDofInfo info;
	info.m_eq = ieq;
	info.m_node = -1;
	info.m_dof = -1;
	info.szdof = "";

	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	for (int i = 0; i < mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
		vector<int>& id = node.m_ID;
		for (int j = 0; j < id.size(); ++j)
		{
			if (id[j] == ieq)
			{
				info.m_node = node.GetID();
				info.m_dof = j;
				DOFS& Dofs = GetFEModel()->GetDOFS();
				info.szdof = Dofs.GetDOFName(info.m_dof);
				if (info.szdof == nullptr) info.szdof = "???";
				return info;
			}
		}
	}
	return info;
}
