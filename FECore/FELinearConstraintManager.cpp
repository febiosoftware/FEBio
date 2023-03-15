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
#include "FELinearConstraintManager.h"
#include "FEMesh.h"
#include "FEModel.h"
#include "FEAnalysis.h"
#include "DumpStream.h"
#include "FEDomain.h"
#include "FEGlobalMatrix.h"

//-----------------------------------------------------------------------------
FELinearConstraintManager::FELinearConstraintManager(FEModel* fem) : m_fem(fem)
{
}

//-----------------------------------------------------------------------------
void FELinearConstraintManager::Clear()
{
	for (size_t i = 0; i < m_LinC.size(); ++i) delete m_LinC[i];
	m_LinC.clear();
}

//-----------------------------------------------------------------------------
FELinearConstraintManager::~FELinearConstraintManager()
{
	Clear();
}

//-----------------------------------------------------------------------------
void FELinearConstraintManager::CopyFrom(const FELinearConstraintManager& lcm)
{
	Clear();
	for (int i=0; i<lcm.LinearConstraints(); ++i)
	{
		FELinearConstraint* lc = fecore_alloc(FELinearConstraint, m_fem);
		lc->CopyFrom(&(const_cast<FELinearConstraint&>(lcm.LinearConstraint(i))));
		m_LinC.push_back(lc);
	}
	InitTable();
}

//-----------------------------------------------------------------------------
void FELinearConstraintManager::AddLinearConstraint(FELinearConstraint* lc)
{
	m_LinC.push_back(lc);
}

//-----------------------------------------------------------------------------
int FELinearConstraintManager::LinearConstraints() const
{
	return (int)m_LinC.size();
}

//-----------------------------------------------------------------------------
const FELinearConstraint& FELinearConstraintManager::LinearConstraint(int i) const
{
	return *m_LinC[i];
}

//-----------------------------------------------------------------------------
FELinearConstraint& FELinearConstraintManager::LinearConstraint(int i)
{
	return *m_LinC[i];
}

//-----------------------------------------------------------------------------
//! remove a linear constraint
void FELinearConstraintManager::RemoveLinearConstraint(int i)
{
	FELinearConstraint& lc = *m_LinC[i];
	if (lc.IsActive()) lc.Deactivate();
	m_LinC.erase(m_LinC.begin() + i);
}

//-----------------------------------------------------------------------------
void FELinearConstraintManager::Serialize(DumpStream& ar)
{
	if (ar.IsShallow()) return;

	if (ar.IsSaving())
	{
		ar << m_LinC;

		int nr = m_LCT.rows();
		int nc = m_LCT.columns();
		ar << nr << nc;
		if (nr*nc > 0)
		{
			ar.write(&m_LCT(0,0), sizeof(int), nr*nc);
		}
	}
	else
	{
		FEModel& fem = ar.GetFEModel();

		// linear constraints
		ar >> m_LinC;

		int nr, nc;
		ar >> nr >> nc;
		if (nr*nc > 0)
		{
			m_LCT.resize(nr, nc);
			ar.read(&m_LCT(0,0), sizeof(int), nr*nc);
		}
	}
}


//-----------------------------------------------------------------------------
void FELinearConstraintManager::BuildMatrixProfile(FEGlobalMatrix& G)
{
	int nlin = (int)m_LinC.size();
	if (nlin == 0) return;

	FEAnalysis* pstep = m_fem->GetCurrentStep();
	FEMesh& mesh = m_fem->GetMesh();

	// Add linear constraints to the profile
	// TODO: we need to add a function build_add(lmi, lmj) for
	// this type of "elements". Now we allocate too much memory
	vector<int> lm, elm;

	// do the cross-term
	// TODO: I have to make this easier. For instance,
	// keep a list that stores for each node the list of
	// elements connected to that node.
	// loop over all solid elements
	for (int nd = 0; nd<pstep->Domains(); ++nd)
	{
		FEDomain& dom = *pstep->Domain(nd);
		for (int i = 0; i<dom.Elements(); ++i)
		{
			FEElement& el = dom.ElementRef(i);
			dom.UnpackLM(el, elm);
			int ne = (int)elm.size();

			// keep a list of the constraints this element connects to
			vector<int> constraintList;

			// see if this element connects to the 
			// parent node of a linear constraint ...
			int m = el.Nodes();
			for (int j = 0; j<m; ++j)
			{	
				int ncols = m_LCT.columns();
				for (int k = 0; k<ncols; ++k)
				{
					int n = m_LCT(el.m_node[j], k);
					if (n >= 0)
					{
						// ... it does so we need to connect the 
						// element to the linear constraint
						FELinearConstraint* plc = m_LinC[n];
						constraintList.push_back(n);
						
						int ns = (int)plc->Size();

						lm.resize(ne + ns);
						for (int l = 0; l<ne; ++l) lm[l] = elm[l];

						FELinearConstraint::dof_iterator is = plc->begin();
						for (int l = ne; l<ne + ns; ++l, ++is) 
						{
							int neq = mesh.Node((*is)->node).m_ID[(*is)->dof];
							lm[l] = neq;
						}

						G.build_add(lm);
					}
				}
			}

			// This replaces the commented out section below, which sets up the connectivity 
			// for the constraint block of the stiffness matrix. 
			// The problem is that this will only work for linear constraints that are connected
			// the element, so not for constraints connected to contact "elements". Howerver, I
			// don't think this was working anyway, so this is okay for now.
			// TODO: Replace this with a more generic algorithm that looks at the matrix profile 
			// directly instead of element by element. 
			if (constraintList.empty() == false)
			{
				// do the constraint term
				int n = 0;
				for (int i = 0; i<constraintList.size(); ++i) n += (int)m_LinC[constraintList[i]]->Size();
				lm.resize(n);
				n = 0;
				for (int i = 0; i<constraintList.size(); ++i)
				{
					FELinearConstraint& lc = *m_LinC[constraintList[i]];
					int ni = (int)lc.Size();
					for (int j = 0; j<ni; ++j)
					{
						const FELinearConstraintDOF& sj = lc.GetChildDof(j);
						int neq = mesh.Node(sj.node).m_ID[sj.dof];
						lm[n++] = neq;
					}
				}
				G.build_add(lm);
			}
		}
	}

	// do the constraint term
	// NOTE: This block was replaced by the section above which reduces the size
	// of the stiffness matrix, but might be less generic (altough not entirely sure
	// about that). 
/*	vector<FELinearConstraint>::iterator ic = m_LinC.begin();
	int n = 0;
	for (int i = 0; i<nlin; ++i, ++ic) n += (int) ic->m_childDof.size();
	lm.resize(n);
	ic = m_LinC.begin();
	n = 0;
	for (int i = 0; i<nlin; ++i, ++ic)
	{
		int ni = (int)ic->m_childDof.size();
		vector<FELinearConstraintDOF>::iterator is = ic->m_childDof.begin();
		for (int j = 0; j<ni; ++j, ++is) 
		{
			int neq = mesh.Node(is->node).m_ID[is->dof];
			lm[n++] = neq;
		}
	}
	G.build_add(lm);
*/
}

//-----------------------------------------------------------------------------
//! This function initializes the linear constraint table (LCT). This table
//! contains for each dof the linear constraint it belongs to. (or -1 if it is
//! not constraint)
bool FELinearConstraintManager::Initialize()
{
	return true;
}

void FELinearConstraintManager::PrepStep()
{
	FEMesh& mesh = m_fem->GetMesh();

	for (int i=0; i<m_LinC.size(); ++i)
	{
		FELinearConstraint& lc = *m_LinC[i];
		if (lc.IsActive())
		{
			FENode& node = mesh.Node(lc.GetParentNode());
			double u = node.get(lc.GetParentDof());

			double v = 0;
			for (int j = 0; j < lc.Size(); ++j)
			{
				const FELinearConstraintDOF& dofj = lc.GetChildDof(j);
				FENode& nj = mesh.Node(dofj.node);
				v += dofj.val * nj.get(dofj.dof);
			}

			m_up[i] = v + lc.GetOffset() - u;
		}
	}
}

void FELinearConstraintManager::InitTable()
{
	FEMesh& mesh = m_fem->GetMesh();

	// create the linear constraint table
	DOFS& fedofs = m_fem->GetDOFS();
	int MAX_NDOFS = fedofs.GetTotalDOFS();
	m_LCT.resize(mesh.Nodes(), MAX_NDOFS, -1);
	m_LCT.set(-1);

	vector<FELinearConstraint*>::iterator ic = m_LinC.begin();
	int nlin = LinearConstraints();
	for (int i = 0; i<nlin; ++i, ++ic)
	{
		FELinearConstraint& lc = *(*ic);
		if (lc.IsActive())
		{
			int n = lc.GetParentNode();
			int m = lc.GetParentDof();

			m_LCT(n, m) = i;
		}
	}
}

//-----------------------------------------------------------------------------
// This gets called during model activation, i.e. activation of permanent model components.
bool FELinearConstraintManager::Activate()
{
	int nlin = LinearConstraints();
	if (nlin == 0) return true;

	// initialize the lookup table
	InitTable();

	// ensure that none of the parent nodes are child nodes in any of the active linear constraints
	for (int i = 0; i<nlin; ++i)
	{
		FELinearConstraint& lci = *m_LinC[i];
		if (lci.IsActive())
		{
			int n = (int)lci.Size();
			for (int k = 0; k<n; ++k)
			{
				const FELinearConstraintDOF& childDOF = lci.GetChildDof(k);
				int n = m_LCT(childDOF.node, childDOF.dof);
				if (n != -1)
				{
					return false;
				}
			}
		}
	}

	// set the prescribed value array
	m_up.assign(m_LinC.size(), 0.0);
	if (m_LinC.size())
	{
		vector<FELinearConstraint*>::iterator il = m_LinC.begin();
		for (int l = 0; l<(int)m_LinC.size(); ++l, ++il) 
		{
			if ((*il)->IsActive()) (*il)->Activate();
		}
	}

	return true;
}

//-----------------------------------------------------------------------------
void FELinearConstraintManager::AssembleResidual(vector<double>& R, vector<int>& en, vector<int>& elm, vector<double>& fe)
{
	FEMesh& mesh = m_fem->GetMesh();

	int ndof = (int)fe.size();
	int ndn = ndof / (int)en.size();
	const int nodes = (int)en.size();

	// loop over all degrees of freedom of this element
	for (int i = 0; i<ndof; ++i)
	{
		int nodei = i / ndn;
		if (nodei < nodes) {
			// see if this dof belongs to a linear constraint
			int l = m_LCT(en[nodei], i%ndn);
			if (l >= 0)
			{
				// if so, get the linear constraint
				FELinearConstraint& lc = *m_LinC[l];
				assert(elm[i] == -1);

				// now loop over all child nodes and
				// add the contribution to the residual
				int ns = (int)lc.Size();
				FELinearConstraint::dof_iterator is = lc.begin();
				for (int j = 0; j < ns; ++j, ++is)
				{
					int I = mesh.Node((*is)->node).m_ID[(*is)->dof];
					if (I >= 0)
					{
						double A = (*is)->val;
#pragma omp atomic
						R[I] += A*fe[i];
					}
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------
void FELinearConstraintManager::AssembleStiffness(FEGlobalMatrix& G, vector<double>& R, vector<double>& ui, const vector<int>& en, const vector<int>& lmi, const vector<int>& lmj, const matrix& ke)
{
	FEMesh& mesh = m_fem->GetMesh();

	int ndof = ke.rows();
	int ndn = ndof / (int)en.size();
	const int nodes = (int)en.size();

	SparseMatrix& K = *(&G);

	// loop over all stiffness components 
	// and correct for linear constraints
	for (int i = 0; i<ndof; ++i)
	{
		int nodei = i / ndn;
		int li = (nodei < nodes ? m_LCT(en[nodei], i%ndn) : -1);
		for (int j = 0; j < ndof; ++j)
		{
			int nodej = j / ndn;
			int lj = (nodej < nodes ? m_LCT(en[nodej], j%ndn) : -1);
			if ((li >= 0) && (lj < 0))
			{
				// dof i is constrained
				FELinearConstraint& Li = *m_LinC[li];

				assert(lmi[i] == -1);

				FELinearConstraint::dof_iterator is = Li.begin();
				for (int k = 0; k < (int)Li.Size(); ++k, ++is)
				{
					int I = mesh.Node((*is)->node).m_ID[(*is)->dof];
					int J = lmj[j];
					double kij = (*is)->val*ke[i][j];
					if ((J >= 0) && (I >= 0)) K.add(I, J, kij);
					else
					{
						// adjust for prescribed dofs
						J = -J - 2;
						if ((J >= 0) && (I >= 0)) R[I] -= kij*ui[J];
					}
				}
			}
			else if ((lj >= 0) && (li < 0))
			{
				// dof j is constrained
				FELinearConstraint& Lj = *m_LinC[lj];

				assert(lmj[j] == -1);

				FELinearConstraint::dof_iterator js = Lj.begin();
				for (int k = 0; k < (int)Lj.Size(); ++k, ++js)
				{
					int I = lmi[i];
					int J = mesh.Node((*js)->node).m_ID[(*js)->dof];
					double kij = (*js)->val*ke[i][j];
					if ((J >= 0) && (I >= 0)) K.add(I, J, kij);
					else
					{
						// adjust for prescribed dofs
						J = -J - 2;
						if ((J >= 0) && (I >= 0)) R[I] -= kij*ui[J];
					}
				}

				// adjust right-hand side for inhomogeneous linear constraints
				if (Lj.GetOffset() != 0.0)
				{
					double ri = ke[i][j] * m_up[lj];
					int I = lmi[i];
					if (I >= 0) R[I] -= ri;
				}
			}
			else if ((li >= 0) && (lj >= 0))
			{
				// both dof i and j are constrained
				FELinearConstraint& Li = *m_LinC[li];
				FELinearConstraint& Lj = *m_LinC[lj];

				FELinearConstraint::dof_iterator is = Li.begin();
				FELinearConstraint::dof_iterator js = Lj.begin();

				assert(lmi[i] == -1);
				assert(lmj[j] == -1);

				for (int k = 0; k < (int)Li.Size(); ++k, ++is)
				{
					js = Lj.begin();
					for (int l = 0; l < (int)Lj.Size(); ++l, ++js)
					{
						int I = mesh.Node((*is)->node).m_ID[(*is)->dof];
						int J = mesh.Node((*js)->node).m_ID[(*js)->dof];;
						double kij = ke[i][j] * (*is)->val*(*js)->val;

						if ((J >= 0) && (I >= 0)) K.add(I, J, kij);
						else
						{
							// adjust for prescribed dofs
							J = -J - 2;
							if ((J >= 0) && (I >= 0)) R[I] -= kij*ui[J];
						}
					}
				}

				// adjust for inhomogeneous linear constraints
				if (Lj.GetOffset() != 0.0)
				{
					is = Li.begin();
					for (int k = 0; k < (int)Li.Size(); ++k, ++is)
					{
						int I = mesh.Node((*is)->node).m_ID[(*is)->dof];
						double ri = (*is)->val * ke[i][j] * m_up[lj];
						if (I >= 0) R[I] -= ri;
					}
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------
// This updates the nodal degrees of freedom of the parent nodes.
void FELinearConstraintManager::Update()
{
	FEMesh& mesh = m_fem->GetMesh();

	int nlin = LinearConstraints();
	for (int n = 0; n<nlin; ++n)
	{
		FELinearConstraint& lc = LinearConstraint(n);
		if (lc.IsActive())
		{
			// evaluate the linear constraint
			double d = 0;
			int ns = (int)lc.Size();
			FELinearConstraint::dof_iterator si = lc.begin();
			for (int i = 0; i < ns; ++i, ++si)
			{
				FENode& childNode = mesh.Node((*si)->node);
				d += (*si)->val * childNode.get((*si)->dof);
			}

			// assign to parent node
			FENode& parentNode = mesh.Node(lc.GetParentNode());
			parentNode.set(lc.GetParentDof(), d + lc.GetOffset());
		}
	}

	m_up.assign(m_LinC.size(), 0.0);
}
