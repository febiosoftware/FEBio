#include "stdafx.h"
#include "FELinearConstraintManager.h"
#include "FEMesh.h"
#include "FEModel.h"
#include "FEAnalysis.h"
#include "DumpStream.h"

//-----------------------------------------------------------------------------
FELinearConstraintManager::FELinearConstraintManager(FEModel* fem) : m_fem(fem)
{
}

//-----------------------------------------------------------------------------
void FELinearConstraintManager::AddLinearConstraint(FELinearConstraint& lc)
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
	return m_LinC[i];
}

//-----------------------------------------------------------------------------
void FELinearConstraintManager::Serialize(DumpStream& ar)
{
	if (ar.IsShallow()) return;

	if (ar.IsSaving())
	{
		ar << (int)m_LinC.size();
		vector<FELinearConstraint>::iterator it = m_LinC.begin();
		for (int i = 0; i<(int)m_LinC.size(); ++i, ++it) it->Serialize(ar);

		ar << m_LCT;
	}
	else
	{
		FEModel& fem = ar.GetFEModel();
		// linear constraints
		int n = 0;
		ar >> n;
		FELinearConstraint LC(&fem);
		for (int i = 0; i<n; ++i)
		{
			LC.Serialize(ar);
			m_LinC.push_back(LC);
		}

		ar >> m_LCT;
	}
}


//-----------------------------------------------------------------------------
void FELinearConstraintManager::BuildMatrixProfile(FEGlobalMatrix& G)
{
	FEAnalysis* pstep = m_fem->GetCurrentStep();

	DOFS& fedofs = m_fem->GetDOFS();
	int MAX_NDOFS = fedofs.GetTotalDOFS();

	// Add linear constraints to the profile
	// TODO: we need to add a function build_add(lmi, lmj) for
	// this type of "elements". Now we allocate too much memory
	if (m_LinC.size() > 0)
	{
		int nlin = (int)m_LinC.size();
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

				// see if this element connects to the 
				// master node of a linear constraint ...
				int m = el.Nodes();
				for (int j = 0; j<m; ++j)
				{
					for (int k = 0; k<MAX_NDOFS; ++k)
					{
						int n = m_LCT[el.m_node[j] * MAX_NDOFS + k];
						if (n >= 0)
						{
							// ... it does so we need to connect the 
							// element to the linear constraint
							FELinearConstraint* plc = &m_LinC[n];

							int ns = (int)plc->slave.size();

							lm.resize(ne + ns);
							for (int l = 0; l<ne; ++l) lm[l] = elm[l];

							list<FELinearConstraint::SlaveDOF>::iterator is = plc->slave.begin();
							for (int l = ne; l<ne + ns; ++l, ++is) lm[l] = is->neq;

							G.build_add(lm);
							break;
						}
					}
				}
			}
		}

		// TODO: do the same thing for shell elements

		// do the constraint term
		int ni;
		vector<FELinearConstraint>::iterator ic = m_LinC.begin();
		int n = 0;
		for (int i = 0; i<nlin; ++i, ++ic) n += (int) ic->slave.size();
		lm.resize(n);
		ic = m_LinC.begin();
		n = 0;
		for (int i = 0; i<nlin; ++i, ++ic)
		{
			ni = (int)ic->slave.size();
			list<FELinearConstraint::SlaveDOF>::iterator is = ic->slave.begin();
			for (int j = 0; j<ni; ++j, ++is) lm[n++] = is->neq;
		}
		G.build_add(lm);
	}
}

//-----------------------------------------------------------------------------
//! This function initializes the linear constraint table (LCT). This table
//! contains for each dof the linear constraint it belongs to. (or -1 if it is
//! not constraint)
bool FELinearConstraintManager::Initialize()
{
	int nlin = LinearConstraints();
	if (nlin == 0) return true;
	int i;

	FEMesh& mesh = m_fem->GetMesh();

	// set the equation numbers for the linear constraints
	vector<FELinearConstraint>::iterator it = m_LinC.begin();
	for (i = 0; i<nlin; ++i, ++it)
	{
		FELinearConstraint& lc = *it;
		lc.master.neq = mesh.Node(lc.master.node).m_ID[lc.master.bc];

		// make sure the master did not get assigned an equation
		assert(lc.master.neq == -1);

		// set the slave equation numbers
		list<FELinearConstraint::SlaveDOF>::iterator is = lc.slave.begin();
		int nn = (int) lc.slave.size();
		for (int n = 0; n<nn; ++n, ++is)
		{
			FELinearConstraint::SlaveDOF& sn = *is;
			sn.neq = mesh.Node(sn.node).m_ID[sn.bc];
		}
	}

	// create the linear constraint table
	DOFS& fedofs = m_fem->GetDOFS();
	int MAX_NDOFS = fedofs.GetTotalDOFS();
	m_LCT.assign(mesh.Nodes()*MAX_NDOFS, -1);

	vector<FELinearConstraint>::iterator ic = m_LinC.begin();
	for (i = 0; i<nlin; ++i, ++ic)
	{
		FELinearConstraint& lc = *ic;
		int n = lc.master.node;
		int m = lc.master.bc;

		m_LCT[n*MAX_NDOFS + m] = i;
	}

	return true;
}

//-----------------------------------------------------------------------------
void FELinearConstraintManager::Activate()
{
	if (m_LinC.size())
	{
		vector<FELinearConstraint>::iterator il = m_LinC.begin();
		for (int l = 0; l<(int)m_LinC.size(); ++l, ++il) il->Activate();
	}
}

//-----------------------------------------------------------------------------
void FELinearConstraintManager::AssembleResidual(vector<double>& R, vector<int>& en, vector<int>& elm, vector<double>& fe)
{
	// get nodal DOFS
	DOFS& fedofs = m_fem->GetDOFS();
	int MAX_NDOFS = fedofs.GetTotalDOFS();

	int ndof = (int)fe.size();
	int ndn = ndof / (int)en.size();

	// loop over all degrees of freedom of this element
	for (int i = 0; i<ndof; ++i)
	{
		// see if this dof belongs to a linear constraint
		int n = MAX_NDOFS*(en[i / ndn]) + i%ndn;
		int l = m_LCT[n];

		if (l >= 0)
		{
			// if so, get the linear constraint
			FELinearConstraint& lc = m_LinC[l];
			assert(elm[i] == -1);

			// now loop over all "slave" nodes and
			// add the contribution to the residual
			int ns = (int)lc.slave.size();
			list<FELinearConstraint::SlaveDOF>::iterator is = lc.slave.begin();
			for (int j = 0; j<ns; ++j, ++is)
			{
				int I = is->neq;
				if (I >= 0)
				{
					double A = is->val;
					R[I] += A*fe[i];
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------
void FELinearConstraintManager::AssembleStiffness(FEGlobalMatrix& G, vector<double>& R, vector<double>& ui, vector<int>& en, vector<int>& elm, matrix& ke)
{
	// get nodal DOFS
	DOFS& fedofs = m_fem->GetDOFS();
	int MAX_NDOFS = fedofs.GetTotalDOFS();

	int i, j, l;

	int ndof = ke.rows();
	int ndn = ndof / (int)en.size();

	SparseMatrix& K = *(&G);

	// loop over all stiffness components 
	// and correct for linear constraints
	int ni, nj, li, lj, I, J, k;
	double kij;
	for (i = 0; i<ndof; ++i)
	{
		ni = MAX_NDOFS*(en[i / ndn]) + i%ndn;
		li = m_LCT[ni];
		for (j = 0; j<ndof; ++j)
		{
			nj = MAX_NDOFS*(en[j / ndn]) + j%ndn;
			lj = m_LCT[nj];

			if ((li >= 0) && (lj < 0))
			{
				// dof i is constrained
				FELinearConstraint& Li = m_LinC[li];

				assert(elm[i] == -1);

				list<FELinearConstraint::SlaveDOF>::iterator is = Li.slave.begin();
				for (k = 0; k<(int)Li.slave.size(); ++k, ++is)
				{
					I = is->neq;
					J = elm[j];
					kij = is->val*ke[i][j];
					if ((J >= I) && (I >= 0)) K.add(I, J, kij);
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
				FELinearConstraint& Lj = m_LinC[lj];

				assert(elm[j] == -1);

				list<FELinearConstraint::SlaveDOF>::iterator js = Lj.slave.begin();

				for (k = 0; k<(int)Lj.slave.size(); ++k, ++js)
				{
					I = elm[i];
					J = js->neq;
					kij = js->val*ke[i][j];
					if ((J >= I) && (I >= 0)) K.add(I, J, kij);
					else
					{
						// adjust for prescribed dofs
						J = -J - 2;
						if ((J >= 0) && (I >= 0)) R[I] -= kij*ui[J];
					}
				}
			}
			else if ((li >= 0) && (lj >= 0))
			{
				// both dof i and j are constrained
				FELinearConstraint& Li = m_LinC[li];
				FELinearConstraint& Lj = m_LinC[lj];

				list<FELinearConstraint::SlaveDOF>::iterator is = Li.slave.begin();
				list<FELinearConstraint::SlaveDOF>::iterator js = Lj.slave.begin();

				assert(elm[i] == -1);
				assert(elm[j] == -1);

				for (k = 0; k<(int)Li.slave.size(); ++k, ++is)
				{
					js = Lj.slave.begin();
					for (l = 0; l<(int)Lj.slave.size(); ++l, ++js)
					{
						I = is->neq;
						J = js->neq;
						kij = ke[i][j] * is->val*js->val;

						if ((J >= I) && (I >= 0)) K.add(I, J, kij);
						else
						{
							// adjust for prescribed dofs
							J = -J - 2;
							if ((J >= 0) && (I >= 0)) R[I] -= kij*ui[J];
						}
					}
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------
void FELinearConstraintManager::Update()
{
	FEMesh& mesh = m_fem->GetMesh();

	int nlin = LinearConstraints();
	for (int n = 0; n<nlin; ++n)
	{
		const FELinearConstraint& lc = LinearConstraint(n);

		// evaluate the linear constraint
		double d = 0;
		int ns = (int)lc.slave.size();
		list<FELinearConstraint::SlaveDOF>::const_iterator si = lc.slave.begin();
		for (int i = 0; i<ns; ++i, ++si)
		{
			FENode& slaveNode = mesh.Node(si->node);
			d += si->val*slaveNode.get(si->bc);
		}

		// assign to master node
		FENode& masterNode = mesh.Node(lc.master.node);
		masterNode.set(lc.master.bc, d);
	}
}
