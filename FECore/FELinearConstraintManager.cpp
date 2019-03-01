#include "stdafx.h"
#include "FELinearConstraintManager.h"
#include "FEMesh.h"
#include "FEModel.h"
#include "FEAnalysis.h"
#include "DumpStream.h"
#include "FEDomain.h"

//-----------------------------------------------------------------------------
FELinearConstraintManager::FELinearConstraintManager(FEModel* fem) : m_fem(fem)
{
}

//-----------------------------------------------------------------------------
void FELinearConstraintManager::Clear()
{
	m_LinC.clear();
}

//-----------------------------------------------------------------------------
void FELinearConstraintManager::CopyFrom(const FELinearConstraintManager& lcm)
{
	m_LinC.clear();
	for (int i=0; i<lcm.LinearConstraints(); ++i)
	{
		FELinearConstraint lc(m_fem);
		lc.CopyFrom(lcm.LinearConstraint(i));
		m_LinC.push_back(lc);
	}
	InitTable();
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
FELinearConstraint& FELinearConstraintManager::LinearConstraint(int i)
{
	return m_LinC[i];
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
		m_LCT.resize(nr, nc);
		if (nr*nc > 0)
		{
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

			// see if this element connects to the 
			// master node of a linear constraint ...
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
						FELinearConstraint* plc = &m_LinC[n];
						
						int ns = (int)plc->slave.size();

						lm.resize(ne + ns);
						for (int l = 0; l<ne; ++l) lm[l] = elm[l];

						vector<FELinearConstraint::DOF>::iterator is = plc->slave.begin();
						for (int l = ne; l<ne + ns; ++l, ++is) 
						{
							int neq = mesh.Node(is->node).m_ID[is->dof];
							lm[l] = neq;
						}

						G.build_add(lm);
					}
				}
			}
		}
	}

	// do the constraint term
	vector<FELinearConstraint>::iterator ic = m_LinC.begin();
	int n = 0;
	for (int i = 0; i<nlin; ++i, ++ic) n += (int) ic->slave.size();
	lm.resize(n);
	ic = m_LinC.begin();
	n = 0;
	for (int i = 0; i<nlin; ++i, ++ic)
	{
		int ni = (int)ic->slave.size();
		vector<FELinearConstraint::DOF>::iterator is = ic->slave.begin();
		for (int j = 0; j<ni; ++j, ++is) 
		{
			int neq = mesh.Node(is->node).m_ID[is->dof];
			lm[n++] = neq;
		}
	}
	G.build_add(lm);
}

//-----------------------------------------------------------------------------
//! This function initializes the linear constraint table (LCT). This table
//! contains for each dof the linear constraint it belongs to. (or -1 if it is
//! not constraint)
bool FELinearConstraintManager::Initialize()
{
	int nlin = LinearConstraints();
	if (nlin == 0) return true;

	// ensure that none of the master nodes are slave nodes in any of the active linear constraints
	for (int i=0; i<nlin; ++i)
	{
		FELinearConstraint& lci = m_LinC[i];
		if (lci.IsActive())
		{
			FELinearConstraint::DOF& masterDOF = lci.master;
			for (int j=0; j<nlin; j++)
			{
				FELinearConstraint& lcj = m_LinC[j];
				if (lcj.IsActive())
				{
					int n = (int)lcj.slave.size();
					for (int k=0; k<n; ++k)
					{
						FELinearConstraint::DOF& slaveDOF = lcj.slave[k];
						if ((slaveDOF.node == masterDOF.node) && (slaveDOF.dof == masterDOF.dof)) 
						{
							return false;
						}
					}

					// also make sure the master dof is not repeated
					if (i != j)
					{
						if ((lci.master.node == lcj.master.node) && (lci.master.dof == lcj.master.dof)) 
						{
							return false;
						}
					}
				}
			}
		}
	}

	// initialize the lookup table
	InitTable();

	// set the prescribed value array
	m_up.assign(m_LinC.size(), 0.0);

	return true;
}

void FELinearConstraintManager::PrepStep()
{
	FEMesh& mesh = m_fem->GetMesh();

	for (int i=0; i<m_LinC.size(); ++i)
	{
		FELinearConstraint& lc = m_LinC[i];
		FENode& node = mesh.Node(lc.master.node);
		double u = node.get(lc.master.dof);

		double v = 0;
		for (int j=0; j<lc.slave.size(); ++j)
		{
			FENode& nj = mesh.Node(lc.slave[j].node);
			v += lc.slave[j].val* nj.get(lc.slave[j].dof);
		}

		m_up[i] = v + lc.m_off - u;
	}
}

void FELinearConstraintManager::InitTable()
{
	FEMesh& mesh = m_fem->GetMesh();

	// create the linear constraint table
	DOFS& fedofs = m_fem->GetDOFS();
	int MAX_NDOFS = fedofs.GetTotalDOFS();
	m_LCT.resize(mesh.Nodes(), MAX_NDOFS, -1);

	vector<FELinearConstraint>::iterator ic = m_LinC.begin();
	int nlin = LinearConstraints();
	for (int i = 0; i<nlin; ++i, ++ic)
	{
		FELinearConstraint& lc = *ic;
		int n = lc.master.node;
		int m = lc.master.dof;

		m_LCT(n, m) = i;
	}
}

//-----------------------------------------------------------------------------
// This gets called during model activation, i.e. activation of permanent model components.
void FELinearConstraintManager::Activate()
{
	if (m_LinC.size())
	{
		vector<FELinearConstraint>::iterator il = m_LinC.begin();
		for (int l = 0; l<(int)m_LinC.size(); ++l, ++il) 
		{
			if (il->IsActive()) il->Activate();
		}
	}
}

//-----------------------------------------------------------------------------
void FELinearConstraintManager::AssembleResidual(vector<double>& R, vector<int>& en, vector<int>& elm, vector<double>& fe)
{
	FEMesh& mesh = m_fem->GetMesh();

	int ndof = (int)fe.size();
	int ndn = ndof / (int)en.size();

	// loop over all degrees of freedom of this element
	for (int i = 0; i<ndof; ++i)
	{
		// see if this dof belongs to a linear constraint
		int l = m_LCT(en[i / ndn], i%ndn);
		if (l >= 0)
		{
			// if so, get the linear constraint
			FELinearConstraint& lc = m_LinC[l];
			assert(elm[i] == -1);

			// now loop over all "slave" nodes and
			// add the contribution to the residual
			int ns = (int)lc.slave.size();
			vector<FELinearConstraint::DOF>::iterator is = lc.slave.begin();
			for (int j = 0; j<ns; ++j, ++is)
			{
				int I = mesh.Node(is->node).m_ID[is->dof];
				if (I >= 0)
				{
					double A = is->val;
#pragma omp atomic
					R[I] += A*fe[i];
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------
void FELinearConstraintManager::AssembleStiffness(FEGlobalMatrix& G, vector<double>& R, vector<double>& ui, vector<int>& en, vector<int>& elm, matrix& ke)
{
	FEMesh& mesh = m_fem->GetMesh();

	int ndof = ke.rows();
	int ndn = ndof / (int)en.size();

	SparseMatrix& K = *(&G);

	// loop over all stiffness components 
	// and correct for linear constraints
	for (int i = 0; i<ndof; ++i)
	{
		int li = m_LCT(en[i / ndn], i%ndn);
		for (int j = 0; j<ndof; ++j)
		{
			int lj = m_LCT(en[j / ndn], j%ndn);

			if ((li >= 0) && (lj < 0))
			{
				// dof i is constrained
				FELinearConstraint& Li = m_LinC[li];

				assert(elm[i] == -1);

				vector<FELinearConstraint::DOF>::iterator is = Li.slave.begin();
				for (int k = 0; k<(int)Li.slave.size(); ++k, ++is)
				{
					int I = mesh.Node(is->node).m_ID[is->dof];
					int J = elm[j];
					double kij = is->val*ke[i][j];
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
				FELinearConstraint& Lj = m_LinC[lj];

				assert(elm[j] == -1);

				vector<FELinearConstraint::DOF>::iterator js = Lj.slave.begin();

				for (int k = 0; k<(int)Lj.slave.size(); ++k, ++js)
				{
					int I = elm[i];
					int J = mesh.Node(js->node).m_ID[js->dof];
					double kij = js->val*ke[i][j];
					if ((J >= 0) && (I >= 0)) K.add(I, J, kij);
					else
					{
						// adjust for prescribed dofs
						J = -J - 2;
						if ((J >= 0) && (I >= 0)) R[I] -= kij*ui[J];
					}
				}

				// adjust right-hand side for inhomogeneous linear constraints
				if (Lj.m_off != 0.0)
				{
					double ri = ke[i][j]*m_up[lj];
					int I = elm[i];
					if (I >= 0) R[i] -= ri;
				}
			}
			else if ((li >= 0) && (lj >= 0))
			{
				// both dof i and j are constrained
				FELinearConstraint& Li = m_LinC[li];
				FELinearConstraint& Lj = m_LinC[lj];

				vector<FELinearConstraint::DOF>::iterator is = Li.slave.begin();
				vector<FELinearConstraint::DOF>::iterator js = Lj.slave.begin();

				assert(elm[i] == -1);
				assert(elm[j] == -1);

				for (int k = 0; k<(int)Li.slave.size(); ++k, ++is)
				{
					js = Lj.slave.begin();
					for (int l = 0; l<(int)Lj.slave.size(); ++l, ++js)
					{
						int I = mesh.Node(is->node).m_ID[is->dof];
						int J = mesh.Node(js->node).m_ID[js->dof];;
						double kij = ke[i][j] * is->val*js->val;

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
				if (Lj.m_off != 0.0)
				{
					is = Li.slave.begin();
					for (int k = 0; k<(int)Li.slave.size(); ++k, ++is)
					{
						int I = mesh.Node(is->node).m_ID[is->dof];
						double ri = is->val * ke[i][j] * m_up[lj];
						if (I >= 0) R[i] -= ri;
					}
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------
// This updates the nodal degrees of freedom of the master nodes.
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
		vector<FELinearConstraint::DOF>::const_iterator si = lc.slave.begin();
		for (int i = 0; i<ns; ++i, ++si)
		{
			FENode& slaveNode = mesh.Node(si->node);
			d += si->val*slaveNode.get(si->dof);
		}

		// assign to master node
		FENode& masterNode = mesh.Node(lc.master.node);
		masterNode.set(lc.master.dof, d + lc.m_off);
	}

	m_up.assign(m_LinC.size(), 0.0);
}
