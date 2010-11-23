#include "stdafx.h"
#include "FEAugLagLinearConstraint.h"
#include "fem.h"
#include "FESolidSolver.h"
#include "log.h"

//-----------------------------------------------------------------------------
void FEAugLagLinearConstraint::Serialize(DumpFile& ar)
{
	if (ar.IsSaving())
	{
		ar << m_lam;

		ar << (int) m_dof.size();
		list<DOF>::iterator it = m_dof.begin();
		for (int i=0; i<(int) m_dof.size(); ++i, ++it) ar << it->bc << it->neq << it->node << it->val;
	}
	else
	{
		ar >> m_lam;

		int n;
		ar >> n;
		DOF dof;
		m_dof.clear();
		for (int i=0; i<n; ++i)
		{
			ar >> dof.bc >> dof.neq >> dof.node >> dof.val;
			m_dof.push_back(dof);
		}
	}
}

//-----------------------------------------------------------------------------
FELinearConstraintSet::FELinearConstraintSet(FEM* pfem)
{
	static int nc = 1;
	m_nID = nc++;

	m_eps = 1;
	m_tol = 0.1;
	m_naugmax = 50;

	// store a pointer to the FEM data
	m_pfem = pfem;
}

void FELinearConstraintSet::Init()
{
	// set the equation numbers for the linear constraints
	list<FEAugLagLinearConstraint*>::iterator it = m_LC.begin();
	int N = m_LC.size();
	FEMesh& mesh = m_pfem->m_mesh;
	for (int i=0; i<N; ++i, ++it)
	{
		FEAugLagLinearConstraint& lc = *(*it);

		// set the slave equation numbers
		FEAugLagLinearConstraint::Iterator is = lc.m_dof.begin();
		int nn = lc.m_dof.size();
		for (int n=0; n<nn; ++n, ++is)
		{
			FEAugLagLinearConstraint::DOF& sn = *is;
			sn.neq = mesh.Node(sn.node).m_ID[sn.bc];
		}		
	}
}

//-----------------------------------------------------------------------------
//! This function calculates the current value of the constraint.

double FELinearConstraintSet::constraint(FEAugLagLinearConstraint& LC)
{
	int n = LC.m_dof.size();
	double c = 0;
	list<FEAugLagLinearConstraint::DOF>::iterator it = LC.m_dof.begin();
	double u;
	FEMesh* pm = &m_pfem->m_mesh;
	for (int i=0; i<n; ++i, ++it) 
	{
		FENode& node = pm->Node(it->node);
		switch (it->bc)
		{
		case 0: u = node.m_rt.x - node.m_r0.x; break;
		case 1: u = node.m_rt.y - node.m_r0.y; break;
		case 2: u = node.m_rt.z - node.m_r0.z; break;
		default:
			assert(false);
		}
		c += it->val*u;
	}

	return c;
}

//-----------------------------------------------------------------------------
//! This function performs an augmentation, if the Lagrange multiplier 
//! has not converged

bool FELinearConstraintSet::Augment(int naug)
{
	int M = m_LC.size(), i;
	list<FEAugLagLinearConstraint*>::iterator im = m_LC.begin();

	// calculate lag multipliers
	double L0 = 0, L1 = 0;
	for (i=0; i<M; ++i, ++im)
	{
		FEAugLagLinearConstraint& LC = *(*im);
		double c = constraint(LC);
		double lam = LC.m_lam + m_eps*c;

		L0 += LC.m_lam*LC.m_lam;
		L1 += lam*lam;
	}

	L0 = sqrt(L0);
	L1 = sqrt(L1);

	double p;
	if (L1 != 0)
		p = fabs((L1 - L0)/L1);
	else p = fabs(L1 - L0);

	// get the logfile
	Logfile& log = GetLogfile();

	log.printf("linear constraint set %d: %15.7lg %15.7lg %15.7lg", m_nID, L0, fabs(L1 - L0), fabs(m_tol*L1));

	if ((m_naugmax >= 0) && (naug >= m_naugmax)) return true;

	if (p<= m_tol) 
	{
//		log.printf("(conv)\n");
		return true;
	}
	else 
	{
		im = m_LC.begin();
		for (i=0; i<M; ++i, ++im)
		{
			FEAugLagLinearConstraint& LC = *(*im);
			double c = constraint(LC);
			LC.m_lam += m_eps*c;
		}
//		log.printf("\n");
		return false;
	}

	return true;
}

//-----------------------------------------------------------------------------
//! This function calculates the contribution to the residual.

void FELinearConstraintSet::Residual(vector<double> &R)
{
	int M = m_LC.size();
	list<FEAugLagLinearConstraint*>::iterator  im = m_LC.begin();
	for (int m=0; m<M; ++m, ++im)
	{
		FEAugLagLinearConstraint& LC = *(*im);
		int n = LC.m_dof.size();
		double c = constraint(LC);
		FEAugLagLinearConstraint::Iterator it = LC.m_dof.begin();
		for (int i=0; i<n; ++i, ++it)
		{
			if (it->neq >= 0)
			{
				R[it->neq] -= (LC.m_lam+m_eps*c)*it->val;
			}		
		}
	}
}

//-----------------------------------------------------------------------------
//! This function calculates the contribution to the stiffness matrix.

void FELinearConstraintSet::Stiffness()
{
	FESolidSolver* psolver = dynamic_cast<FESolidSolver*>(m_pfem->m_pStep->m_psolver);

	vector<int> en;
	vector<int> elm;
	matrix ke;

	int M = m_LC.size();
	list<FEAugLagLinearConstraint*>::iterator im = m_LC.begin();
	for (int m=0; m<M; ++m, ++im)
	{
		FEAugLagLinearConstraint& LC = *(*im);
		int n = LC.m_dof.size(), i, j;
		ke.Create(n, n);
		FEAugLagLinearConstraint::Iterator it = LC.m_dof.begin(), jt;
		for (i=0; i<n; ++i, ++it)
		{
			jt = LC.m_dof.begin();
			for (j=0; j<n; ++j, ++jt)
			{
				ke[i][j] = m_eps*it->val*jt->val;
			}
		}

		en.resize(n);
		elm.resize(n);
		it = LC.m_dof.begin();
		for (i=0; i<n; ++i, ++it)
		{
			en[i] = it->node;
			elm[i] = it->neq;
		}

		psolver->AssembleStiffness(en, elm, ke);
	}
}

//-----------------------------------------------------------------------------

void FELinearConstraintSet::Serialize(DumpFile& ar)
{
	if (ar.IsSaving())
	{
		ar << m_tol << m_eps << m_naugmax << m_nID;
		ar << (int) m_LC.size();
		list<FEAugLagLinearConstraint*>::iterator it = m_LC.begin();
		for (int i=0; i<(int) m_LC.size(); ++i, ++it) (*it)->Serialize(ar);
	}
	else
	{
		ar >> m_tol >> m_eps >> m_naugmax >> m_nID;
		int n;
		ar >> n;
		m_LC.clear();
		for (int i=0; i<n; ++i)
		{
			FEAugLagLinearConstraint* plc = new FEAugLagLinearConstraint;
			plc->Serialize(ar);
			m_LC.push_back(plc);
		}
	}
}
