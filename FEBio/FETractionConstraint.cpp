#include "stdafx.h"
#include "FETractionConstraint.h"
#include "fem.h"
#include "FESolidSolver.h"
#include "log.h"

//-----------------------------------------------------------------------------
void FETractionConstraint::Serialize(Archive& ar)
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
FETractionConstraintSet::FETractionConstraintSet(FEM* pfem)
{
	static int nc = 1;
	m_nID = nc++;

	m_eps = 1;
	m_tol = 0.1;
	m_naugmax = 50;

	// store a pointer to the FEM data
	m_pfem = pfem;
}

void FETractionConstraintSet::Init()
{
	// set the equation numbers for the traction constraints
	list<FETractionConstraint*>::iterator it = m_RC.begin();
	int N = m_RC.size();
	FEMesh& mesh = m_pfem->m_mesh;
	for (int i=0; i<N; ++i, ++it)
	{
		FETractionConstraint& rc = *(*it);

		// set the slave equation numbers
		FETractionConstraint::Iterator is = rc.m_dof.begin();
		int nn = rc.m_dof.size();
		for (int n=0; n<nn; ++n, ++is)
		{
			FETractionConstraint::DOF& sn = *is;
			sn.neq = mesh.Node(sn.node).m_ID[sn.bc];
		}
	}
}

//-----------------------------------------------------------------------------
//! This function calculates the current value of the constraint.

double FETractionConstraintSet::constraint(FETractionConstraint& RC)
{
	FESolidSolver* psolver = dynamic_cast<FESolidSolver*>(m_pfem->m_pStep->m_psolver);
	vector<double>& T = psolver->m_Ti;

	int n = RC.m_dof.size();
	double c = 0;
	list<FETractionConstraint::DOF>::iterator it = RC.m_dof.begin();
	for (int i=0; i<n; ++i, ++it)
	{
		if (it->neq >= 0) c += it->val*T[it->neq];
	}

	return c;
}

//-----------------------------------------------------------------------------
//! This function calculates the contribution to the residual.

void FETractionConstraintSet::Residual(vector<double> &R)
{
	int M = m_RC.size();
	list<FETractionConstraint*>::iterator  im = m_RC.begin();
	for (int m=0; m<M; ++m, ++im)
	{
		FETractionConstraint& RC = *(*im);
		int n = RC.m_dof.size();
		double c = constraint(RC);
		FETractionConstraint::Iterator it = RC.m_dof.begin();
		for (int i=0; i<n; ++i, ++it)
		{
			if (it->neq >= 0)
			{
				R[it->neq] += (RC.m_lam+m_eps*c)*it->val;
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! This function calculates the contribution to the stiffness matrix.

void FETractionConstraintSet::Stiffness()
{

}

//-----------------------------------------------------------------------------
//! This function performs an augmentation, if the Lagrange multiplier
//! has not converged

bool FETractionConstraintSet::Augment(int naug)
{
	int M = m_RC.size(), i;
	list<FETractionConstraint*>::iterator im = m_RC.begin();

	// calculate lag multipliers
	double L0 = 0, L1 = 0;
	for (i=0; i<M; ++i, ++im)
	{
		FETractionConstraint& RC = *(*im);
		double c = constraint(RC);
		double lam = RC.m_lam + m_eps*c;

		L0 += RC.m_lam*RC.m_lam;
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

	log.printf("traction constraint set %d: %15.7lg %15.7lg %15.7lg\n", m_nID, L0, fabs(L1 - L0), fabs(m_tol*L1));

	if ((m_naugmax >= 0) && (naug >= m_naugmax)) return true;

	if (p<= m_tol)
	{
//		log.printf("(conv)\n");
		im = m_RC.begin();
		for (i=0; i<M; ++i, ++im)
		{
			double c = constraint(*(*im));
			log.printf("constraint %d: %lg\n", i+1, c);
		}
		return true;
	}
	else
	{
		im = m_RC.begin();
		for (i=0; i<M; ++i, ++im)
		{
			FETractionConstraint& RC = *(*im);
			double c = constraint(RC);
			RC.m_lam += m_eps*c;
		}
//		log.printf("\n");
		return false;
	}

	return true;
}

//-----------------------------------------------------------------------------

void FETractionConstraintSet::Serialize(Archive& ar)
{
	if (ar.IsSaving())
	{
		ar << m_tol << m_eps << m_naugmax << m_nID;
		ar << (int) m_RC.size();
		list<FETractionConstraint*>::iterator it = m_RC.begin();
		for (int i=0; i<(int) m_RC.size(); ++i, ++it) (*it)->Serialize(ar);
	}
	else
	{
		ar >> m_tol >> m_eps >> m_naugmax >> m_nID;
		int n;
		ar >> n;
		m_RC.clear();
		for (int i=0; i<n; ++i)
		{
			FETractionConstraint* prc = new FETractionConstraint;
			prc->Serialize(ar);
			m_RC.push_back(prc);
		}
	}
}
