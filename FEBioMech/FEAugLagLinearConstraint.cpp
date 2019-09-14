/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, The Trustees of Columbia University in 
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
#include "FEAugLagLinearConstraint.h"
#include "FECore/FEModel.h"
#include "FECore/log.h"
#include <FECore/FELinearSystem.h>

//-----------------------------------------------------------------------------
void FEAugLagLinearConstraint::Serialize(DumpStream& ar)
{
	if (ar.IsSaving())
	{
		ar << m_lam;

		ar << (int) m_dof.size();
		list<DOF>::iterator it = m_dof.begin();
		for (int i=0; i<(int) m_dof.size(); ++i, ++it) ar << it->bc << it->node << it->val;
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
			ar >> dof.bc >> dof.node >> dof.val;
			m_dof.push_back(dof);
		}
	}
}

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FELinearConstraintSet, FESurfaceConstraint)
	ADD_PARAMETER(m_laugon , "laugon");
	ADD_PARAMETER(m_tol    , "tol");
	ADD_PARAMETER(m_eps    , "penalty");
    ADD_PARAMETER(m_rhs    , "rhs");
	ADD_PARAMETER(m_naugmin, "minaug");
	ADD_PARAMETER(m_naugmax, "maxaug");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FELinearConstraintSet::FELinearConstraintSet(FEModel* pfem) : FESurfaceConstraint(pfem)
{
	static int nc = 1;
	m_nID = nc++;

	m_laugon = true;
	m_eps = 1;
	m_tol = 0.1;
    m_rhs = 0;
	m_naugmax = 50;
	m_naugmin = 0;
}

//-----------------------------------------------------------------------------
void FELinearConstraintSet::BuildMatrixProfile(FEGlobalMatrix& M)
{
	FEMesh& mesh = GetFEModel()->GetMesh();
	list<FEAugLagLinearConstraint*>& LC = m_LC;
	vector<int> lm;
	int N = (int)LC.size();
	list<FEAugLagLinearConstraint*>::iterator it = LC.begin();
	for (int i=0; i<N; ++i, ++it)
	{
		int n = (int)(*it)->m_dof.size();
		lm.resize(n);
		FEAugLagLinearConstraint::Iterator is = (*it)->m_dof.begin();
		for (int j = 0; j<n; ++j, ++is) lm[j] = mesh.Node(is->node).m_ID[is->bc];;
		M.build_add(lm);
	}
}

//-----------------------------------------------------------------------------
//! This function calculates the current value of the constraint.

double FELinearConstraintSet::constraint(FEAugLagLinearConstraint& LC)
{
	int n = (int)LC.m_dof.size();
	double c = 0;
	list<FEAugLagLinearConstraint::DOF>::iterator it = LC.m_dof.begin();
	double u;
	FEMesh& mesh = GetFEModel()->GetMesh();
	for (int i=0; i<n; ++i, ++it) 
	{
		FENode& node = mesh.Node(it->node);
		switch (it->bc)
		{
		case 0: u = node.m_rt.x - node.m_r0.x; break;
		case 1: u = node.m_rt.y - node.m_r0.y; break;
		case 2: u = node.m_rt.z - node.m_r0.z; break;
		default:
                u = node.get(it->bc);
		}
		c += it->val*u;
	}

	return c;
}

//-----------------------------------------------------------------------------
//! This function performs an augmentation, if the Lagrange multiplier 
//! has not converged

bool FELinearConstraintSet::Augment(int naug, const FETimeInfo& tp)
{	
	if (m_laugon == false) return true;

	int M = (int)m_LC.size(), i;
	list<FEAugLagLinearConstraint*>::iterator im = m_LC.begin();

	// calculate lag multipliers
	double L0 = 0, L1 = 0;
	for (i=0; i<M; ++i, ++im)
	{
		FEAugLagLinearConstraint& LC = *(*im);
		double c = constraint(LC) - m_rhs;
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

	feLog("linear constraint set %d: %15.7lg %15.7lg %15.7lg\n", m_nID, L0, fabs(L1 - L0), fabs(m_tol*L1));

	bool bconv = false;
	if (p <= m_tol) bconv = true;
	if ((m_naugmax >= 0) && (naug >= m_naugmax)) bconv = true;
	if (naug < m_naugmin) bconv = false;

	if (bconv == false)
	{
		im = m_LC.begin();
		for (i=0; i<M; ++i, ++im)
		{
			FEAugLagLinearConstraint& LC = *(*im);
			double c = constraint(LC) - m_rhs;
			LC.m_lam += m_eps*c;
		}
	}

	return bconv;
}

//-----------------------------------------------------------------------------
//! This function calculates the contribution to the residual.

void FELinearConstraintSet::LoadVector(FEGlobalVector& R, const FETimeInfo& tp)
{
	FEMesh& mesh = GetFEModel()->GetMesh();

	int M = (int)m_LC.size();
	list<FEAugLagLinearConstraint*>::iterator  im = m_LC.begin();
	for (int m=0; m<M; ++m, ++im)
	{
		FEAugLagLinearConstraint& LC = *(*im);
		int n = (int)LC.m_dof.size();
		double c = constraint(LC);
		FEAugLagLinearConstraint::Iterator it = LC.m_dof.begin();
		for (int i=0; i<n; ++i, ++it)
		{
			int neq = mesh.Node(it->node).m_ID[it->bc];
			if (neq >= 0)
			{
				R[neq] -= (LC.m_lam+m_eps*c)*it->val;
			}		
		}
	}
}

//-----------------------------------------------------------------------------
//! This function calculates the contribution to the stiffness matrix.

void FELinearConstraintSet::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp)
{
	FEMesh& mesh = GetFEModel()->GetMesh();

	vector<int> en;
	vector<int> elm;
	FEElementMatrix ke;

	int M = (int)m_LC.size();
	list<FEAugLagLinearConstraint*>::iterator im = m_LC.begin();
	for (int m=0; m<M; ++m, ++im)
	{
		FEAugLagLinearConstraint& LC = *(*im);
		int n = (int)LC.m_dof.size(), i, j;
		ke.resize(n, n);
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
			int neq = mesh.Node(it->node).m_ID[it->bc];
			elm[i] = neq;
		}
		ke.SetNodes(en);
		ke.SetIndices(elm);
		LS.Assemble(ke);
	}
}

//-----------------------------------------------------------------------------
void FELinearConstraintSet::Serialize(DumpStream& ar)
{
	if (ar.IsShallow())
	{
		if (ar.IsSaving())
		{
			list<FEAugLagLinearConstraint*>::iterator it = m_LC.begin();
			for (int i=0; i<(int) m_LC.size(); ++i, ++it) ar << (*it)->m_lam;
		}
		else
		{
			list<FEAugLagLinearConstraint*>::iterator it = m_LC.begin();
			for (int i=0; i<(int) m_LC.size(); ++i, ++it) ar >> (*it)->m_lam;
		}
	}
	else
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
}
