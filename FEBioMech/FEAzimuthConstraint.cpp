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
#include "FEAzimuthConstraint.h"
#include "FEBioMech.h"
#include <FECore/FELinearSystem.h>
#include <FECore/FENode.h>
#include <FECore/log.h>

BEGIN_FECORE_CLASS(FEAzimuthConstraint, FENodeSetConstraint)
	ADD_PARAMETER(m_laugon, "laugon");
	ADD_PARAMETER(m_tol, "augtol");
	ADD_PARAMETER(m_eps, "penalty");
	ADD_PARAMETER(m_minaug, "minaug");
	ADD_PARAMETER(m_maxaug, "maxaug");
	END_FECORE_CLASS();

FEAzimuthConstraint::FEAzimuthConstraint(FEModel* fem) : FENodeSetConstraint(fem), m_dofU(fem), m_nodeSet(fem)
{
	m_laugon = 0;
	m_tol = 0.1;
	m_eps = 0.0;

	m_minaug = 0;
	m_maxaug = 0;

	// TODO: Can this be done in Init, since there is no error checking
	if (fem)
	{
		m_dofU.AddVariable(FEBioMech::GetVariableName(FEBioMech::DISPLACEMENT));
	}
}

//-----------------------------------------------------------------------------
bool FEAzimuthConstraint::Init()
{
	if (FENLConstraint::Init() == false) return false;

	// make sure there is a node set
	if (m_nodeSet.Size() == 0) return false;

	return true;
}

//-----------------------------------------------------------------------------
void FEAzimuthConstraint::Serialize(DumpStream& ar)
{
	FENLConstraint::Serialize(ar);
	ar & m_Lm;
}

//-----------------------------------------------------------------------------
//! return node set
FENodeSet* FEAzimuthConstraint::GetNodeSet()
{
	return &m_nodeSet;
}

//-----------------------------------------------------------------------------
void FEAzimuthConstraint::BuildMatrixProfile(FEGlobalMatrix& M)
{
	int N = m_nodeSet.Size();
	for (int i = 0; i < N; ++i)
	{
		vector<int> lm(2);
		FENode& n0 = *m_nodeSet.Node(i);
		lm[0] = n0.m_ID[m_dofU[0]];
		lm[1] = n0.m_ID[m_dofU[1]];

		if (m_laugon == 2)
		{
			lm.push_back(m_eq[2*i  ]);
			lm.push_back(m_eq[2*i+1]);
		}

		M.build_add(lm);
	}
}

//-----------------------------------------------------------------------------
//! return number of equations to be allocated for Lagrange multipliers
int FEAzimuthConstraint::InitEquations(int neq)
{
	int N = m_nodeSet.Size();
	m_Lm.assign(N, vec3d(0, 0, 0));

	// make sure we want to use Lagrange Multiplier method
	if (m_laugon != 2) return 0;

	// allocate two equations per node
	m_eq.resize(2 * N);
	for (int i = 0; i < N; ++i)
	{
		m_eq[2*i  ] = neq++;
		m_eq[2*i+1] = neq++;
	}
	return 2*N;
}

//! Calculate the constraint force
void FEAzimuthConstraint::LoadVector(FEGlobalVector& R, const FETimeInfo& tp)
{
	int N = m_nodeSet.Size();
	for (int i = 0; i < N; ++i)
	{
		FENode& node = *m_nodeSet.Node(i);

		vec3d e0 = node.m_r0; e0.z = 0.0; e0.unit();
		vec3d et = node.m_rt; et.z = 0.0; double l = et.unit();

		// projection mat3ds
		mat3dd I(1.0);
		mat3ds P = (I - dyad(et)) / l;

		// evaluate constraint 
		vec3d c = et - e0;

		// calculate force
		vector<double> F;
		vector<int> lm;
		if (m_laugon == 2)
		{
			// Lagrange multiplier formulation
			vec3d Pl = P * m_Lm[i];
			F.resize(4);
			F[0] = -Pl.x;
			F[1] = -Pl.y;
			F[2] = -c.x;
			F[3] = -c.y;

			lm.resize(4);
			lm[0] = node.m_ID[m_dofU[0]];
			lm[1] = node.m_ID[m_dofU[1]];
			lm[2] = m_eq[2*i  ];
			lm[3] = m_eq[2*i+1];
		}
		else
		{
			// Augmented Lagrangian formulation
			vec3d lam = m_Lm[i] + c * m_eps;
			vec3d Pl = P * lam;

			F.resize(2);
			F[0] = -Pl.x;
			F[1] = -Pl.y;

			lm.resize(2);
			lm[0] = node.m_ID[m_dofU[0]];
			lm[1] = node.m_ID[m_dofU[1]];
		}

		R.Assemble(lm, F);
	}
}

//! calculate the constraint stiffness
void FEAzimuthConstraint::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp)
{
	int N = m_nodeSet.Size();
	for (int i = 0; i < N; ++i)
	{
		FENode* node = m_nodeSet.Node(i); assert(node);

		vec3d e0 = node->m_r0; e0.z = 0.0; e0.unit();
		vec3d et = node->m_rt; et.z = 0.0; double l = et.unit();

		// projection mat3ds
		mat3dd I(1.0);
		mat3ds P = (I - dyad(et)) / l;

		// evaluate constraint 
		vec3d c = et - e0;

		mat3ds Ploe = dyads(P*m_Lm[i], et);
		mat3ds Q = (Ploe + P * (et*m_Lm[i])) / (-l);

		// stiffness
		if (m_laugon == 2)
		{
			FEElementMatrix ke(4, 4); ke.zero();
			ke[0][0] = Q(0, 0); ke[0][1] = Q(0, 1);
			ke[1][0] = Q(1, 0); ke[1][1] = Q(1, 1);

			ke[0][2] = P(0, 0); ke[0][3] = P(0, 1);
			ke[1][2] = P(1, 0); ke[1][3] = P(1, 1);

			ke[2][0] = P(0, 0); ke[3][0] = P(0, 1);
			ke[2][1] = P(1, 0); ke[3][1] = P(1, 1);

			vector<int> lm(4);
			lm[0] = node->m_ID[m_dofU[0]];
			lm[1] = node->m_ID[m_dofU[1]];
			lm[2] = m_eq[2*i  ];
			lm[3] = m_eq[2*i+1];
			ke.SetIndices(lm);

			vector<int> en(1);
			en[0] = m_nodeSet[i];
			ke.SetNodes(en);

			LS.Assemble(ke);
		}
		else
		{
			vec3d lam = m_Lm[i] + c * m_eps;

			mat3ds K = Q + (P.sqr())*m_eps;

			FEElementMatrix ke(2, 2); ke.zero();
			ke[0][0] = K(0, 0); ke[0][1] = K(0, 1);
			ke[1][0] = K(1, 0); ke[1][1] = K(1, 1);

			vector<int> lm(2);
			lm[0] = node->m_ID[m_dofU[0]];
			lm[1] = node->m_ID[m_dofU[1]];

			vector<int> en(1);
			en[0] = m_nodeSet[i];
			ke.SetIndices(lm);
			ke.SetNodes(en);
			LS.Assemble(ke);
		}
	}
}

void FEAzimuthConstraint::Update(const std::vector<double>& ui)
{
	if (m_laugon == 2)
	{
		int N = m_nodeSet.Size();
		for (int i = 0; i < N; ++i)
		{
			m_Lm[i].x += ui[m_eq[2*i  ]];
			m_Lm[i].y += ui[m_eq[2*i+1]];
		}
	}
}

//! augmentations
bool FEAzimuthConstraint::Augment(int naug, const FETimeInfo& tp)
{
	bool bconv = true;
	int N = m_nodeSet.Size();
	for (int i = 0; i < N; ++i)
	{
		FENode& node = *m_nodeSet.Node(i);
		int nodeId = node.GetID();

		vec3d e0 = node.m_r0; e0.z = 0.0; e0.unit();
		vec3d et = node.m_rt; et.z = 0.0; double l = et.unit();

		// projection mat3ds
		mat3dd I(1.0);
		mat3ds P = (I - dyad(et)) / l;

		// Lagrange multiplier
		vec3d c = et - e0;
		vec3d lam1 = m_Lm[i] + c * m_eps;

		double g = c.norm();
		double L0 = m_Lm[i].norm();
		double L1 = lam1.norm();

		if (m_laugon == 0)
		{
			// penalty-formulation
			feLog("\t%d: %lg, %lg\n", nodeId, g, L1);
		}
		else if (m_laugon == 2)
		{
			// Lagrange multiplier
			feLog("\t%d: %lg, %lg\n", nodeId, g, L0);
		}
		else
		{
			// aug-lag constraint
			bool bconvi = true;
			if (m_Lm[i].norm() == 0.0)
			{
				feLog("\t%d: %lg, %lg\n", nodeId, g, L1);
				bconvi = (L1 == 0.0);
			}
			else
			{
				double dl = fabs((L1 - L0) / L1);
				feLog("\t%d: %lg, %lg (%lg)\n", nodeId, g, L1, dl);
				if (dl >= m_tol)
				{
					bconvi = false;
				}
			}

			// check augmentation range
			naug++;
			if ((m_minaug > 0) && (naug < m_minaug)) bconvi = false;
			if ((m_maxaug > 0) && (naug >= m_maxaug)) bconvi = true;

			// update LM if not converged
			if (bconvi == false)
			{
				m_Lm[i] = lam1;
			}

			bconv = bconv & bconvi;
		}
	}

	return bconv;
}
