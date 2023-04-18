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
#include "FEMortarTiedContact.h"
#include "FECore/mortar.h"
#include "FECore/FEGlobalMatrix.h"
#include "FECore/log.h"
#include <FECore/FELinearSystem.h>

//=============================================================================
// FEMortarTiedSurface
//=============================================================================

//-----------------------------------------------------------------------------
FEMortarTiedSurface::FEMortarTiedSurface(FEModel* pfem) : FEMortarContactSurface(pfem) {}

//-----------------------------------------------------------------------------
bool FEMortarTiedSurface::Init()
{
	// always intialize base class first!
	if (FEMortarContactSurface::Init() == false) return false;

	// get the number of nodes
	int NN = Nodes();

	// allocate data structures
	m_L.resize(NN, vec3d(0,0,0));

	return true;
}

//=============================================================================
// FEMortarTiedContact
//=============================================================================

//-----------------------------------------------------------------------------
// Define sliding interface parameters
BEGIN_FECORE_CLASS(FEMortarTiedContact, FEMortarInterface)
	ADD_PARAMETER(m_atol   , "tolerance"    );
	ADD_PARAMETER(m_eps    , "penalty"      );
	ADD_PARAMETER(m_naugmin, "minaug"       );
	ADD_PARAMETER(m_naugmax, "maxaug"       );
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEMortarTiedContact::FEMortarTiedContact(FEModel* pfem) : FEMortarInterface(pfem), m_ss(pfem), m_ms(pfem)
{
	m_dofX = GetDOFIndex("x");
	m_dofY = GetDOFIndex("y");
	m_dofZ = GetDOFIndex("z");
}

//-----------------------------------------------------------------------------
bool FEMortarTiedContact::Init()
{
	// initialize surfaces
	if (m_ms.Init() == false) return false;
	if (m_ss.Init() == false) return false;

	return true;
}

//-----------------------------------------------------------------------------
void FEMortarTiedContact::Activate()
{
	//! don't forget the base class
	FEContactInterface::Activate();

	m_ss.UpdateNodalAreas();

	// update the mortar weights
	// For tied interfaces, this is only done once, during activation
	UpdateMortarWeights(m_ss, m_ms);

	// update the nodal gaps
	// (must be done after mortar eights are updated)
	UpdateNodalGaps(m_ss, m_ms);
}

//-----------------------------------------------------------------------------
//! build the matrix profile for use in the stiffness matrix
void FEMortarTiedContact::BuildMatrixProfile(FEGlobalMatrix& K)
{
	// For now we'll assume that each node on the primary side is connected to the secondary side
	// This is obviously too much, but we'll worry about improving this later
	int NS = m_ss.Nodes();
	int NM = m_ms.Nodes();
	vector<int> LM(3*(NS+NM));
	for (int i=0; i<NS; ++i)
	{
		FENode& ni = m_ss.Node(i);
		LM[3*i  ] = ni.m_ID[0];
		LM[3*i+1] = ni.m_ID[1];
		LM[3*i+2] = ni.m_ID[2];
	}
	for (int i=0; i<NM; ++i)
	{
		FENode& ni = m_ms.Node(i);
		LM[3*NS + 3*i  ] = ni.m_ID[0];
		LM[3*NS + 3*i+1] = ni.m_ID[1];
		LM[3*NS + 3*i+2] = ni.m_ID[2];
	}
	K.build_add(LM);
}

//-----------------------------------------------------------------------------
//! calculate contact forces
void FEMortarTiedContact::LoadVector(FEGlobalVector& R, const FETimeInfo& tp)
{
	int NS = m_ss.Nodes();
	int NM = m_ms.Nodes();

	// loop over all primary nodes
	for (int A=0; A<NS; ++A)
	{
		double eps = m_eps*m_ss.m_A[A];
		vec3d gA = m_ss.m_gap[A];
		vec3d tA = m_ss.m_L[A] + gA*eps;

		// loop over all primary nodes
		vector<int> en(1);
		vector<int> lm(3);
		vector<double> fe(3);
		for (int B=0; B<NS; ++B)
		{
			FENode& nodeB = m_ss.Node(B);
			en[0] = m_ss.NodeIndex(B);
			lm[0] = nodeB.m_ID[m_dofX];
			lm[1] = nodeB.m_ID[m_dofY];
			lm[2] = nodeB.m_ID[m_dofZ];

			double nAB = -m_n1[A][B];
			if (nAB != 0.0)
			{
				fe[0] = tA.x*nAB;
				fe[1] = tA.y*nAB;
				fe[2] = tA.z*nAB;

				R.Assemble(en, lm, fe);
			}
		}

		// loop over secondary side
		for (int C=0; C<NM; ++C)
		{
			FENode& nodeC = m_ms.Node(C);
			en[0] = m_ms.NodeIndex(C);
			lm[0] = nodeC.m_ID[m_dofX];
			lm[1] = nodeC.m_ID[m_dofY];
			lm[2] = nodeC.m_ID[m_dofZ];

			double nAC = m_n2[A][C];
			if (nAC != 0.0)
			{
				fe[0] = tA.x*nAC;
				fe[1] = tA.y*nAC;
				fe[2] = tA.z*nAC;

				R.Assemble(en, lm, fe);
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! calculate contact stiffness
void FEMortarTiedContact::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp)
{
	int NS = m_ss.Nodes();
	int NM = m_ms.Nodes();

	// A. Linearization of the gap function
	vector<int> lmi(3), lmj(3);
	FEElementMatrix ke;
	ke.resize(3, 3);
	for (int A=0; A<NS; ++A)
	{
		double eps = m_eps*m_ss.m_A[A];

		// loop over all primary nodes
		for (int B=0; B<NS; ++B)
		{
			FENode& nodeB = m_ss.Node(B);
			lmi[0] = nodeB.m_ID[0];
			lmi[1] = nodeB.m_ID[1];
			lmi[2] = nodeB.m_ID[2];

			double nAB = m_n1[A][B]*eps;
			if (nAB != 0.0)
			{
				// loop over primary nodes
				for (int C=0; C<NS; ++C)
				{
					FENode& nodeC = m_ss.Node(C);
					lmj[0] = nodeC.m_ID[0];
					lmj[1] = nodeC.m_ID[1];
					lmj[2] = nodeC.m_ID[2];

					double nAC = m_n1[A][C]*nAB;
					if (nAC != 0.0)
					{
						ke[0][0] = nAC; ke[0][1] = 0.0; ke[0][2] = 0.0;
						ke[1][0] = 0.0; ke[1][1] = nAC; ke[1][2] = 0.0;
						ke[2][0] = 0.0; ke[2][1] = 0.0; ke[2][2] = nAC;

						ke.SetIndices(lmi, lmj);
						LS.Assemble(ke);
					}
				}

				// loop over secondary nodes
				for (int C=0; C<NM; ++C)
				{
					FENode& nodeC = m_ms.Node(C);
					lmj[0] = nodeC.m_ID[0];
					lmj[1] = nodeC.m_ID[1];
					lmj[2] = nodeC.m_ID[2];

					double nAC = -m_n2[A][C]*nAB;
					if (nAC != 0.0)
					{
						ke[0][0] = nAC; ke[0][1] = 0.0; ke[0][2] = 0.0;
						ke[1][0] = 0.0; ke[1][1] = nAC; ke[1][2] = 0.0;
						ke[2][0] = 0.0; ke[2][1] = 0.0; ke[2][2] = nAC;

						ke.SetIndices(lmi, lmj);
						LS.Assemble(ke);
					}
				}
			}
		}

		// loop over all secondary nodes
		for (int B=0; B<NM; ++B)
		{
			FENode& nodeB = m_ms.Node(B);
			lmi[0] = nodeB.m_ID[0];
			lmi[1] = nodeB.m_ID[1];
			lmi[2] = nodeB.m_ID[2];

			double nAB = -m_n2[A][B]*eps;
			if (nAB != 0.0)
			{
				// loop over primary nodes
				for (int C=0; C<NS; ++C)
				{
					FENode& nodeC = m_ss.Node(C);
					lmj[0] = nodeC.m_ID[0];
					lmj[1] = nodeC.m_ID[1];
					lmj[2] = nodeC.m_ID[2];

					double nAC = m_n1[A][C]*nAB;
					if (nAC != 0.0)
					{
						ke[0][0] = nAC; ke[0][1] = 0.0; ke[0][2] = 0.0;
						ke[1][0] = 0.0; ke[1][1] = nAC; ke[1][2] = 0.0;
						ke[2][0] = 0.0; ke[2][1] = 0.0; ke[2][2] = nAC;

						ke.SetIndices(lmi, lmj);
						LS.Assemble(ke);
					}
				}

				// loop over secondary nodes
				for (int C=0; C<NM; ++C)
				{
					FENode& nodeC = m_ms.Node(C);
					lmj[0] = nodeC.m_ID[0];
					lmj[1] = nodeC.m_ID[1];
					lmj[2] = nodeC.m_ID[2];

					double nAC = -m_n2[A][C]*nAB;
					if (nAC != 0.0)
					{
						ke[0][0] = nAC; ke[0][1] = 0.0; ke[0][2] = 0.0;
						ke[1][0] = 0.0; ke[1][1] = nAC; ke[1][2] = 0.0;
						ke[2][0] = 0.0; ke[2][1] = 0.0; ke[2][2] = nAC;

						ke.SetIndices(lmi, lmj);
						LS.Assemble(ke);
					}
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! calculate Lagrangian augmentations
bool FEMortarTiedContact::Augment(int naug, const FETimeInfo& tp)
{
	if (m_laugon != 1) return true;

	double max_err = 0.0;
	int NS = m_ss.Nodes();
	// loop over all primary nodes
	for (int A=0; A<NS; ++A)
	{
		double eps = m_eps*m_ss.m_A[A];
		vec3d gA = m_ss.m_gap[A];
		vec3d Lold = m_ss.m_L[A];
		vec3d Lnew = Lold + gA*eps;

		double uold = Lold.norm();
		double unew = Lnew.norm();

		double err = fabs((uold - unew)/(uold + unew));
		if (err > max_err) max_err = err;
	}

	bool bconv = true;
	if ((m_atol > 0) && (max_err > m_atol)) bconv = false;
	if (m_naugmin > naug) bconv = false;
	if (m_naugmax <= naug) bconv = true;

	feLog(" mortar interface # %d\n", GetID());
	feLog("                        CURRENT        REQUIRED\n");
	feLog("    normal force : %15le", max_err);
	if (m_atol > 0) feLog("%15le\n", m_atol); else feLog("       ***\n");

	if (bconv == false)
	{
		// loop over all primary nodes
		for (int A=0; A<NS; ++A)
		{
			double eps = m_eps*m_ss.m_A[A];
			vec3d gA = m_ss.m_gap[A];
			vec3d Lold = m_ss.m_L[A];
			vec3d Lnew = Lold + gA*eps;
			m_ss.m_L[A] = Lnew;
		}
	}

	return bconv;
}

//-----------------------------------------------------------------------------
//! update interface data
void FEMortarTiedContact::Update()
{
	UpdateNodalGaps(m_ss, m_ms);
}

//-----------------------------------------------------------------------------
//! serialize data to archive
void FEMortarTiedContact::Serialize(DumpStream& ar)
{
}
