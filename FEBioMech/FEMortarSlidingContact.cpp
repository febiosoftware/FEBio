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
#include "FEMortarSlidingContact.h"
#include "FECore/mortar.h"
#include "FECore/FEGlobalMatrix.h"
#include "FECore/log.h"
#include <FECore/FELinearSystem.h>
#include <FECore/FEDataExport.h>

//=============================================================================
// FEMortarSlidingSurface
//=============================================================================

//-----------------------------------------------------------------------------
FEMortarSlidingSurface::FEMortarSlidingSurface(FEModel* pfem) : FEMortarContactSurface(pfem) 
{
	// class exports
	EXPORT_DATA(PLT_VEC3F, FMT_NODE, &m_gap, "mortar-gap vector");
	EXPORT_DATA(PLT_VEC3F, FMT_NODE, &m_nu , "mortar-normal"    );
}

//-----------------------------------------------------------------------------
bool FEMortarSlidingSurface::Init()
{
	// always intialize base class first!
	if (FEMortarContactSurface::Init() == false) return false;

	// get the number of nodes
	int NN = Nodes();

	// allocate data structures
	m_p.resize(NN, 0.0);
	m_L.resize(NN, 0.0);
	m_nu.resize(NN, vec3d(0,0,0));
	m_norm0.resize(NN, 0.0);

	return true;
}

//-----------------------------------------------------------------------------
void FEMortarSlidingSurface::UpdateNormals(bool binit)
{
	const int NN = Nodes();

	// reset normals
	m_nu.assign(NN, vec3d(0,0,0));

	// recalculate normals
	const int NF = Elements();
	for (int i=0; i<NF; ++i)
	{
		FESurfaceElement& el = Element(i);
		int neln = el.Nodes();
		for (int j=0; j<neln; ++j)
		{
			vec3d r0 = Node(el.m_lnode[j]).m_rt;
			vec3d r1 = Node(el.m_lnode[(j+1)%neln]).m_rt;
			vec3d r2 = Node(el.m_lnode[(j+neln-1)%neln]).m_rt;

			vec3d n = (r1 - r0)^(r2 - r0);
			m_nu[el.m_lnode[j]] += n;
		}
	}

	// normalize
	if (binit)
	{
		for (int i=0; i<NN; ++i) 
		{
			double d = m_nu[i].unit();
			m_norm0[i] = 1.0/d;
		}
	}
	else
	{
		for (int i=0; i<NN; ++i) m_nu[i] *= m_norm0[i];
	}
}

//=============================================================================
// FEMortarSlidingContact
//=============================================================================

//-----------------------------------------------------------------------------
// Define sliding interface parameters
BEGIN_FECORE_CLASS(FEMortarSlidingContact, FEMortarInterface)
	ADD_PARAMETER(m_atol   , "tolerance"    );
	ADD_PARAMETER(m_eps    , "penalty"      );
	ADD_PARAMETER(m_naugmin, "minaug"       );
	ADD_PARAMETER(m_naugmax, "maxaug"       );
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEMortarSlidingContact::FEMortarSlidingContact(FEModel* pfem) : FEMortarInterface(pfem), m_ss(pfem), m_ms(pfem)
{
	m_dofX = GetDOFIndex("x");
	m_dofY = GetDOFIndex("y");
	m_dofZ = GetDOFIndex("z");
}

//-----------------------------------------------------------------------------
FEMortarSlidingContact::~FEMortarSlidingContact()
{
}

//-----------------------------------------------------------------------------
bool FEMortarSlidingContact::Init()
{
	// initialize surfaces
	if (m_ms.Init() == false) return false;
	if (m_ss.Init() == false) return false;

	return true;
}

//-----------------------------------------------------------------------------
void FEMortarSlidingContact::Activate()
{
	//! don't forget the base class
	FEContactInterface::Activate();

	// update the normals on the primary surface
	m_ss.UpdateNormals(true);

	// update nodal areas
	m_ss.UpdateNodalAreas();

	// update the mortar weights
	UpdateMortarWeights(m_ss, m_ms);

	// update the nodal gaps
	// (must be done after mortar eights are updated)
	UpdateNodalGaps(m_ss, m_ms);
}

//-----------------------------------------------------------------------------
//! build the matrix profile for use in the stiffness matrix
void FEMortarSlidingContact::BuildMatrixProfile(FEGlobalMatrix& K)
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
void FEMortarSlidingContact::LoadVector(FEGlobalVector& R, const FETimeInfo& tp)
{
	int NS = m_ss.Nodes();
	int NM = m_ms.Nodes();

	// loop over all primary nodes
	for (int A=0; A<NS; ++A)
	{
		vec3d nuA = m_ss.m_nu[A];
		vec3d gA = m_ss.m_gap[A];
		double eps = m_eps*m_ss.m_A[A];
		double gap = gA*nuA;
		double pA = m_ss.m_L[A] + eps*gap;
//		if (gap < 0.0) pA = 0.0;
		
		vec3d tA = nuA*(pA);

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
void FEMortarSlidingContact::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp)
{
	ContactGapStiffness(LS);
	ContactNormalStiffness(LS);
}

//-----------------------------------------------------------------------------
//! calculate contact stiffness
void FEMortarSlidingContact::ContactGapStiffness(FELinearSystem& LS)
{
	int NS = m_ss.Nodes();
	int NM = m_ms.Nodes();

	// A. Linearization of the gap function
	vector<int> lmi(3), lmj(3);
	matrix kA(3, 3), kG(3, 3);
	FEElementMatrix ke;
	ke.resize(3, 3);
	for (int A=0; A<NS; ++A)
	{
		vec3d nuA = m_ss.m_nu[A];
		double eps = m_eps*m_ss.m_A[A];

		// loop over all primary nodes
		for (int B=0; B<NS; ++B)
		{
			FENode& nodeB = m_ss.Node(B);
			lmi[0] = nodeB.m_ID[0];
			lmi[1] = nodeB.m_ID[1];
			lmi[2] = nodeB.m_ID[2];

			double nAB = m_n1[A][B];
			if (nAB != 0.0)
			{
				kA[0][0] = eps*nAB*(nuA.x*nuA.x); kA[0][1] = eps*nAB*(nuA.x*nuA.y); kA[0][2] = eps*nAB*(nuA.x*nuA.z);
				kA[1][0] = eps*nAB*(nuA.y*nuA.x); kA[1][1] = eps*nAB*(nuA.y*nuA.y); kA[1][2] = eps*nAB*(nuA.y*nuA.z);
				kA[2][0] = eps*nAB*(nuA.z*nuA.x); kA[2][1] = eps*nAB*(nuA.z*nuA.y); kA[2][2] = eps*nAB*(nuA.z*nuA.z);

				// loop over primary nodes
				for (int C=0; C<NS; ++C)
				{
					FENode& nodeC = m_ss.Node(C);
					lmj[0] = nodeC.m_ID[0];
					lmj[1] = nodeC.m_ID[1];
					lmj[2] = nodeC.m_ID[2];

					double nAC = m_n1[A][C];
					if (nAC != 0.0)
					{
						kG[0][0] = nAC; kG[0][1] = 0.0; kG[0][2] = 0.0;
						kG[1][0] = 0.0; kG[1][1] = nAC; kG[1][2] = 0.0;
						kG[2][0] = 0.0; kG[2][1] = 0.0; kG[2][2] = nAC;

						ke = kA*kG;

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

					double nAC = -m_n2[A][C];
					if (nAC != 0.0)
					{
						kG[0][0] = nAC; kG[0][1] = 0.0; kG[0][2] = 0.0;
						kG[1][0] = 0.0; kG[1][1] = nAC; kG[1][2] = 0.0;
						kG[2][0] = 0.0; kG[2][1] = 0.0; kG[2][2] = nAC;

						ke = kA*kG;

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

			double nAB = -m_n2[A][B];
			if (nAB != 0.0)
			{
				kA[0][0] = eps*nAB*(nuA.x*nuA.x); kA[0][1] = eps*nAB*(nuA.x*nuA.y); kA[0][2] = eps*nAB*(nuA.x*nuA.z);
				kA[1][0] = eps*nAB*(nuA.y*nuA.x); kA[1][1] = eps*nAB*(nuA.y*nuA.y); kA[1][2] = eps*nAB*(nuA.y*nuA.z);
				kA[2][0] = eps*nAB*(nuA.z*nuA.x); kA[2][1] = eps*nAB*(nuA.z*nuA.y); kA[2][2] = eps*nAB*(nuA.z*nuA.z);

				// loop over primary nodes
				for (int C=0; C<NS; ++C)
				{
					FENode& nodeC = m_ss.Node(C);
					lmj[0] = nodeC.m_ID[0];
					lmj[1] = nodeC.m_ID[1];
					lmj[2] = nodeC.m_ID[2];

					double nAC = m_n1[A][C];
					if (nAC != 0.0)
					{
						kG[0][0] = nAC; kG[0][1] = 0.0; kG[0][2] = 0.0;
						kG[1][0] = 0.0; kG[1][1] = nAC; kG[1][2] = 0.0;
						kG[2][0] = 0.0; kG[2][1] = 0.0; kG[2][2] = nAC;

						ke = kA*kG;

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

					double nAC = -m_n2[A][C];
					if (nAC != 0.0)
					{
						kG[0][0] = nAC; kG[0][1] = 0.0; kG[0][2] = 0.0;
						kG[1][0] = 0.0; kG[1][1] = nAC; kG[1][2] = 0.0;
						kG[2][0] = 0.0; kG[2][1] = 0.0; kG[2][2] = nAC;

						ke = kA*kG;

						ke.SetIndices(lmi, lmj);
						LS.Assemble(ke);
					}
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! calculate contact stiffness
void FEMortarSlidingContact::ContactNormalStiffness(FELinearSystem& LS)
{
	int NS = m_ss.Nodes();
	int NM = m_ms.Nodes();

	vector<int> lm1(3);
	vector<int> lm2(3);
	FEElementMatrix ke;
	ke.resize(3, 3);
	int NF = m_ss.Elements();
	for (int i=0; i<NF; ++i)
	{
		FESurfaceElement& f = m_ss.Element(i);
		int nn = f.Nodes();
		for (int j=0; j<nn; ++j)
		{
			int jp1 = (j+1)%nn;
			int jm1 = (j+nn-1)%nn;
			int A = f.m_lnode[j];

			vec3d vA = m_ss.m_nu[A];
			vec3d gA = m_ss.m_gap[A];
			double eps = m_eps*m_ss.m_A[A];
			double pA = m_ss.m_L[A] + eps*(gA*vA);
			double normA = m_ss.m_norm0[A];

			mat3d kA = ((vA&gA)*eps + mat3dd(pA));

			FENode& nodej1 = m_ss.Node(f.m_lnode[jp1]);
			vec3d& x1 = nodej1.m_rt;
			mat3da k1(x1);
			lm1[0] = nodej1.m_ID[0];
			lm1[1] = nodej1.m_ID[1];
			lm1[2] = nodej1.m_ID[2];

			FENode& nodej2 = m_ss.Node(f.m_lnode[jm1]);
			vec3d& x2 = nodej2.m_rt;
			mat3da k2(x2);
			lm2[0] = nodej2.m_ID[0];
			lm2[1] = nodej2.m_ID[1];
			lm2[2] = nodej2.m_ID[2];

			// loop over primary nodes
			for (int B=0; B<NS; ++B)
			{
				FENode& nodeB = m_ss.Node(B);
				
				double nAB = m_n1[A][B];
				if (nAB != 0.0)
				{
					vector<int> lmi(3);
					lmi[0] = nodeB.m_ID[0];
					lmi[1] = nodeB.m_ID[1];
					lmi[2] = nodeB.m_ID[2];

					mat3d kab = (kA*k1)*(nAB*normA);
					ke[0][0] = kab(0,0); ke[0][1] = kab(0,1); ke[0][2] = kab(0,2);
					ke[1][0] = kab(1,0); ke[1][1] = kab(1,1); ke[1][2] = kab(1,2);
					ke[2][0] = kab(2,0); ke[2][1] = kab(2,1); ke[2][2] = kab(2,2);

					ke.SetIndices(lmi, lm2);
					LS.Assemble(ke);

					kab = (kA*k2)*(-nAB*normA);
					ke[0][0] = kab(0,0); ke[0][1] = kab(0,1); ke[0][2] = kab(0,2);
					ke[1][0] = kab(1,0); ke[1][1] = kab(1,1); ke[1][2] = kab(1,2);
					ke[2][0] = kab(2,0); ke[2][1] = kab(2,1); ke[2][2] = kab(2,2);

					ke.SetIndices(lmi, lm1);
					LS.Assemble(ke);
				}
			}

			// loop over secondary nodes
			for (int B=0; B<NM; ++B)
			{
				FENode& nodeB = m_ms.Node(B);
				
				double nAB = m_n2[A][B];
				if (nAB != 0.0)
				{
					vector<int> lmi(3);
					lmi[0] = nodeB.m_ID[0];
					lmi[1] = nodeB.m_ID[1];
					lmi[2] = nodeB.m_ID[2];

					mat3d kab = (kA*k1)*(nAB*normA);
					ke[0][0] = kab(0,0); ke[0][1] = kab(0,1); ke[0][2] = kab(0,2);
					ke[1][0] = kab(1,0); ke[1][1] = kab(1,1); ke[1][2] = kab(1,2);
					ke[2][0] = kab(2,0); ke[2][1] = kab(2,1); ke[2][2] = kab(2,2);

					ke.SetIndices(lmi, lm2);
					LS.Assemble(ke);

					kab = (kA*k2)*(-nAB*normA);
					ke[0][0] = kab(0,0); ke[0][1] = kab(0,1); ke[0][2] = kab(0,2);
					ke[1][0] = kab(1,0); ke[1][1] = kab(1,1); ke[1][2] = kab(1,2);
					ke[2][0] = kab(2,0); ke[2][1] = kab(2,1); ke[2][2] = kab(2,2);

					ke.SetIndices(lmi, lm1);
					LS.Assemble(ke);
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! calculate Lagrangian augmentations
bool FEMortarSlidingContact::Augment(int naug, const FETimeInfo& tp)
{
	if (m_laugon != 1) return true;

	double max_err = 0.0;
	int NS = m_ss.Nodes();
	// loop over all primary nodes
	for (int A=0; A<NS; ++A)
	{
		vec3d vA = m_ss.m_nu[A];
		vec3d gA = m_ss.m_gap[A];
		double gap = gA*vA;
		double eps = m_eps*m_ss.m_A[A];
		
		double Lold = m_ss.m_L[A];
		double Lnew = Lold + eps*gap;

		double err = fabs((Lold - Lnew)/(Lold + Lnew));
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
			vec3d vA = m_ss.m_nu[A];
			vec3d gA = m_ss.m_gap[A];
			double gap = gA*vA;
			double eps = m_eps*m_ss.m_A[A];
		
			double Lold = m_ss.m_L[A];
			double Lnew = Lold + eps*gap;
			m_ss.m_L[A] = Lnew;
		}
	}

	return bconv;
}

//-----------------------------------------------------------------------------
//! update interface data
void FEMortarSlidingContact::Update()
{
	m_ss.UpdateNormals(false);
	UpdateMortarWeights(m_ss, m_ms);
	UpdateNodalGaps(m_ss, m_ms);
}

//-----------------------------------------------------------------------------
//! serialize data to archive
void FEMortarSlidingContact::Serialize(DumpStream& ar)
{
}
