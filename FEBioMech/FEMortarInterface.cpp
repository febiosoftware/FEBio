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
#include "FEMortarInterface.h"
#include "FECore/mortar.h"
#include "FECore/log.h"
#include <FECore/FEMesh.h>

//-----------------------------------------------------------------------------
FEMortarInterface::FEMortarInterface(FEModel* pfem) : FEContactInterface(pfem)
{
	// set the integration rule
	m_pT = dynamic_cast<FESurfaceElementTraits*>(FEElementLibrary::GetElementTraits(FE_TRI3G7));
}

//-----------------------------------------------------------------------------
void FEMortarInterface::UpdateMortarWeights(FESurface& ss, FESurface& ms)
{
	// allocate sturcture for the integration weights
	int NS = ss.Nodes();
	int NM = ms.Nodes();
	m_n1.resize(NS,NS);
	m_n2.resize(NS,NM);

	// clear weights
	m_n1.zero();
	m_n2.zero();

	// number of integration points
	const int MAX_INT = 11;
	const int nint = m_pT->m_nint;
	vector<double>& gw = m_pT->gw;
	vector<double>& gr = m_pT->gr;
	vector<double>& gs = m_pT->gs;

	// calculate the mortar surface
	MortarSurface mortar;
	CalculateMortarSurface(ss, ms, mortar);

	// These arrays will store the shape function values of the projection points 
	// on the primary and secondary side when evaluating the integral over a pallet
	double Ns[MAX_INT][4], Nm[MAX_INT][4];

	// loop over the mortar patches
	int NP = mortar.Patches();
	for (int i=0; i<NP; ++i)
	{
		// get the next patch
		Patch& pi = mortar.GetPatch(i);

		// get the facet ID's that generated this patch
		int k = pi.GetPrimaryFacetID();
		int l = pi.GetSecondaryFacetID();

		// get the non-mortar surface element
		FESurfaceElement& se = ss.Element(k);
		// get the mortar surface element
		FESurfaceElement& me = ms.Element(l);

		// loop over all patch triangles
		int np = pi.Size();
		for (int j=0; j<np; ++j)
		{
			// get the next facet
			Patch::FACET& fj = pi.Facet(j);

			// calculate the patch area
			// (We multiply by two because the sum of the integration weights in FEBio sum up to the area
			// of the triangle in natural coordinates (=0.5)).
			double Area = fj.Area()*2.0;
			if (Area > 1e-15)
			{
				// loop over integration points
				for (int n=0; n<nint; ++n)
				{
					// evaluate the spatial position of the integration point on the patch
					vec3d xp = fj.Position(gr[n], gs[n]);

					// evaluate the integration points on the primary and secondary surfaces
					// i.e. determine rs, rm
					double r1 = 0, s1 = 0, r2 = 0, s2 = 0;
					vec3d xs = ss.ProjectToSurface(se, xp, r1, s1);
					vec3d xm = ms.ProjectToSurface(me, xp, r2, s2);

//					assert((r1>=0.0)&&(s1>=0)&&(r1+s1<1.0));
//					assert((r2>=0.0)&&(s2>=0)&&(r2+s2<1.0));

					// evaluate shape functions
					se.shape_fnc(Ns[n], r1, s1);
					me.shape_fnc(Nm[n], r2, s2);
				}

				// Evaluate the contributions to the integrals
				int ns = se.Nodes();
				int nm = me.Nodes();
				for (int A=0; A<ns; ++A)
				{
					int a = se.m_lnode[A];

					// loop over all the nodes on the primary facet
					for (int B=0; B<ns; ++B)
					{
						double n1 = 0;
						for (int n=0; n<nint; ++n)
						{
							n1 += gw[n]*Ns[n][A]*Ns[n][B];
						}
						n1 *= Area;

						int b = se.m_lnode[B];
						m_n1[a][b] += n1;
					}

					// loop over all the nodes on the secondary facet
					for (int C = 0; C<nm; ++C)
					{
						double n2 = 0;
						for (int n=0; n<nint; ++n)
						{
							n2 += gw[n]*Ns[n][A]*Nm[n][C];
						}
						n2 *= Area;

						int c = me.m_lnode[C];
						m_n2[a][c] += n2;
					}
				}
			}		
		}
	}

#ifdef _DEBUG
	// Sanity check: sum should add up to contact area
	// This is for a hardcoded problem. Remove or generalize this!
	double sum1 = 0.0;
	for (int A=0; A<NS; ++A)
		for (int B=0; B<NS; ++B) sum1 += m_n1[A][B];

	double sum2 = 0.0;
	for (int A=0; A<NS; ++A)
		for (int C=0; C<NM; ++C) sum2 += m_n2[A][C];

	if (fabs(sum1 - 1.0) > 1e-5) feLog("WARNING: Mortar weights are not correct (%lg).\n", sum1);
	if (fabs(sum2 - 1.0) > 1e-5) feLog("WARNING: Mortar weights are not correct (%lg).\n", sum2);
#endif
}

//-----------------------------------------------------------------------------
//! Update the nodal gaps
void FEMortarInterface::UpdateNodalGaps(FEMortarContactSurface& ss, FEMortarContactSurface& ms)
{
	// reset nodal gaps
	vector<vec3d>& gap = ss.m_gap;
	zero(ss.m_gap);

	int NS = ss.Nodes();
	int NM = ms.Nodes();

	// loop over all primary nodes
	for (int A=0; A<NS; ++A)
	{
		// loop over all primary nodes
		for (int B=0; B<NS; ++B)
		{
			FENode& nodeB = ss.Node(B);
			vec3d& xB = nodeB.m_rt;
			double nAB = m_n1[A][B];
			gap[A] += xB*nAB;
		}

		// loop over secondary side
		for (int C=0; C<NM; ++C)
		{
			FENode& nodeC = ms.Node(C);
			vec3d& xC = nodeC.m_rt;
			double nAC = m_n2[A][C];
			gap[A] -= xC*nAC;
		}
	}
}
