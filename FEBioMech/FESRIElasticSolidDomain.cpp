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
#include "FESRIElasticSolidDomain.h"
#include "FEElasticMaterial.h"
#include "FECore/FEElementLibrary.h"

//-----------------------------------------------------------------------------
FESRIElasticSolidDomain::FESRIElasticSolidDomain(FEModel* pfem) : FEElasticSolidDomain(pfem)
{
}

//-----------------------------------------------------------------------------
void FESRIElasticSolidDomain::ElementInternalForce(FESolidElement& el, vector<double>& fe)
{
	// jacobian matrix, inverse jacobian matrix and determinants
	double Ji[3][3];

	// first, evaluate the deviatoric component using the full-integration rule
	int nint = el.GaussPoints();
	int neln = el.Nodes();

	double*	gw = el.GaussWeights();

	// repeat for all integration points
	for (int n=0; n<nint; ++n)
	{
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());

		// calculate the jacobian
		double detJt = invjact(el, Ji, n);

		detJt *= gw[n];

		// get the stress vector for this integration point
		// we only use the deviatoric component
		mat3ds s = pt.m_s.dev();

		double* Gr = el.Gr(n);
		double* Gs = el.Gs(n);
		double* Gt = el.Gt(n);

		for (int i=0; i<neln; ++i)
		{
			// calculate global gradient of shape functions
			// note that we need the transposed of Ji, not Ji itself !
			double Gx = Ji[0][0]*Gr[i]+Ji[1][0]*Gs[i]+Ji[2][0]*Gt[i];
			double Gy = Ji[0][1]*Gr[i]+Ji[1][1]*Gs[i]+Ji[2][1]*Gt[i];
			double Gz = Ji[0][2]*Gr[i]+Ji[1][2]*Gs[i]+Ji[2][2]*Gt[i];

			// calculate internal force
			// the '-' sign is so that the internal forces get subtracted
			// from the global residual vector
			fe[3*i  ] -= ( Gx*s.xx() +
				           Gy*s.xy() +
					       Gz*s.xz() )*detJt;

			fe[3*i+1] -= ( Gy*s.yy() +
				           Gx*s.xy() +
					       Gz*s.yz() )*detJt;

			fe[3*i+2] -= ( Gz*s.zz() +
				           Gy*s.yz() +
					       Gx*s.xz() )*detJt;
		}
	}

	// next, the pressure component is evaluated using the reduced integration rule
	const int NME = FEElement::MAX_NODES;
	double Gr[NME], Gs[NME], Gt[NME];
	FESRISolidElementTraits* prt_ri = dynamic_cast<FESRISolidElementTraits*>(el.GetTraits());

	// if this element does not have a reduced rule, we use the full integration rule
	FESolidElementTraits* prt = (prt_ri ? prt_ri->m_pTRI : dynamic_cast<FESolidElementTraits*>(el.GetTraits()));

	nint = prt->m_nint;
	vector<double>& gr = prt->gr;
	vector<double>& gs = prt->gs;
	vector<double>& gt = prt->gt;
	gw = &(prt->gw[0]);
	for (int n=0; n<nint; ++n)
	{
		double r = gr[n];
		double s = gs[n];
		double t = gt[n];

		// setup the material point
		FEElasticMaterialPoint pt;
		FEMaterialPoint mp(&pt);

		pt.m_J = defgrad(el, pt.m_F, r, s, t);

		// calculate the jacobian
		double detJt = invjact(el, Ji, r, s, t);

		// get the shape function derivatives
		el.shape_deriv(Gr, Gs, Gt, r, s, t);

		// calculate the stress
		mat3ds stress = m_pMat->Stress(mp);

		// evaluate the pressure
		double p = -(stress.tr()/3.0);

		detJt *= gw[n];
		for (int i=0; i<neln; ++i)
		{
			// calculate global gradient of shape functions
			// note that we need the transposed of Ji, not Ji itself !
			double Gx = Ji[0][0]*Gr[i]+Ji[1][0]*Gs[i]+Ji[2][0]*Gt[i];
			double Gy = Ji[0][1]*Gr[i]+Ji[1][1]*Gs[i]+Ji[2][1]*Gt[i];
			double Gz = Ji[0][2]*Gr[i]+Ji[1][2]*Gs[i]+Ji[2][2]*Gt[i];

			// calculate internal force
			// the '-' sign is so that the internal forces get subtracted
			// from the global residual vector
			fe[3*i  ] += ( Gx*p )*detJt;
			fe[3*i+1] += ( Gy*p )*detJt;
			fe[3*i+2] += ( Gz*p )*detJt;
		}
	}
}
