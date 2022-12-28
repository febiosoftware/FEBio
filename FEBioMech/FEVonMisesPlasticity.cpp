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
#include "FEVonMisesPlasticity.h"

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FEVonMisesPlasticity, FESolidMaterial)
	ADD_PARAMETER(m_E, FE_RANGE_GREATER(0.0), "E");
	ADD_PARAMETER(m_v, FE_RANGE_RIGHT_OPEN(-1.0, 0.5), "v");
	ADD_PARAMETER(m_Y, FE_RANGE_GREATER(0.0), "Y");
	ADD_PARAMETER(m_H, FE_RANGE_GREATER_OR_EQUAL(0.0), "H");
END_FECORE_CLASS();


//-----------------------------------------------------------------------------
FEVonMisesPlasticity::FEVonMisesPlasticity(FEModel* pfem) : FESolidMaterial(pfem)
{
	m_E = m_v = m_Y = m_H = 0;
	m_K = m_G = 0;
}

//-----------------------------------------------------------------------------
bool FEVonMisesPlasticity::Init()
{
	if (FESolidMaterial::Init() == false) return false;

	m_K = m_E/(3.0*(1.0 - 2*m_v));
	m_G = m_E/(2.0*(1.0 +   m_v));

	return true;
}

//-----------------------------------------------------------------------------
FEMaterialPointData* FEVonMisesPlasticity::CreateMaterialPointData()
{
	FEJ2PlasticMaterialPoint* pt = new FEJ2PlasticMaterialPoint;
	pt->Y0 = m_Y;
	return pt;
}

//-----------------------------------------------------------------------------
mat3ds FEVonMisesPlasticity::Stress(FEMaterialPoint &mp)
{
	FEJ2PlasticMaterialPoint& pp = *mp.ExtractData<FEJ2PlasticMaterialPoint>();
	mat3d& F = pp.m_F;
	
	// get the current strain
	mat3ds e = F.sym() - mat3dd(1.0);

	// calculate strain increment
	mat3ds de = (e - pp.e0);

	// get the trial stress
	mat3ds strial = pp.sn + (de.dev()*(2.0*m_G) + de.iso()*(3.0*m_K));
	mat3ds dev_strial = strial.dev();
	double devs_norm = dev_strial.norm();

	// get current yield strenght
	double Y = pp.Y0;

	double k = Y / sqrt(3.0);
	double fac = devs_norm / (sqrt(2.0)*k);

	mat3ds s;
	if (fac<=1)
	{
		s = strial;
		pp.b = false;
	}
	else
	{
		// calculate plastic strain rate
		double L = (devs_norm - sqrt(2.0/3.0)*Y)/(2*m_G + m_H);

		// update yield strength
		pp.Y1 = Y + m_H*L/sqrt(2.0/3.0);

		// update stress
		s = strial.iso() + dev_strial*(1.0 - 2.0*m_G*L/devs_norm);
		pp.b = true;
	}

	// store the current strain measure
	pp.e1 = e;

	return s;
}

//-----------------------------------------------------------------------------
tens4ds FEVonMisesPlasticity::Tangent(FEMaterialPoint &mp)
{
	FEJ2PlasticMaterialPoint& pp = *mp.ExtractData<FEJ2PlasticMaterialPoint>();

	// lame parameters
	double lam = m_K - m_G*2.0/3.0;
	double mu  = m_G;

	double D[6][6] = {0};
	D[0][0] = lam+2.*mu; D[0][1] = lam      ; D[0][2] = lam      ;
	D[1][0] = lam      ; D[1][1] = lam+2.*mu; D[1][2] = lam      ;
	D[2][0] = lam      ; D[2][1] = lam      ; D[2][2] = lam+2.*mu;
	D[3][3] = mu;
	D[4][4] = mu;
	D[5][5] = mu;
	tens4ds C(D);

	// see if we are in plastic flow mode
	if (pp.b)
	{
		// get the stress
		mat3ds s = pp.m_s;
		mat3ds n = s.dev()*2.0;

		mat3ds A = C.dot(n);
		double G = n.dotdot(A) + m_H;

		C -= dyad4s(A)/G;
	}

	return C;
}
