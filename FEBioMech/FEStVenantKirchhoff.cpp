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
#include "FEStVenantKirchhoff.h"

// define the material parameters
BEGIN_FECORE_CLASS(FEStVenantKirchhoff, FEElasticMaterial)
	ADD_PARAMETER(m_E, FE_RANGE_GREATER(0.0), "E");
	ADD_PARAMETER(m_v, FE_RANGE_RIGHT_OPEN(-1.0, 0.5), "v");
END_FECORE_CLASS();

//////////////////////////////////////////////////////////////////////
// FEStVenantKirchhoff
//////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
mat3ds FEStVenantKirchhoff::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// lame parameters
	double lam = m_v*m_E/((1+m_v)*(1-2*m_v));
	double mu  = 0.5*m_E/(1+m_v);

	// calculate left Cauchy-Green tensor (ie. b-matrix)
	mat3ds b = pt.LeftCauchyGreen();
	mat3ds b2 = b.sqr();

	// calculate trace of Green-Lagrance strain tensor
	double trE = 0.5*(b.tr()-3);

	// inverse jacobian
	double Ji = 1.0 / pt.m_J;

	// calculate stress
	// s = (lam*trE*b - mu*(b2 - b))/J;
	return b*(lam*trE*Ji) + (b2 - b)*(mu*Ji);
}

//-----------------------------------------------------------------------------
tens4ds FEStVenantKirchhoff::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// jacobian
	double J = pt.m_J;

	// lame parameters
	double lam = m_v*m_E/((1+m_v)*(1-2*m_v));
	double mu  = 0.5*m_E/(1+m_v);

	double lam1 = lam / J;
	double mu1  = mu / J;

	// left cauchy-green matrix (i.e. the 'b' matrix)
	mat3ds b = pt.LeftCauchyGreen();

	tens4ds c = dyad1s(b)*lam1 + dyad4s(b)*(2.0*mu1);

	return c;
}

//-----------------------------------------------------------------------------
double FEStVenantKirchhoff::StrainEnergyDensity(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
	// lame parameters
	double lam = m_v*m_E/((1+m_v)*(1-2*m_v));
	double mu  = 0.5*m_E/(1+m_v);
    
	// calculate right Cauchy-Green tensor
	mat3ds C = pt.RightCauchyGreen();
	mat3ds C2 = C.sqr();
    
	double trE = 0.5*(C.tr()-3);
    double trE2 = 0.25*(C2.tr() - 2*C.tr() + 3);
    
	// calculate strain energy density
	return lam*trE*trE/2.0 + mu*trE2;
}
