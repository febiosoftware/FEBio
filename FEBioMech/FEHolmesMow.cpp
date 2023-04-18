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
#include "FEHolmesMow.h"

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FEHolmesMow, FEElasticMaterial)
	ADD_PARAMETER(m_E, FE_RANGE_GREATER(0.0), "E")->setUnits(UNIT_PRESSURE)->setLongName("Young's modulus");
	ADD_PARAMETER(m_v, FE_RANGE_RIGHT_OPEN(-1.0, 0.5), "v")->setLongName("Poisson's ratio");
	ADD_PARAMETER(m_b, FE_RANGE_GREATER_OR_EQUAL(0.0), "beta")->setLongName("power exponent");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
bool FEHolmesMow::Validate()
{
	if (FEElasticMaterial::Validate() == false) return false;
	
	// Lame coefficients
	lam = m_v*m_E/((1+m_v)*(1-2*m_v));
	mu  = 0.5*m_E/(1+m_v);
	Ha = lam + 2*mu;	

	return true;
}

//-----------------------------------------------------------------------------
mat3ds FEHolmesMow::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	double detF = pt.m_J;
	double detFi = 1.0/detF;
	
	// calculate left Cauchy-Green tensor
	mat3ds b = pt.LeftCauchyGreen(); //(F*F.transpose()).sym();
	mat3ds b2 = b.sqr();
	mat3ds identity(1.,1.,1.,0.,0.,0.);

	// calculate invariants of B
	double I1 = b.tr();
	double I2 = (I1*I1 - b2.tr())/2.;
	double I3 = b.det();

	// Exponential term
	double eQ = exp(m_b*((2*mu-lam)*(I1-3) + lam*(I2-3))/Ha)/pow(I3,m_b);
	
	// calculate stress
	mat3ds s = 0.5*detFi*eQ*((2*mu+lam*(I1-1))*b - lam*b2 - Ha*identity);

	return s;
}

//-----------------------------------------------------------------------------
tens4ds FEHolmesMow::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	double detF = pt.m_J;
	double detFi = 1.0/detF;
	
	// calculate left Cauchy-Green tensor
	mat3ds b = pt.LeftCauchyGreen(); //(F*F.transpose()).sym();
	mat3ds b2 = b.sqr();
	mat3ds identity(1.,1.,1.,0.,0.,0.);
	
	// calculate invariants of B
	double I1 = b.tr();
	double I2 = (I1*I1 - b2.tr())*0.5;
	double I3 = b.det();
	
	// Exponential term
	double eQ = exp(m_b*((2*mu-lam)*(I1-3) + lam*(I2-3))/Ha)/pow(I3,m_b);
	
	// calculate stress
	mat3ds s = 0.5*detFi*eQ*((2*mu+lam*(I1-1))*b - lam*b2 - Ha*identity);
	
	// calculate elasticity tensor
	tens4ds c = 4.*m_b/Ha*detF/eQ*dyad1s(s) 
	+ detFi*eQ*(lam*(dyad1s(b) - dyad4s(b)) + Ha*dyad4s(identity));
	
	return c;
}

//-----------------------------------------------------------------------------
double FEHolmesMow::StrainEnergyDensity(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
	// calculate left Cauchy-Green tensor
	mat3ds b = pt.LeftCauchyGreen(); //(F*F.transpose()).sym();
	mat3ds b2 = b.sqr();
    
	// calculate invariants of B
	double I1 = b.tr();
	double I2 = (I1*I1 - b2.tr())/2.;
	double I3 = b.det();
    
	// Exponential term
	double eQ = exp(m_b*((2*mu-lam)*(I1-3) + lam*(I2-3))/Ha)/pow(I3,m_b);
	
	// calculate strain energy density
	double sed = Ha/(4*m_b)*(eQ - 1);
	
	return sed;
}
