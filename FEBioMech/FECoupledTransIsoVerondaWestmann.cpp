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
#include "FECoupledTransIsoVerondaWestmann.h"
#include <FECore/log.h>
#include <FECore/expint_Ei.h>
#include <FECore/FEConstValueVec3.h>

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FECoupledTransIsoVerondaWestmann, FEElasticMaterial)
	ADD_PARAMETER(m_c1  , FE_RANGE_GREATER(0.0), "c1")->setUnits(UNIT_PRESSURE);
	ADD_PARAMETER(m_c2  , "c2");
	ADD_PARAMETER(m_c3  , "c3")->setUnits(UNIT_PRESSURE);
	ADD_PARAMETER(m_c4  , "c4");
	ADD_PARAMETER(m_c5  , "c5");
	ADD_PARAMETER(m_flam, FE_RANGE_GREATER_OR_EQUAL(1.0), "lambda");
	ADD_PARAMETER(m_K   , FE_RANGE_GREATER(0.0), "k");
	
	ADD_PROPERTY(m_fiber, "fiber");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FECoupledTransIsoVerondaWestmann::FECoupledTransIsoVerondaWestmann(FEModel* pfem) : FEElasticMaterial(pfem)
{
	m_c1 = 0.0;
	m_c2 = 0.0;
	m_c3 = 0.0;
	m_c4 = 0.0;
	m_c5 = 0.0;
	m_flam = 1.0;
	m_K = 0.0;

	m_fiber = nullptr;
}

//-----------------------------------------------------------------------------
//! Calculate the Cauchy stress
mat3ds FECoupledTransIsoVerondaWestmann::Stress(FEMaterialPoint& mp)
{
	// get the material point data
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// calculate left Cauchy-Green tensor
	mat3ds B = pt.LeftCauchyGreen();

	// calculate square of B
	mat3ds B2 = B.sqr();

	// get the material fiber axis
	vec3d a0 = m_fiber->unitVector(mp);

	// get the spatial fiber axis
	vec3d a = pt.m_F*a0;
	double l = a.unit();

	// Invariants of B (= invariants of C)
	double I1 = B.tr();
	double J = pt.m_J;

	// some useful tensors
	mat3dd I(1.0);
	mat3ds A = dyad(a);

	// a. define the matrix stress
	//-----------------------------
	// W = C1*(exp(C2*(I1-3)-1)-0.5*C1*C2*(I2 - 3)
	// Wi = dW/dIi
	double W1 = m_c1*m_c2*exp(m_c2*(I1-3));
	double W2 = -0.5*m_c1*m_c2;

	mat3ds s = (B*(W1 + I1*W2) - B2*W2)*(2.0/J);

	// b. define fiber stress
	//-------------------------
	if (l > 1.0)
	{
		double Wl = 0.0;
		if (l < m_flam)
		{
			Wl = m_c3*(exp(m_c4*(l - 1.0)) - 1.0);
		}
		else
		{
			double c6 = m_c3*(exp(m_c4*(m_flam - 1.0)) - 1.0) - m_c5*m_flam;
			Wl = m_c5*l + c6;
		}
		s += A*(Wl/J);
	}

	// c. define dilational stress
	//------------------------------
	// U(J) = 1/2*k*(lnJ)^2
	double UJ = m_K*log(J)/J;
	s += I*UJ;

	return s;
}

//-----------------------------------------------------------------------------
//! Calculate the spatial elasticity tangent
tens4ds FECoupledTransIsoVerondaWestmann::Tangent(FEMaterialPoint& mp)
{
	// get material point data
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// calculate left Cauchy-Green tensor: B = F*Ft
	mat3ds B = pt.LeftCauchyGreen();

	// get the material fiber axis
	vec3d a0 = m_fiber->unitVector(mp);

	// get the spatial fiber axis
	vec3d a = pt.m_F*a0;
	double l = a.unit();

	// Invariants of B (= invariants of C)
	double I1 = B.tr();
	double J = pt.m_J;
	double I4 = l*l;

	// some useful tensors
	mat3dd I(1.0);
	mat3ds A = dyad(a);
	tens4ds IxI = dyad1s(I);
	tens4ds IoI = dyad4s(I);
	tens4ds BxB = dyad1s(B);
	tens4ds BoB = dyad4s(B);
	tens4ds AxA = dyad1s(A);

	// a. matrix tangent
	//----------------------------------
	double W1 = m_c1*m_c2*exp(m_c2*(I1-3));
	double W11 = m_c2*W1;
	double W2 = -0.5*m_c1*m_c2;
	tens4ds c = BxB*(4.0*(W11+W2)/J) - BoB*(4.0*W2/J);

	// b. fiber tangent
	// ---------------------------------
	if (l > 1.0)
	{
		double Fl = 0.0, Fll = 0.0;
		if (l < m_flam)
		{
			Fl = m_c3*(exp(m_c4*(l - 1.0)) - 1.0)/l;
			Fll = -m_c3*(exp(m_c4*(l-1.0)) - 1.0)/(l*l) + m_c3*m_c4*exp(m_c4*(l-1.0))/l;
		}
		else
		{
			double c6 = m_c3*(exp(m_c4*(m_flam - 1.0)) - 1.0) - m_c5*m_flam;
			Fl = m_c5 + c6 / l;
			Fll = -c6/(l*l);
		}

		double W44 = (Fll - Fl/l)/(4*l*l);

		c += AxA*(4.0*W44*I4*I4/J);
	}

	// c. dilational tangent
	// ---------------------------------
	double UJ = m_K*log(J)/J;
	double UJJ = m_K*(1 - log(J))/(J*J);
	c += IxI*(UJ + J*UJJ) - IoI*(2.0*UJ);

	return c;
}

//-----------------------------------------------------------------------------
//! calculate strain energy density at material point
double FECoupledTransIsoVerondaWestmann::StrainEnergyDensity(FEMaterialPoint& mp)
{
	// get the material point data
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
	// calculate left Cauchy-Green tensor
	mat3ds B = pt.LeftCauchyGreen();
    
	// calculate square of B
	mat3ds B2 = B.sqr();
    
	// get the material fiber axis
	vec3d a0 = m_fiber->unitVector(mp);

	// get the spatial fiber axis
	vec3d a = pt.m_F*a0;
	double l = a.unit();
    
	// Invariants of B (= invariants of C)
	double I1 = B.tr();
    double I2 = 0.5*(I1*I1 - B2.tr());
	double J = pt.m_J;
    double lnJ = log(J);
    
	// a. define the matrix strain energy density
	//-----------------------------
	// W = C1*(exp(C2*(I1-3)-1)-0.5*C1*C2*(I2 - 3)
    double sed = m_c1*(exp(m_c2*(I1-3))-1) - m_c1*m_c2*(I2-3)/2;
    
	// b. define fiber strain energy density
	//-------------------------
	if (l > 1.0)
	{
		if (l < m_flam)
		{
            sed += m_c3*(exp(-m_c4)*(expint_Ei(m_c4*l) - expint_Ei(m_c4))-log(l));
		}
		else
		{
			double c6 = m_c3*(exp(m_c4*(m_flam - 1.0)) - 1.0) - m_c5*m_flam;
			sed += m_c5*(l-1) +c6*log(l);
		}
	}
    
	// c. add dilational strain energy density
	//------------------------------
	// U(J) = 1/2*k*(lnJ)^2
    sed += m_K*lnJ*lnJ/2;
    
    return sed;
}
