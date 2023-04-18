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
#include "FECoupledMooneyRivlin.h"

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FECoupledMooneyRivlin, FEElasticMaterial)
	ADD_PARAMETER(m_c1, FE_RANGE_GREATER(0.0), "c1")->setUnits(UNIT_PRESSURE);
	ADD_PARAMETER(m_c2, "c2")->setUnits(UNIT_PRESSURE);
	ADD_PARAMETER(m_K , FE_RANGE_GREATER(0.0), "k" )->setUnits(UNIT_PRESSURE);
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! calculate stress at material point
mat3ds FECoupledMooneyRivlin::Stress(FEMaterialPoint& mp)
{
	// get the elastic material point data
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// determinant of deformation gradient
	double J = pt.m_J;

	// calculate left Cauchy-Green tensor
	mat3ds B = pt.LeftCauchyGreen();

	// calculate square of B
	mat3ds B2 = B.sqr();

	// Invariants of B (= invariants of C)
	double I1 = B.tr();

	// identity tensor
	mat3dd I(1.0);

	// calculate stress
	return (B*(m_c1+I1*m_c2) - B2*m_c2 - I*(m_c1+2.0*m_c2))*(2.0/J) + I*(m_K*log(J)/J);
}

//-----------------------------------------------------------------------------
//! calculate tangent at material point
tens4ds FECoupledMooneyRivlin::Tangent(FEMaterialPoint& mp)
{
	// get material point data
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// calculate left Cauchy-Green tensor: B = F*Ft
	mat3ds B = pt.LeftCauchyGreen();

	// Invariants of B (= invariants of C)
	double J = pt.m_J;

	// some useful tensors
	mat3dd I(1.0);
	tens4ds IxI = dyad1s(I);
	tens4ds IoI = dyad4s(I);
	tens4ds BxB = dyad1s(B);
	tens4ds BoB = dyad4s(B);

	// strain energy derivates
	double W2 = m_c2;

	// spatial tangent
	tens4ds c = BxB*(4.0*W2/J) - BoB*(4.0*W2/J) + IoI*(4.0*(m_c1+2.0*m_c2)/J) + IxI*(m_K/J) - IoI*(2.0*m_K*log(J)/J);

	return c;
}

//-----------------------------------------------------------------------------
//! calculate strain energy density at material point
double FECoupledMooneyRivlin::StrainEnergyDensity(FEMaterialPoint& mp)
{
	// get the elastic material point data
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
	// determinant of deformation gradient
	double J = pt.m_J;
    double lnJ = log(J);
    
	// calculate left Cauchy-Green tensor
	mat3ds B = pt.LeftCauchyGreen();
    
	// calculate square of B
	mat3ds B2 = B.sqr();
    
	// Invariants of B (= invariants of C)
	double I1 = B.tr();
    double I2 = (I1*I1 - B2.tr())/2.;
    
    double sed = m_c1*(I1-3) + m_c2*(I2-3)
    - 2*(m_c1+2*m_c2)*lnJ + m_K*lnJ*lnJ/2.;
    
    return sed;
}
