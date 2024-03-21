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
#include "FEYeoh.h"

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FEYeoh, FEUncoupledMaterial)
	ADD_PARAMETER(m_c[0], "c1")->setUnits(UNIT_PRESSURE);
	ADD_PARAMETER(m_c[1], "c2")->setUnits(UNIT_PRESSURE);
    ADD_PARAMETER(m_c[2], "c3")->setUnits(UNIT_PRESSURE);
    ADD_PARAMETER(m_c[3], "c4")->setUnits(UNIT_PRESSURE);
    ADD_PARAMETER(m_c[4], "c5")->setUnits(UNIT_PRESSURE);
    ADD_PARAMETER(m_c[5], "c6")->setUnits(UNIT_PRESSURE);
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Calculate the deviatoric stress
mat3ds FEYeoh::DevStress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

    double c[MAX_TERMS];
    
	// get material parameters
    for (int i=0; i<MAX_TERMS; ++i) c[i] = m_c[i](mp);

	// determinant of deformation gradient
	double J = pt.m_J;

	// calculate deviatoric left Cauchy-Green tensor
	mat3ds B = pt.DevLeftCauchyGreen();

	// Invariant of B tilde (= invariants of C tilde)
	double I1 = B.tr();

    double sum = 0;
    for (int i=0; i<MAX_TERMS; ++i)
        sum += c[i]*pow(I1-3,i);
    
	// calculate sigma tilde
	mat3ds sig = B*(sum*2.0/J);

	return sig.dev();
}

//-----------------------------------------------------------------------------
//! Calculate the deviatoric tangent
tens4ds FEYeoh::DevTangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

    double c[MAX_TERMS];
    
    // get material parameters
    for (int i=0; i<MAX_TERMS; ++i) c[i] = m_c[i](mp);
    
    // determinant of deformation gradient
    double J = pt.m_J;
    
    // calculate deviatoric left Cauchy-Green tensor
    mat3ds B = pt.DevLeftCauchyGreen();
    
    // Invariant of B tilde (= invariants of C tilde)
    double I1 = B.tr();
    
    double sum = 0;
    for (int i=0; i<MAX_TERMS; ++i)
        sum += c[i]*pow(I1-3,i);
    double csum = 0;
    for (int i=1; i<MAX_TERMS; ++i)
        csum += c[i]*pow(I1-3,i-1);

    // calculate sigma tilde
    mat3ds sd = (B*(sum*2.0/J));
    
    // identity tensor
    mat3dd I(1);
    tens4ds IxI = dyad1s(I);
    tens4ds I4 = dyad4s(I);
    
	tens4ds BxB = dyad1s(B);

	tens4ds C = BxB*(csum*4.0/J);

    C += - 1./3.*(ddots(C,IxI) - IxI*(C.tr()/3.))
    + 2./3.*((I4-IxI/3.)*sd.tr()-dyad1s(sd.dev(),I));

	return C;
}

//-----------------------------------------------------------------------------
//! calculate deviatoric strain energy density
double FEYeoh::DevStrainEnergyDensity(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    double c[MAX_TERMS];
    
    // get material parameters
    for (int i=0; i<MAX_TERMS; ++i) c[i] = m_c[i](mp);
    
    // determinant of deformation gradient
    double J = pt.m_J;
    
    // calculate deviatoric left Cauchy-Green tensor
    mat3ds B = pt.DevLeftCauchyGreen();
    
    // Invariant of B tilde (= invariants of C tilde)
    double I1 = B.tr();
    
    double sed = 0;
    for (int i=0; i<MAX_TERMS; ++i)
        sed += c[i]*pow(I1-3,i+1);

    return sed;
}
