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
#include "FEIncompNeoHookean.h"

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FEIncompNeoHookean, FEUncoupledMaterial)
	ADD_PARAMETER(m_G, FE_RANGE_GREATER(0.0), "G")->setUnits(UNIT_PRESSURE)->setLongName("shear modulus");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEIncompNeoHookean::FEIncompNeoHookean(FEModel* pfem) : FEUncoupledMaterial(pfem) 
{
	m_G = 0.0;
}

//-----------------------------------------------------------------------------
//! Calculate deviatoric stress
mat3ds FEIncompNeoHookean::DevStress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// shear modulus
	double G = m_G(mp);

	// deformation gradient
	double J = pt.m_J;

	// calculate left Cauchy-Green tensor
	mat3ds B = pt.LeftCauchyGreen();

	// calculate deviatoric stress
	return (B - mat3dd(B.tr()/3.))*(G/pow(J, 5.0/3.0));
}

//-----------------------------------------------------------------------------
//! Calculate deviatoric tangent
tens4ds FEIncompNeoHookean::DevTangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// shear modulus
	double G = m_G(mp);

	// deformation gradient
	double J = pt.m_J;

	// left cauchy-green matrix (i.e. the 'b' matrix)
	mat3ds B = pt.LeftCauchyGreen();

	// trace of b
	double Ib = B.tr();

	double muJ = G/pow(J, 5.0/3.0);

	mat3ds I(1,1,1,0,0,0);	// Identity

	tens4ds IxI = dyad1s(I);
	tens4ds I4  = dyad4s(I);
	tens4ds BxI = dyad1s(B, I); // = BxI + IxB

	return (I4*Ib -BxI + IxI*(Ib/3))*(2.0*muJ/3.0);
}

//-----------------------------------------------------------------------------
//! Calculate deviatoric strain energy density
double FEIncompNeoHookean::DevStrainEnergyDensity(FEMaterialPoint& mp)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // shear modulus
    double G = m_G(mp);
    
    mat3ds C = pt.DevRightCauchyGreen();
    double I1 = C.tr();
    
    // calculate deviatoric stress
    return (I1-3)*G/2;
}

