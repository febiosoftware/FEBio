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
#include "FEIsotropicElastic.h"

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FEIsotropicElastic, FEElasticMaterial)
	ADD_PARAMETER(m_E, FE_RANGE_GREATER(0.0), "E")->setUnits(UNIT_PRESSURE)->setLongName("Young's modulus");
	ADD_PARAMETER(m_v, FE_RANGE_RIGHT_OPEN(-1.0, 0.5), "v")->setLongName("Poisson's ratio");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
mat3ds FEIsotropicElastic::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	mat3d &F = pt.m_F;
	double Ji = 1.0 / pt.m_J;

	double E = m_E(mp);
	double v = m_v(mp);

	// lame parameters
	double lam = Ji*(v*E/((1+v)*(1-2*v)));
	double mu  = Ji*(0.5*E/(1+v));

	// calculate left Cauchy-Green tensor (ie. b-matrix)
	mat3ds b = pt.LeftCauchyGreen();

	// calculate trace of Green-Lagrance strain tensor
	double trE = 0.5*(b.tr()-3);

	// calculate square of b-matrix
	// (we commented out the matrix components we do not need)
	mat3ds b2 = b.sqr();

	// calculate stress
	mat3ds s = b*(lam*trE - mu) + b2*mu;

	return s;
}

//-----------------------------------------------------------------------------
tens4ds FEIsotropicElastic::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	double E = m_E(mp);
	double v = m_v(mp);

	// deformation gradient
	mat3d& F = pt.m_F;
	double Ji = 1.0 / pt.m_J;

	// lame parameters
	double lam = Ji*(v*E/((1+v)*(1-2*v)));
	double mu  = Ji*(0.5*E/(1+v));

	// left cauchy-green matrix (i.e. the 'b' matrix)
	mat3ds b = pt.LeftCauchyGreen();

	return dyad1s(b)*lam + dyad4s(b)*(2.0*mu);
}

//-----------------------------------------------------------------------------
double FEIsotropicElastic::StrainEnergyDensity(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	double mE = m_E(mp);
	double mv = m_v(mp);

    mat3ds E = (pt.RightCauchyGreen() - mat3dd(1))/2;
    
    double lam = mv*mE/((1+mv)*(1-2*mv));
	double mu  = 0.5*mE/(1+mv);

    double trE = E.tr();
    double Enorm = E.norm();
    
    double sed = lam*trE*trE/2 + mu*Enorm*Enorm;
    
	return sed;
}

//-----------------------------------------------------------------------------
mat3ds FEIsotropicElastic::PK2Stress(FEMaterialPoint& pt, const mat3ds E)
{
	double mE = m_E(pt);
	double mv = m_v(pt);

    // lame parameters
    double lam = (mv*mE/((1+mv)*(1-2*mv)));
    double mu  = (0.5*mE/(1+mv));
    
    // Identity
    mat3dd I(1);
    
    // calculate stress
    mat3ds S = I*(E.tr()*lam) + E*(2*mu);
    
    return S;
}

//-----------------------------------------------------------------------------
tens4dmm FEIsotropicElastic::MaterialTangent(FEMaterialPoint& pt, const mat3ds E)
{
	double mE = m_E(pt);
	double mv = m_v(pt);

    // lame parameters
    double lam = (mv*mE/((1+mv)*(1-2*mv)));
    double mu  = (0.5*mE/(1+mv));
    
    // Identity
    mat3dd I(1);
    
    tens4dmm c = dyad1s(I)*lam + dyad4s(I)*(2*mu);
    
    return c;
}
