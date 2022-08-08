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
#include "FENeoHookean.h"

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FENeoHookean, FEElasticMaterial)
	ADD_PARAMETER(m_E, FE_RANGE_GREATER(0.0), "E")->setUnits(UNIT_PRESSURE)->setLongName("Young's modulus");
	ADD_PARAMETER(m_v, FE_RANGE_RIGHT_OPEN(-1, 0.5), "v")->setLongName("Poisson's ratio");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FENeoHookean::FENeoHookean(FEModel* pfem) : FEElasticMaterial(pfem) {}

//-----------------------------------------------------------------------------
mat3ds FENeoHookean::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	double detF = pt.m_J;
	double detFi = 1.0/detF;
	double lndetF = log(detF);

	// get the material parameters
	double E = m_E(mp);
	double v = m_v(mp);

	// calculate left Cauchy-Green tensor
	mat3ds b = pt.LeftCauchyGreen();

	// lame parameters
	double lam = v*E/((1+v)*(1-2*v));
	double mu  = 0.5*E/(1+v);

	// Identity
	mat3dd I(1);

	// calculate stress
	mat3ds s = (b - I)*(mu*detFi) + I*(lam*lndetF*detFi);

	return s;
}

//-----------------------------------------------------------------------------
tens4ds FENeoHookean::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// deformation gradient
	double detF = pt.m_J;

	// get the material parameters
	double E = m_E(mp);
	double v = m_v(mp);

	// lame parameters
	double lam = v*E/((1+v)*(1-2*v));
	double mu  = 0.5*E/(1+v);

	double lam1 = lam / detF;
	double mu1  = (mu - lam*log(detF)) / detF;
	
    mat3dd I(1);

	return dyad1s(I)*lam1 + dyad4s(I)*(2*mu1);
}

//-----------------------------------------------------------------------------
double FENeoHookean::StrainEnergyDensity(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	double J = pt.m_J;
	double lnJ = log(J);
	
	// get the material parameters
	double E = m_E(mp);
	double v = m_v(mp);

	// calculate left Cauchy-Green tensor
	mat3ds b = pt.LeftCauchyGreen();
	double I1 = b.tr();
	
	// lame parameters
	double lam = v*E/((1+v)*(1-2*v));
	double mu  = 0.5*E/(1+v);
	
	double sed = mu*((I1-3)/2.0 - lnJ)+lam*lnJ*lnJ/2.0;
	
	return sed;
}

//-----------------------------------------------------------------------------
mat3ds FENeoHookean::PK2Stress(FEMaterialPoint& pt, const mat3ds ES)
{
    // Identity
    mat3dd I(1);
    
    // calculate right Cauchy-Green tensor
    mat3ds C = I + ES*2;
    mat3ds Ci = C.inverse();
    
    double detF = sqrt(C.det());
    double lndetF = log(detF);
    
	// get the material parameters
	double E = m_E(pt);
	double v = m_v(pt);
    
    // lame parameters
    double lam = v*E/((1+v)*(1-2*v));
    double mu  = 0.5*E/(1+v);
    
    // calculate stress
    mat3ds S = (I - Ci)*mu + Ci*(lam*lndetF);
    
    return S;
}

//-----------------------------------------------------------------------------
tens4dmm FENeoHookean::MaterialTangent(FEMaterialPoint& pt, const mat3ds ES)
{
    // calculate right Cauchy-Green tensor
    mat3ds C = mat3dd(1) + ES*2;
    mat3ds Ci = C.inverse();
    double J = sqrt(C.det());

	// get the material parameters
	double E = m_E(pt);
	double v = m_v(pt);

    // lame parameters
    double lam = v*E/((1+v)*(1-2*v));
    double mu  = 0.5*E/(1+v);
    
    tens4dmm c = dyad1s(Ci)*lam + dyad4s(Ci)*(2*(mu-lam*log(J)));
    
    return c;
}
