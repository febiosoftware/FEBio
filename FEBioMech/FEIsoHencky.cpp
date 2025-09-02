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
#include "FEIsoHencky.h"

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FEIsoHencky, FEElasticMaterial)
	ADD_PARAMETER(m_E, FE_RANGE_GREATER(0.0), "E")->setUnits(UNIT_PRESSURE)->setLongName("Young's modulus");
	ADD_PARAMETER(m_v, FE_RANGE_RIGHT_OPEN(-1, 0.5), "v")->setLongName("Poisson's ratio");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEIsoHencky::FEIsoHencky(FEModel* pfem) : FEElasticMaterial(pfem) {}

//-----------------------------------------------------------------------------
mat3ds FEIsoHencky::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	double J = pt.m_J;
	double lnJ = log(J);

	// get the material parameters
	double E = m_E(mp);
	double v = m_v(mp);

	// calculate left Hencky tensor
	mat3ds h = pt.LeftHencky();

	// lame parameters
	double lam = v*E/((1+v)*(1-2*v));
	double mu  = 0.5*E/(1+v);

	// Identity
	mat3dd I(1);

	// calculate stress
	mat3ds s = h*(2*mu/J) + I*(lam*lnJ/J);

	return s;
}

//-----------------------------------------------------------------------------
tens4ds FEIsoHencky::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// deformation gradient
	double J = pt.m_J;

	// get the material parameters
	double E = m_E(mp);
	double v = m_v(mp);

	// lame parameters
	double lam = v*E/((1+v)*(1-2*v));
	double mu  = 0.5*E/(1+v);

	double lam1 = lam / J;
	double mu1  = mu / J;
	
    mat3dd I(1);

	return dyad1s(I)*lam1 + dyad4s(I)*(2*mu1);
}

//-----------------------------------------------------------------------------
double FEIsoHencky::StrainEnergyDensity(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	double J = pt.m_J;
	double lnJ = log(J);
	
	// get the material parameters
	double E = m_E(mp);
	double v = m_v(mp);

	// calculate right Hencky tensor
	mat3ds H = pt.RightHencky();
	double I1 = H.tr();
    double I2 = H.dotdot(H);
	
	// lame parameters
	double lam = v*E/((1+v)*(1-2*v));
	double mu  = 0.5*E/(1+v);
	
	double sed = (lam/2)*pow(I1,2) + mu*I2;
	
	return sed;
}
