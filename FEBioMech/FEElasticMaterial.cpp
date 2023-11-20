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
#include "FEElasticMaterial.h"
#include "FECore/FEModel.h"

//BEGIN_FECORE_CLASS(FEElasticMaterial, FESolidMaterial)
//END_FECORE_CLASS();

FEElasticMaterial::FEElasticMaterial(FEModel* pfem) : FESolidMaterial(pfem)
{ 
	m_density = 1;
	AddDomainParameter(new FEElasticStress());
}

//-----------------------------------------------------------------------------
FEElasticMaterial::~FEElasticMaterial()
{ 
	
}

//-----------------------------------------------------------------------------
FEMaterialPointData* FEElasticMaterial::CreateMaterialPointData()
{ 
	return new FEElasticMaterialPoint;
}

//-----------------------------------------------------------------------------
//! calculate spatial tangent stiffness at material point, using secant method
mat3ds FEElasticMaterial::SecantStress(FEMaterialPoint& mp, bool PK2)
{
	// extract the deformation gradient
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	mat3d F = pt.m_F;
	double J = pt.m_J;
	mat3ds E = pt.Strain();
	mat3dd I(1);
	mat3d FiT = F.transinv();

	// calculate the 2nd P-K stress at the current deformation gradient
	double W = StrainEnergyDensity(mp);

	// create deformation gradient increment
	double eps = 1e-9;
	vec3d e[3];
	e[0] = vec3d(1, 0, 0); e[1] = vec3d(0, 1, 0); e[2] = vec3d(0, 0, 1);
	mat3ds S(0.0);
	for (int k = 0; k < 3; ++k) {
		for (int l = k; l < 3; ++l) {
			// evaluate incremental stress
			mat3d dF = FiT * ((e[k] & e[l]))*(eps*0.5);
			mat3d F1 = F + dF;
			pt.m_F = F1;
			pt.m_J = pt.m_F.det();

			double dW = StrainEnergyDensity(mp) - W;

			// evaluate the secant modulus
			S(k, l) = 2.0 * dW / eps;
		}
	}

    // restore values
    pt.m_F = F;
    pt.m_J = J;
    
    if (PK2) return S;
    else {
        // push from material to spatial frame
        mat3ds s = pt.push_forward(S);
        
        // return secant stress
        return s;
    }
}

//-----------------------------------------------------------------------------
//! return the strain energy density
double FEElasticMaterial::StrainEnergyDensity(FEMaterialPoint& pt) { return 0; }

//-----------------------------------------------------------------------------
FEElasticStress::FEElasticStress() : FEDomainParameter("stress")
{

}

FEParamValue FEElasticStress::value(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();
	return ep.m_s;
}
