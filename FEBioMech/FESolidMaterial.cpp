/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, Columbia University, and others.

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
#include "FESolidMaterial.h"
#include "FEElasticMaterial.h"

// Material parameters for FEElasticMaterial
BEGIN_FECORE_CLASS(FESolidMaterial, FEMaterial)
	ADD_PARAMETER(m_density, "density");
END_FECORE_CLASS();

FESolidMaterial::FESolidMaterial(FEModel* pfem) : FEMaterial(pfem) {}

//! return the material density
FEParamDouble& FESolidMaterial::Density()
{ 
	return m_density; 
}

//! set the material density
void FESolidMaterial::SetDensity(const double d)
{ 
	m_density = d;
}

//! evaluate density
double FESolidMaterial::Density(FEMaterialPoint& pt)
{
	return m_density(pt);
}

//-----------------------------------------------------------------------------
//! calculate the 2nd Piola-Kirchhoff stress at material point, using prescribed Lagrange strain
//! needed for EAS analyses where the compatible strain (calculated from displacements) is enhanced
mat3ds FESolidMaterial::PK2Stress(FEMaterialPoint& mp, const mat3ds E)
{
    // Evaluate right Cauchy-Green tensor from E
    mat3ds C = mat3dd(1) + E*2;
    
    // Evaluate right stretch tensor U from C
    vec3d v[3];
    double lam[3];
    C.eigen2(lam, v);
    lam[0] = sqrt(lam[0]); lam[1] = sqrt(lam[1]); lam[2] = sqrt(lam[2]);
    mat3ds U = dyad(v[0])*lam[0] + dyad(v[1])*lam[1] + dyad(v[2])*lam[2];
    double J = lam[0]*lam[1]*lam[2];
    
    // temporarily replace F in material point with U
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    mat3d Fsafe = pt.m_F;
    double Jsafe = pt.m_J;
    pt.m_F = U;
    pt.m_J = J;
    
    // Evaluate Cauchy stress
    mat3ds s = Stress(mp);
    
    // Restore original F
    pt.m_F = Fsafe;
    pt.m_J = Jsafe;
    
    // Convert Cauchy stress to 2nd P-K stress
    mat3ds Ui = dyad(v[0])/lam[0] + dyad(v[1])/lam[1] + dyad(v[2])/lam[2];
    mat3ds S = Ui*s*Ui*J;
    
    return S;
}

//-----------------------------------------------------------------------------
//! calculate material tangent stiffness at material point, using prescribed Lagrange strain
//! needed for EAS analyses where the compatible strain (calculated from displacements) is enhanced
tens4ds FESolidMaterial::MaterialTangent(FEMaterialPoint& mp, const mat3ds E)
{
    // Evaluate right Cauchy-Green tensor from E
    mat3ds C = mat3dd(1) + E*2;
    
    // Evaluate right stretch tensor U from C
    vec3d v[3];
    double lam[3];
    C.eigen2(lam, v);
    lam[0] = sqrt(lam[0]); lam[1] = sqrt(lam[1]); lam[2] = sqrt(lam[2]);
    mat3d U = dyad(v[0])*lam[0] + dyad(v[1])*lam[1] + dyad(v[2])*lam[2];
    double J = lam[0]*lam[1]*lam[2];
    
    // temporarily replace F in material point with U
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    mat3d Fsafe = pt.m_F;
    double Jsafe = pt.m_J;
    pt.m_F = U;
    pt.m_J = J;
    
    // Evaluate Cauchy stress
    tens4ds c = Tangent(mp);
    
    // Restore original F
    pt.m_F = Fsafe;
    pt.m_J = Jsafe;
    
    // Convert spatial tangent to material tangent
    mat3d Ui = dyad(v[0])/lam[0] + dyad(v[1])/lam[1] + dyad(v[2])/lam[2];
    tens4ds Cm = c.pp(Ui)*J;
    
    return Cm;
}
