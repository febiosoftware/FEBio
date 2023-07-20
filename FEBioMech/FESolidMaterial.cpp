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
#include "FESolidMaterial.h"
#include "FEElasticMaterial.h"

// Material parameters for FEElasticMaterial
BEGIN_FECORE_CLASS(FESolidMaterial, FEMaterial)
	ADD_PARAMETER(m_density, "density")->setUnits(UNIT_DENSITY)->MakeTopLevel(true);
END_FECORE_CLASS();

FESolidMaterial::FESolidMaterial(FEModel* pfem) : FEMaterial(pfem)
{
    m_density = 1.0;
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

tens4dmm FESolidMaterial::SolidTangent(FEMaterialPoint& mp)
{
	return (UseSecantTangent() ? SecantTangent(mp) : Tangent(mp));
}

//-----------------------------------------------------------------------------
mat3ds FESolidMaterial::SecantStress(FEMaterialPoint& pt, bool PK2)
{
    assert(false);
    return mat3ds(0.0);
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
    mat3ds S = (Ui*s*Ui).sym()*J;
    
    return S;
}

//-----------------------------------------------------------------------------
//! calculate material tangent stiffness at material point, using prescribed Lagrange strain
//! needed for EAS analyses where the compatible strain (calculated from displacements) is enhanced
tens4dmm FESolidMaterial::MaterialTangent(FEMaterialPoint& mp, const mat3ds E)
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
    tens4dmm c = SolidTangent(mp);
    
    // Restore original F
    pt.m_F = Fsafe;
    pt.m_J = Jsafe;
    
    // Convert spatial tangent to material tangent
    mat3d Ui = dyad(v[0])/lam[0] + dyad(v[1])/lam[1] + dyad(v[2])/lam[2];
    tens4dmm Cm = c.pp(Ui)*J;
    
    return Cm;
}

//-----------------------------------------------------------------------------
//! calculate spatial tangent stiffness at material point, using secant method
tens4dmm FESolidMaterial::SecantTangent(FEMaterialPoint& mp, bool mat)
{
    // extract the deformation gradient
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    mat3d F = pt.m_F;
    double J = pt.m_J;
    mat3ds E = pt.Strain();
    mat3dd I(1);

    // calculate the 2nd P-K stress at the current deformation gradient
    mat3ds S = PK2Stress(mp,E);
    
    // create deformation gradient increment
    double eps = 1e-9;
    vec3d e[3];
    e[0] = vec3d(1,0,0); e[1] = vec3d(0,1,0); e[2] = vec3d(0,0,1);
    tens4dmm C;
    for (int k=0; k<3; ++k) {
        for (int l=k; l<3; ++l) {
            // evaluate incremental stress
            mat3ds dE = dyads(e[k], e[l])*(eps/2);
            mat3ds dS = (PK2Stress(mp,E+dE) - S)/eps;
            
            // evaluate the secant modulus
            C(0,0,k,l) = C(0,0,l,k) = dS.xx();
            C(1,1,k,l) = C(1,1,l,k) = dS.yy();
            C(2,2,k,l) = C(2,2,l,k) = dS.zz();
            C(0,1,k,l) = C(1,0,k,l) = C(0,1,l,k) = C(1,0,l,k) = dS.xy();
            C(1,2,k,l) = C(2,1,k,l) = C(1,2,l,k) = C(2,1,l,k) = dS.yz();
            C(2,0,k,l) = C(0,2,k,l) = C(2,0,l,k) = C(0,2,l,k) = dS.xz();
        }
    }
    
    if (mat) return C;
    else {
        
        // push from material to spatial frame
        tens4dmm c = C.pp(F)/J;
        
        // return secant tangent
        return c;
    }
}
