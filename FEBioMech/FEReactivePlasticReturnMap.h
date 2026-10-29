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



#pragma once
#include <FECore/mat3d.h>

class FEMaterialPoint;
class FEElasticMaterial;
class FEDamageCriterion;

//-----------------------------------------------------------------------------
// Return mapping for a single bond family of the reactive plasticity materials
// (Zimmerman et al., J Mech Phys Solids 155:104534, 2021, Section 3.4).
// Given the total deformation gradient Fs, the inverse plastic deformation gradient
// (Fus)^-1 of the previous generation, and the yield threshold Ky, find (Fvs)^-1 such
// that the yield function phi = Kv - Ky vanishes. On entry the elastic material point
// must hold the total deformation gradient; it is restored on exit.
class FEReactivePlasticReturnMap
{
public:
    struct Result {
        mat3d   Fvsi;   // inverse of plastic deformation gradient of yielded bonds
        double  Kv;     // yield measure evaluated at the final elastic deformation
        double  phi;    // remaining yield function Kv - Ky
        bool    conv;   // true if |phi| <= rtol*Ky
    };

public:
    FEReactivePlasticReturnMap(FEElasticMaterial* pBase, FEDamageCriterion* pCrit, bool isochoric, double rtol);

    // original formulation: Uv = xi*sym(Ue*(I - lam*Nhat)), Fv = R*Uv, with R the rotation
    // of the total deformation gradient and Nhat = dPhi/dU at the elastic predictor
    Result Original(FEMaterialPoint& pt, const mat3d& Fs, const mat3d& Fusi, double Ky);

    // symmetric formulation: (Fvs)^-1 = xi*(Fus)^-1*G, G = I - lam*Nhat, Fv = xi*Fe*G,
    // without polar decomposition. If mandel is true, Nhat is normal to the yield surface
    // in Mandel stress space (classical associated flow); otherwise Nhat = dPhi/dG.
    Result Symmetric(FEMaterialPoint& pt, const mat3d& Fs, const mat3d& Fusi, double Ky, bool mandel);

    // get the tensor M such that dPhi = M : L, where L = dF*F^-1 (Eq. 3.27)
    mat3d YieldStressTangent(FEMaterialPoint& mp);

    // get the yield surface normal dPhi/dU (Eq. 3.26)
    mat3ds YieldSurfaceNormal(FEMaterialPoint& mp);

private:
    FEElasticMaterial*  m_pBase;    // base elastic material
    FEDamageCriterion*  m_pCrit;    // yield criterion
    bool                m_isochrc;  // flag for constraining plastic def grad to be isochoric
    double              m_rtol;     // relative tolerance on the yield function
};
