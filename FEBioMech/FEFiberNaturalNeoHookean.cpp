/*This file is part of the FEBio source code and is licensed under the MIT license
 listed below.
 
 See Copyright-FEBio.txt for details.
 
 Copyright (c) 2020 University of Utah, The Trustees of Columbia University in
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
#include <limits>
#include "FEFiberNaturalNeoHookean.h"

//-----------------------------------------------------------------------------
// FEFiberNH
//-----------------------------------------------------------------------------

// define the material parameters
BEGIN_FECORE_CLASS(FEFiberNaturalNH, FEFiberMaterial)
    ADD_PARAMETER(m_ksi, FE_RANGE_GREATER_OR_EQUAL(0.0), "ksi");
    ADD_PARAMETER(m_lam0 , FE_RANGE_GREATER_OR_EQUAL(1.0), "lam0");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEFiberNaturalNH::FEFiberNaturalNH(FEModel* pfem) : FEFiberMaterial(pfem)
{
    m_ksi = 0;
    m_lam0 = 1;
    m_epsf = 1;
}

//-----------------------------------------------------------------------------
mat3ds FEFiberNaturalNH::FiberStress(FEMaterialPoint& mp, const vec3d& n0)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // deformation gradient
    mat3d &F = pt.m_F;
    double J = pt.m_J;
    
    double ksi = m_ksi(mp);
    double lam0 = m_lam0(mp);
    
    // get stretch ratio lamn along fiber
    vec3d nt = F*n0;
    double lamn = nt.unit();
    double epsn = log(lamn/lam0);
    const double epsf = m_epsf* std::numeric_limits<double>::epsilon();

    // only take fibers in tension into consideration
    mat3ds s;
    if (epsn > epsf)
    {
        // calculate the outer product of nt
        mat3ds N = dyad(nt);
        
        // calculate the fiber stress
        s = N*(ksi*epsn/J);
    }
    else
    {
        s.zero();
    }
    
    return s;
}

//-----------------------------------------------------------------------------
tens4ds FEFiberNaturalNH::FiberTangent(FEMaterialPoint& mp, const vec3d& n0)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // deformation gradient
    mat3d &F = pt.m_F;
    double J = pt.m_J;
    
    double ksi = m_ksi(mp);
    double lam0 = m_lam0(mp);
    
    // get stretch ratio lamn along fiber
    vec3d nt = F*n0;
    double lamn = nt.unit();
    double epsn = log(lamn/lam0);
    const double epsf = m_epsf* std::numeric_limits<double>::epsilon();

    // only take fibers in tension into consideration
    tens4ds c;
    if (epsn > epsf)
    {
        // calculate the outer product of nt
        mat3ds N = dyad(nt);
        tens4ds NxN = dyad1s(N);
        
        // calculate the fiber tangent
        c = NxN*(ksi*(1-2*epsn)/J);
    }
    else
    {
        c.zero();
    }
    
    return c;
}

//-----------------------------------------------------------------------------
//! Strain energy density
double FEFiberNaturalNH::FiberStrainEnergyDensity(FEMaterialPoint& mp, const vec3d& n0)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // deformation gradient
    mat3d &F = pt.m_F;
    
    double ksi = m_ksi(mp);
    double lam0 = m_lam0(mp);
    
    // get stretch ratio lamn along fiber
    vec3d nt = F*n0;
    double lamn = nt.unit();
    double epsn = log(lamn/lam0);
    const double epsf = m_epsf* std::numeric_limits<double>::epsilon();

    // only take fibers in tension into consideration
    double sed = 0.0;
    if (epsn > epsf)
        sed = 0.5*ksi*epsn*epsn;
    
    return sed;
}

// define the material parameters
BEGIN_FECORE_CLASS(FEElasticFiberNaturalNH, FEElasticFiberMaterial)
    ADD_PARAMETER(m_fib.m_ksi, FE_RANGE_GREATER_OR_EQUAL(0.0), "ksi");
    ADD_PARAMETER(m_fib.m_lam0 , FE_RANGE_GREATER_OR_EQUAL(1.0), "lam0");
END_FECORE_CLASS();
