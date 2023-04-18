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
#include "FEFiberKiousisUncoupled.h"
#include <limits>
#include <FECore/log.h>

// define the material parameters
BEGIN_FECORE_CLASS(FEUncoupledFiberKiousis, FEElasticFiberMaterialUC)
    ADD_PARAMETER(m_d1, "d1");
    ADD_PARAMETER(m_d2, "d2");
    ADD_PARAMETER(m_n , "n" );
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEUncoupledFiberKiousis::FEUncoupledFiberKiousis(FEModel* pfem) : FEElasticFiberMaterialUC(pfem)
{
    m_d1 = 0;
    m_d2 = 1;
    m_n = 2;
}

//-----------------------------------------------------------------------------
mat3ds FEUncoupledFiberKiousis::DevFiberStress(FEMaterialPoint& mp, const vec3d& n0)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // deformation gradient
    double J = pt.m_J;
    mat3d F = pt.m_F*pow(J,-1.0/3.0);
    
    // loop over all integration points
    const double eps = 0;
    mat3ds C = pt.DevRightCauchyGreen();
    mat3ds s;
    
    double d1 = m_d1(mp);
    double d2 = m_d2(mp);
    double n = m_n(mp);
    
    // Calculate In = n0*C*n0
    double In_d2 = n0*(C*n0) - d2;

    // only take fibers in tension into consideration
    if (In_d2 >= eps)
    {
        // get the global spatial fiber direction in current configuration
        vec3d nt = F*n0;
        
        // calculate the outer product of nt
        mat3ds N = dyad(nt);
        
        // calculate the fiber stress
        s = N*(2.0*d1*pow(In_d2, n-1)/J);
    }
    else
    {
        s.zero();
    }
    
    return s.dev();
}

//-----------------------------------------------------------------------------
tens4ds FEUncoupledFiberKiousis::DevFiberTangent(FEMaterialPoint& mp, const vec3d& n0)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // deformation gradient
    double J = pt.m_J;
    mat3d F = pt.m_F*pow(J,-1.0/3.0);
    
    // loop over all integration points
    const double eps = 0;
    mat3ds C = pt.DevRightCauchyGreen();
    mat3ds s;
    tens4ds c;
    
    // Calculate In = n0*C*n0
    double In_d2 = n0*(C*n0) - m_d2(mp);
    double n = m_n(mp);

    // only take fibers in tension into consideration
    if (In_d2 >= eps)
    {
        // get the global spatial fiber direction in current configuration
        vec3d nt = F*n0;
        
        // calculate the outer product of nt
        mat3ds N = dyad(nt);
        tens4ds NxN = dyad1s(N);
        
        // calculate the fiber stress
        s = N*(2.0*m_d1(mp)*pow(In_d2, n-1)/J);
        
        // calculate the fiber tangent
        c = NxN*(4.0*(n-1)*m_d1(mp)*pow(In_d2, n-2)/J);
        
        // This is the final value of the elasticity tensor
        mat3dd I(1);
        tens4ds IxI = dyad1s(I);
        tens4ds I4  = dyad4s(I);
        c += ((I4+IxI/3.0)*s.tr() - dyad1s(I,s))*(2./3.)
        - (ddots(IxI, c)-IxI*(c.tr()/3.))/3.;
    }
    else
    {
        c.zero();
    }
    
    return c;
}

//-----------------------------------------------------------------------------
double FEUncoupledFiberKiousis::DevFiberStrainEnergyDensity(FEMaterialPoint& mp, const vec3d& n0)
{
    double sed = 0.0;
    
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // loop over all integration points
    mat3ds C = pt.DevRightCauchyGreen();
    
    // Calculate In = n0*C*n0
    double In_d2 = n0*(C*n0) - m_d2(mp);
    double n = m_n(mp);
    
    // only take fibers in tension into consideration
    const double eps = 0;
    if (In_d2 >= eps)
    {
        // calculate strain energy derivative
        sed = m_d1(mp)/n*pow(In_d2, n);
    }
    
    return sed;
}
