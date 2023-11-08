/*This file is part of the FEBio source code and is licensed under the MIT license
 listed below.
 
 See Copyright-FEBio.txt for details.
 
 Copyright (c) 2023 University of Utah, The Trustees of Columbia University in
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
#include "FEFiberCDFUncoupled.h"
#include "FEFiberCDFMaterialPoint.h"
#include <limits>
#include <FECore/log.h>

//-----------------------------------------------------------------------------
// FEFiberCDFUncoupled
//-----------------------------------------------------------------------------

// define the material parameters
BEGIN_FECORE_CLASS(FEFiberCDFUncoupled, FEFiberMaterialUncoupled)
    ADD_PARAMETER(m_E  , FE_RANGE_GREATER_OR_EQUAL(0.0), "E" )->setUnits(UNIT_PRESSURE)->setLongName("fiber modulus");
    ADD_PROPERTY (m_CDF,  "cdf")->SetLongName("cumulative distribution function");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEFiberCDFUncoupled::FEFiberCDFUncoupled(FEModel* pfem) : FEFiberMaterialUncoupled(pfem)
{
    m_E = 0;
    m_epsf = 1.0;
    m_CDF = nullptr;
}

//-----------------------------------------------------------------------------
// returns a pointer to a new material point object
FEMaterialPointData* FEFiberCDFUncoupled::CreateMaterialPointData()
{
    return new FEFiberCDFMaterialPoint(FEFiberMaterialUncoupled::CreateMaterialPointData());
}

//-----------------------------------------------------------------------------
mat3ds FEFiberCDFUncoupled::DevFiberStress(FEMaterialPoint& mp, const vec3d& n0)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // deformation gradient
    double J = pt.m_J;
    mat3d F = pt.m_F*pow(J, -1./3.);
    
    // loop over all integration points
    mat3ds C = pt.DevRightCauchyGreen();
    mat3ds s;
    
    // Calculate In - 1 = n0*C*n0 - 1
    double In_1 = n0*(C*n0) - 1;
    
    FEFiberCDFMaterialPoint& fp = *mp.ExtractData<FEFiberCDFMaterialPoint>();
    fp.SetFiberStrain(In_1);
    
    // only take fibers in tension into consideration
    const double eps = m_epsf* std::numeric_limits<double>::epsilon();
    if (In_1 >= eps)
    {
        double ksi = m_E(mp)/4;
        
        // get the global spatial fiber direction in current configuration
        vec3d nt = F*n0;
        
        // calculate the outer product of nt
        mat3ds N = dyad(nt);
        
        // perform integration
        Integrate(mp, In_1);

        // calculate strain energy first derivative
        double Wl = ksi*fp.m_ds_t*In_1;
        
        // calculate the fiber stress
        s = N*(2.0*Wl/J);
    }
    else
    {
        s.zero();
    }
    
    return s.dev();
}

//-----------------------------------------------------------------------------
tens4ds FEFiberCDFUncoupled::DevFiberTangent(FEMaterialPoint& mp, const vec3d& n0)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // deformation gradient
    double J = pt.m_J;
    mat3d F = pt.m_F*pow(J, -1./3.);
    
    // loop over all integration points
    mat3ds C = pt.DevRightCauchyGreen();
    tens4ds c;
    
    // Calculate In - 1 = n0*C*n0 - 1
    double In_1 = n0*(C*n0) - 1;
    
    FEFiberCDFMaterialPoint& fp = *mp.ExtractData<FEFiberCDFMaterialPoint>();
    fp.SetFiberStrain(In_1);
    
    // only take fibers in tension into consideration
    const double eps = m_epsf*std::numeric_limits<double>::epsilon();
    if (In_1 >= eps)
    {
        double ksi = m_E(mp)/4;
        
        // get the global spatial fiber direction in current configuration
        vec3d nt = F*n0;
        
        // calculate the outer product of nt
        mat3ds N = dyad(nt);
        tens4ds NxN = dyad1s(N);
        
        // perform integration
        Integrate(mp, In_1);

        // calculate strain energy 2nd derivative
        double Wll = ksi*fp.m_d2s_t;
        
        // calculate the fiber tangent
        c = NxN*(4.0*Wll/J);
        
        // calculate strain energy first derivative
        double Wl = ksi*fp.m_ds_t*In_1;
        
        // calculate the fiber stress
        mat3ds s = N*(2.0*Wl/J);
        
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
double FEFiberCDFUncoupled::DevFiberStrainEnergyDensity(FEMaterialPoint& mp, const vec3d& n0)
{
    double sed = 0.0;
    
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    mat3d F = pt.m_F*pow(pt.m_J, -1./3.);

    // loop over all integration points
    mat3ds C = pt.DevRightCauchyGreen();
    
    // Calculate In - 1 = n0*C*n0 - 1
    double In = n0*(C*n0);
    double In_1 = In - 1.0;
    
    FEFiberCDFMaterialPoint& fp = *mp.ExtractData<FEFiberCDFMaterialPoint>();
    fp.SetFiberStrain(In_1);
    
    // only take fibers in tension into consideration
    const double eps = 0;
    if (In_1 >= eps)
    {
        double ksi = m_E(mp)/4;
        
        // perform integration
        Integrate(mp, In_1);

        // calculate strain energy
        sed = fp.m_sed_t*ksi;
    }
    
    return sed;
}

//-----------------------------------------------------------------------------
void FEFiberCDFUncoupled::Integrate(FEMaterialPoint& mp, const double In_1)
{
    double dImax = 1e-2;
    FEFiberCDFMaterialPoint& fp = *mp.ExtractData<FEFiberCDFMaterialPoint>();
    
    // check increment in strain
    double dIn = In_1 - fp.m_In_1_p;
    
    // if small enough, integrate in a single step
    if (fabs(dIn) <= dImax) {
        fp.Integrate(m_CDF->cdf(mp,In_1));
        return;
    }
    // otherwise, integrate over multiple steps
    double In_1_p = fp.m_In_1_p;
    double d2s_p = fp.m_d2s_p;
    double ds_p = fp.m_ds_p;
    double sed_p = fp.m_sed_p;
    double In_1_t, d2s_t, ds_t, sed_t;
    do {
        In_1_t = In_1_p + dImax;
        if (In_1_t > In_1) In_1_t = In_1;
        dIn = In_1_t - In_1_p;
        d2s_t = m_CDF->cdf(mp,In_1_t);
        sed_t = sed_p + ds_p*dIn + 0.25*dIn*dIn*(d2s_t + d2s_p);
        ds_t = ds_p + 0.5*dIn*(d2s_t + d2s_p);
        fp.m_d2s_t = d2s_t;
        fp.m_ds_t = ds_t;
        fp.m_sed_t = sed_t;
        In_1_p = In_1_t;
        d2s_p = d2s_t;
        ds_p = ds_t;
        sed_p = sed_t;
    } while (In_1_t < In_1);
    return;
}

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FEElasticFiberCDFUncoupled, FEElasticFiberMaterialUC)
    ADD_PARAMETER(m_fib.m_E  , FE_RANGE_GREATER_OR_EQUAL(0.0), "E" )->setUnits(UNIT_PRESSURE)->setLongName("fiber modulus");
    ADD_PROPERTY (m_fib.m_CDF  ,  "cdf")->SetLongName("cumulative distribution function");
END_FECORE_CLASS();

