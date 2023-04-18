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
#include "FEFiberPowLinearUncoupled.h"

// define the material parameters
BEGIN_FECORE_CLASS(FEFiberPowLinearUC, FEFiberMaterialUncoupled)
	ADD_PARAMETER(m_E    , FE_RANGE_GREATER(0.0), "E"    )->setUnits(UNIT_PRESSURE)->setLongName("fiber modulus");
	ADD_PARAMETER(m_lam0 , FE_RANGE_GREATER(1.0), "lam0" )->setLongName("toe stretch ratio");
	ADD_PARAMETER(m_beta , FE_RANGE_GREATER_OR_EQUAL(2.0), "beta" )->setLongName("toe power exponent");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEFiberPowLinearUC::FEFiberPowLinearUC(FEModel* pfem) : FEFiberMaterialUncoupled(pfem)
{ 
    m_E = 0.0;
    m_lam0 = 1.0;
    m_beta = 2.0;
}

//-----------------------------------------------------------------------------
mat3ds FEFiberPowLinearUC::DevFiberStress(FEMaterialPoint& mp, const vec3d& n0)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // deformation gradient
    double J = pt.m_J;
    mat3d F = pt.m_F*pow(J, -1./3.);

    // loop over all integration points
    mat3ds C = pt.DevRightCauchyGreen();
    mat3ds s;
    
    // Calculate In
    double In = n0*(C*n0);
    
    // initialize material constants
    double E = m_E(mp);
    double lam0 = m_lam0(mp);
    double beta = m_beta(mp);
	double I0 = lam0*lam0;
	double ksi = E / 4.0 / (beta - 1)*pow(I0, -3.0 / 2.0)*pow(I0 - 1.0, 2.0 - beta);
	double b = ksi*pow(I0 - 1.0, beta - 1.0) + E / 2.0 / sqrt(I0);

    // only take fibers in tension into consideration
	const double eps = 0;
	if (In - 1 >= eps)
    {
        // get the global spatial fiber direction in current configuration
        vec3d nt = F*n0/sqrt(In);
        
        // calculate the outer product of nt
        mat3ds N = dyad(nt);
        
        // calculate the fiber stress magnitude
        double sn = (In < I0) ?
        2*In*ksi*pow(In-1, beta-1) :
        2*b*In - E*sqrt(In);
        
        // calculate the fiber stress
        s = N*(sn/J);
    }
    else
    {
        s.zero();
    }
    
    return s.dev();
}

//-----------------------------------------------------------------------------
tens4ds FEFiberPowLinearUC::DevFiberTangent(FEMaterialPoint& mp, const vec3d& n0)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // deformation gradient
    double J = pt.m_J;
    mat3d F = pt.m_F*pow(J, -1./3.);
    
    // loop over all integration points
    const double eps = 0;
    mat3ds C = pt.DevRightCauchyGreen();
    mat3ds s;
    tens4ds c;
    
    // Calculate In
    double In = n0*(C*n0);
    
    // only take fibers in tension into consideration
    if (In - 1 >= eps)
    {
        // get the global spatial fiber direction in current configuration
        vec3d nt = F*n0/sqrt(In);
        
        // calculate the outer product of nt
        mat3ds N = dyad(nt);
        tens4ds NxN = dyad1s(N);
        
        // initialize material constants
        double E = m_E(mp);
        double lam0 = m_lam0(mp);
        double beta = m_beta(mp);
        double I0 = lam0*lam0;
        double ksi = E / 4.0 / (beta - 1)*pow(I0, -3.0 / 2.0)*pow(I0 - 1.0, 2.0 - beta);
        double b = ksi*pow(I0 - 1.0, beta - 1.0) + E / 2.0 / sqrt(I0);

        // calculate the fiber stress magnitude
        double sn = (In < I0) ?
        2*In*ksi*pow(In-1, beta-1) :
        2*b*In - E*sqrt(In);
        
        // calculate the fiber stress
        s = N*(sn/J);
        
        // calculate modulus
        double cn = (In < I0) ? 4*In*In*ksi*(beta-1)*pow(In-1, beta-2) : E*sqrt(In);
        
        // calculate the fiber tangent
        c = NxN*(cn/J);

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
double FEFiberPowLinearUC::DevFiberStrainEnergyDensity(FEMaterialPoint& mp, const vec3d& n0)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // loop over all integration points
    const double eps = 0;
    mat3ds C = pt.DevRightCauchyGreen();
    
    // Calculate In = n0*C*n0
    double In = n0*(C*n0);
    
    // only take fibers in tension into consideration
	double sed = 0.0;
	if (In - 1 >= eps)
    {
        // initialize material constants
        double E = m_E(mp);
        double lam0 = m_lam0(mp);
        double beta = m_beta(mp);
        double I0 = lam0*lam0;
        double ksi = E / 4.0 / (beta - 1)*pow(I0, -3.0 / 2.0)*pow(I0 - 1.0, 2.0 - beta);
        double b = ksi*pow(I0 - 1.0, beta - 1.0) + E / 2.0 / sqrt(I0);

        // calculate strain energy density
        sed = (In < I0) ?
        ksi/beta*pow(In-1, beta) :
        b*(In-I0) - E*(sqrt(In) - sqrt(I0)) + ksi/beta*pow(I0-1, beta);
    }
    
    return sed;
}

// define the material parameters
BEGIN_FECORE_CLASS(FEUncoupledFiberPowLinear, FEElasticFiberMaterialUC)
	ADD_PARAMETER(m_fib.m_E    , FE_RANGE_GREATER(0.0), "E"    )->setUnits(UNIT_PRESSURE);
	ADD_PARAMETER(m_fib.m_lam0 , FE_RANGE_GREATER(1.0), "lam0" );
	ADD_PARAMETER(m_fib.m_beta , FE_RANGE_GREATER_OR_EQUAL(2.0), "beta" );
END_FECORE_CLASS();
