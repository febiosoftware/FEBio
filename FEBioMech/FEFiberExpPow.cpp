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
#include "FEFiberExpPow.h"
#include <limits>
#include <FECore/log.h>

// define the material parameters
BEGIN_FECORE_CLASS(FEFiberExpPow, FEFiberMaterial)
	ADD_PARAMETER(m_alpha, FE_RANGE_GREATER_OR_EQUAL(0.0), "alpha");
	ADD_PARAMETER(m_beta , FE_RANGE_GREATER_OR_EQUAL(2.0), "beta");
	ADD_PARAMETER(m_ksi  , FE_RANGE_GREATER_OR_EQUAL(0.0), "ksi" )->setUnits(UNIT_PRESSURE);
	ADD_PARAMETER(m_mu   , FE_RANGE_GREATER_OR_EQUAL(0.0), "mu"  )->setUnits(UNIT_PRESSURE);
    ADD_PARAMETER(m_lam0 , FE_RANGE_GREATER_OR_EQUAL(1.0), "lam0");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
// FEFiberExpPow
//-----------------------------------------------------------------------------

FEFiberExpPow::FEFiberExpPow(FEModel* pfem) : FEFiberMaterial(pfem)
{ 
	m_alpha = 0;
	m_beta = 2;
	m_ksi = 0;
	m_mu = 0;
    m_lam0 = 1;

	m_epsf = 1.0;
}

//-----------------------------------------------------------------------------
mat3ds FEFiberExpPow::FiberStress(FEMaterialPoint& mp, const vec3d& n0)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// deformation gradient
	mat3d &F = pt.m_F;
	double J = pt.m_J;
	
	// loop over all integration points
	mat3ds C = pt.RightCauchyGreen();
	mat3ds s;
	
	// Calculate In - I0 = n0*C*n0 - I0
    double lam0 = m_lam0(mp);
	double In_I0 = n0*(C*n0) - lam0*lam0;

	double ksi = m_ksi(mp);
	double mu = m_mu(mp);
    double alpha = m_alpha(mp);
    double beta = m_beta(mp);
	
	// only take fibers in tension into consideration
	const double eps = m_epsf* std::numeric_limits<double>::epsilon();
	if (In_I0 >= eps)
	{
		// get the global spatial fiber direction in current configuration
		vec3d nt = F*n0;
		
		// calculate the outer product of nt
		mat3ds N = dyad(nt);
		
		// calculate strain energy derivative
		double Wl = ksi*pow(In_I0, beta-1.0)*exp(alpha*pow(In_I0, beta));
		
		// calculate the fiber stress
		s = N*(2.0*Wl/J);

		// add the contribution from shear
		if (mu != 0.0)
		{
			mat3ds BmI = pt.LeftCauchyGreen() - mat3dd(1);
			s += (N*BmI).sym()*(m_mu(mp) / J);
		}
	}
	else
	{
		s.zero();
	}
	
	return s;
}

//-----------------------------------------------------------------------------
tens4ds FEFiberExpPow::FiberTangent(FEMaterialPoint& mp, const vec3d& n0)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// deformation gradient
	mat3d &F = pt.m_F;
	double J = pt.m_J;
	
	// loop over all integration points
	mat3ds C = pt.RightCauchyGreen();
	tens4ds c;
	
    // Calculate In - I0 = n0*C*n0 - I0
    double lam0 = m_lam0(mp);
    double In_I0 = n0*(C*n0) - lam0*lam0;

    double ksi = m_ksi(mp);
    double mu = m_mu(mp);
    double alpha = m_alpha(mp);
    double beta = m_beta(mp);

	// only take fibers in tension into consideration
	const double eps = m_epsf*std::numeric_limits<double>::epsilon();
	if (In_I0 >= eps)
	{
		// get the global spatial fiber direction in current configuration
		vec3d nt = F*n0;
		
		// calculate the outer product of nt
		mat3ds N = dyad(nt);
		tens4ds NxN = dyad1s(N);
		
		// calculate strain energy 2nd derivative
		double tmp = alpha*pow(In_I0, beta);
		double Wll = ksi*pow(In_I0, beta-2.0)*((tmp+1)*beta-1.0)*exp(tmp);
		
		// calculate the fiber tangent
		c = NxN*(4.0*Wll/J);

		// add the contribution from shear
		if (mu != 0.0)
		{
			mat3ds B = pt.LeftCauchyGreen();
			c += dyad4s(N, B)*(mu / J);
		}
	}
	else
	{
		c.zero();
	}
	
	return c;
}

//-----------------------------------------------------------------------------
double FEFiberExpPow::FiberStrainEnergyDensity(FEMaterialPoint& mp, const vec3d& n0)
{
    double sed = 0.0;
    
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// loop over all integration points
	mat3ds C = pt.RightCauchyGreen();
	
    // Calculate In - I0 = n0*C*n0 - I0
    double In = n0*(C*n0);
    double lam0 = m_lam0(mp);
    double In_I0 = In - lam0*lam0;
    
    double ksi = m_ksi(mp);
    double mu = m_mu(mp);
    double alpha = m_alpha(mp);
    double beta = m_beta(mp);

	// only take fibers in tension into consideration
	const double eps = 0;
	if (In_I0 >= eps)
	{
		// calculate strain energy derivative
        if (alpha > 0) {
            sed = ksi/(alpha*beta)*(exp(alpha*pow(In_I0, beta))-1);
        }
        else
            sed = ksi/beta*pow(In_I0, beta);

		// add the contribution from shear
		if (mu != 0.0)
		{
			mat3ds C2 = C.sqr();
			sed += mu * (n0*(C2*n0) - 2 * (In - 1) - 1) / 4.0;
		}

	}

    return sed;
}

// define the material parameters
BEGIN_FECORE_CLASS(FEElasticFiberExpPow, FEElasticFiberMaterial)
	ADD_PARAMETER(m_fib.m_alpha, FE_RANGE_GREATER_OR_EQUAL(0.0), "alpha");
	ADD_PARAMETER(m_fib.m_beta , FE_RANGE_GREATER_OR_EQUAL(2.0), "beta");
	ADD_PARAMETER(m_fib.m_ksi  , FE_RANGE_GREATER_OR_EQUAL(0.0), "ksi" )->setUnits(UNIT_PRESSURE);
	ADD_PARAMETER(m_fib.m_mu   , FE_RANGE_GREATER_OR_EQUAL(0.0), "mu"  )->setUnits(UNIT_PRESSURE);
    ADD_PARAMETER(m_fib.m_lam0 , FE_RANGE_GREATER_OR_EQUAL(1.0), "lam0");
END_FECORE_CLASS();


//-----------------------------------------------------------------------------
// FEFiberExponentialPower
//-----------------------------------------------------------------------------

// define the material parameters
BEGIN_FECORE_CLASS(FEFiberExponentialPower, FEElasticFiberMaterial)
	ADD_PARAMETER(m_alpha, FE_RANGE_GREATER_OR_EQUAL(0.0), "alpha");
	ADD_PARAMETER(m_beta , FE_RANGE_GREATER_OR_EQUAL(2.0), "beta" );
	ADD_PARAMETER(m_ksi  , FE_RANGE_GREATER_OR_EQUAL(0.0), "ksi"  )->setUnits(UNIT_PRESSURE);
    ADD_PARAMETER(m_mu   , FE_RANGE_GREATER_OR_EQUAL(0.0), "mu"   )->setUnits(UNIT_PRESSURE);
    ADD_PARAMETER(m_lam0 , FE_RANGE_GREATER_OR_EQUAL(1.0), "lam0");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEFiberExponentialPower::FEFiberExponentialPower(FEModel* pfem) : FEElasticFiberMaterial(pfem) 
{
	m_alpha = 0; 
	m_beta = 2; 
	m_ksi = 0; 
	m_mu = 0;
    m_lam0 = 1;

	m_epsf = 1.0;	// set to 1 for compatibility with febio 2.10
}

//-----------------------------------------------------------------------------
bool FEFiberExponentialPower::Validate()
{
	// TODO: how validate model parameters?
//	if ((4 * m_ksi + 2 * m_mu) < 0) {
//		feLogError("4*ksi+2*mu must be positive."); return false;
//	}
    return FEElasticFiberMaterial::Validate();
}

//-----------------------------------------------------------------------------
mat3ds FEFiberExponentialPower::FiberStress(FEMaterialPoint& mp, const vec3d& n0)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // deformation gradient
    mat3d &F = pt.m_F;
    double J = pt.m_J;
    
    // loop over all integration points
    mat3ds C = pt.RightCauchyGreen();
    mat3ds s;
    
    // Calculate In - I0 = n0*C*n0 - I0
    double lam0 = m_lam0(mp);
    double In_I0 = n0*(C*n0) - lam0*lam0;
    
    double ksi = m_ksi(mp);
    double mu = m_mu(mp);
    double alpha = m_alpha(mp);
    double beta = m_beta(mp);
    
    // only take fibers in tension into consideration
    const double eps = m_epsf* std::numeric_limits<double>::epsilon();
    if (In_I0 >= eps)
    {
        // get the global spatial fiber direction in current configuration
        vec3d nt = F*n0;
        
        // calculate the outer product of nt
        mat3ds N = dyad(nt);
        
        // calculate strain energy derivative
        double Wl = ksi*pow(In_I0, beta-1.0)*exp(alpha*pow(In_I0, beta));
        
        // calculate the fiber stress
        s = N*(2.0*Wl/J);
        
        // add the contribution from shear
        if (mu != 0.0)
        {
            mat3ds BmI = pt.LeftCauchyGreen() - mat3dd(1);
            s += (N*BmI).sym()*(m_mu(mp) / J);
        }
    }
    else
    {
        s.zero();
    }
    
    return s;
}

//-----------------------------------------------------------------------------
tens4ds FEFiberExponentialPower::FiberTangent(FEMaterialPoint& mp, const vec3d& n0)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // deformation gradient
    mat3d &F = pt.m_F;
    double J = pt.m_J;
    
    // loop over all integration points
    mat3ds C = pt.RightCauchyGreen();
    tens4ds c;
    
    // Calculate In - I0 = n0*C*n0 - I0
    double lam0 = m_lam0(mp);
    double In_I0 = n0*(C*n0) - lam0*lam0;

    double ksi = m_ksi(mp);
    double mu = m_mu(mp);
    double alpha = m_alpha(mp);
    double beta = m_beta(mp);
    
    // only take fibers in tension into consideration
    const double eps = m_epsf*std::numeric_limits<double>::epsilon();
    if (In_I0 >= eps)
    {
        // get the global spatial fiber direction in current configuration
        vec3d nt = F*n0;
        
        // calculate the outer product of nt
        mat3ds N = dyad(nt);
        tens4ds NxN = dyad1s(N);
        
        // calculate strain energy 2nd derivative
        double tmp = alpha*pow(In_I0, beta);
        double Wll = ksi*pow(In_I0, beta-2.0)*((tmp+1)*beta-1.0)*exp(tmp);
        
        // calculate the fiber tangent
        c = NxN*(4.0*Wll/J);
        
        // add the contribution from shear
        if (mu != 0.0)
        {
            mat3ds B = pt.LeftCauchyGreen();
            c += dyad4s(N, B)*(mu / J);
        }
    }
    else
    {
        c.zero();
    }
    
    return c;
}

//-----------------------------------------------------------------------------
//! Strain energy density
double FEFiberExponentialPower::FiberStrainEnergyDensity(FEMaterialPoint& mp, const vec3d& n0)
{
    double sed = 0.0;
    
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // loop over all integration points
    mat3ds C = pt.RightCauchyGreen();
    
    // Calculate In - I0 = n0*C*n0 - I0
    double In = n0*(C*n0);
    double lam0 = m_lam0(mp);
    double In_I0 = In - lam0*lam0;
    
    double ksi = m_ksi(mp);
    double mu = m_mu(mp);
    double alpha = m_alpha(mp);
    double beta = m_beta(mp);
    
    // only take fibers in tension into consideration
    const double eps = 0;
    if (In_I0 >= eps)
    {
        // calculate strain energy derivative
        if (alpha > 0) {
            sed = ksi/(alpha*beta)*(exp(alpha*pow(In_I0, beta))-1);
        }
        else
            sed = ksi/beta*pow(In_I0, beta);
        
        // add the contribution from shear
        if (mu != 0.0)
        {
            mat3ds C2 = C.sqr();
            sed += mu * (n0*(C2*n0) - 2 * (In - 1) - 1) / 4.0;
        }
        
    }
    
    return sed;
}
