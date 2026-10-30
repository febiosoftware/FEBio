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
#include "FEFiberCDF.h"
#include "FEFiberCDFMaterialPoint.h"
#include <limits>
#include <cmath>
#include <FECore/log.h>

//-----------------------------------------------------------------------------
// FEFiberCDF
//-----------------------------------------------------------------------------

// define the material parameters
BEGIN_FECORE_CLASS(FEFiberCDF, FEFiberMaterial)
	ADD_PARAMETER(m_E  , FE_RANGE_GREATER_OR_EQUAL(0.0), "E" )->setUnits(UNIT_PRESSURE)->setLongName("fiber modulus");
    ADD_PROPERTY (m_CDF,  "cdf")->SetLongName("cumulative distribution function");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEFiberCDF::FEFiberCDF(FEModel* pfem) : FEFiberMaterial(pfem)
{ 
	m_E = 0;
	m_epsf = 1.0;
    m_dImax = 0.05;
    m_CDF = nullptr;
}

//-----------------------------------------------------------------------------
// returns a pointer to a new material point object
FEMaterialPointData* FEFiberCDF::CreateMaterialPointData()
{
    return new FEFiberCDFMaterialPoint(FEFiberMaterial::CreateMaterialPointData());
}

//-----------------------------------------------------------------------------
mat3ds FEFiberCDF::FiberStress(FEMaterialPoint& mp, const vec3d& n0)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// deformation gradient
	mat3d &F = pt.m_F;
	double J = pt.m_J;
	
	// loop over all integration points
	mat3ds C = pt.RightCauchyGreen();
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
        double Wl = ksi*fp.m_ds_t;
        
		// calculate the fiber stress
		s = N*(2.0*Wl/J);
	}
	else
	{
		s.zero();
	}
	
	return s;
}

//-----------------------------------------------------------------------------
tens4ds FEFiberCDF::FiberTangent(FEMaterialPoint& mp, const vec3d& n0)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// deformation gradient
	mat3d &F = pt.m_F;
	double J = pt.m_J;
	
	// loop over all integration points
	mat3ds C = pt.RightCauchyGreen();
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
	}
	else
	{
		c.zero();
	}
	
	return c;
}

//-----------------------------------------------------------------------------
double FEFiberCDF::FiberStrainEnergyDensity(FEMaterialPoint& mp, const vec3d& n0)
{
    double sed = 0.0;
    
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// loop over all integration points
	mat3ds C = pt.RightCauchyGreen();
	
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
// Evaluate, with x = In - 1,
//   d2s = cdf(x)
//   ds  = int_0^x cdf(t) dt              (first derivative of sed w.r.t. In)
//   sed = int_0^x (x - t) cdf(t) dt      (= double integral of cdf, by Cauchy's formula)
// using composite 5-point Gauss-Legendre quadrature on [0,x], with segments
// no wider than m_dImax. The result depends only on the current value of x,
// so the response is hyperelastic (path-independent) and needs no history.
void FEFiberCDF::Integrate(FEMaterialPoint& mp, const double In_1)
{
    // 5-point Gauss-Legendre rule on [-1,1] (exact for polynomials of degree 9)
    const int NGP = 5;
    static const double xi[NGP] = {
        -0.906179845938663992797627, -0.538469310105683091036314, 0.0,
         0.538469310105683091036314,  0.906179845938663992797627 };
    static const double wi[NGP] = {
         0.236926885056189087514264,  0.478628670499366468041292, 0.568888888888888888888889,
         0.478628670499366468041292,  0.236926885056189087514264 };

    FEFiberCDFMaterialPoint& fp = *mp.ExtractData<FEFiberCDFMaterialPoint>();

    double ds = 0.0, sed = 0.0;
    if (In_1 > 0.0)
    {
        // number of segments and their half-width
        int nseg = (int)ceil(In_1/m_dImax);
        if (nseg < 1) nseg = 1;
        double hw = 0.5*In_1/nseg;

        for (int k=0; k<nseg; ++k)
        {
            double tc = (2*k + 1)*hw;   // segment midpoint
            for (int i=0; i<NGP; ++i)
            {
                double t = tc + hw*xi[i];
                double f = wi[i]*hw*m_CDF->cdf(mp, t);
                ds  += f;
                sed += f*(In_1 - t);
            }
        }
    }

    fp.m_d2s_t = (In_1 > 0.0) ? m_CDF->cdf(mp, In_1) : 0.0;
    fp.m_ds_t  = ds;
    fp.m_sed_t = sed;
}

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FEElasticFiberCDF, FEElasticFiberMaterial)
    ADD_PARAMETER(m_fib.m_E  , FE_RANGE_GREATER_OR_EQUAL(0.0), "E" )->setUnits(UNIT_PRESSURE)->setLongName("fiber modulus");
    ADD_PROPERTY (m_fib.m_CDF  ,  "cdf")->SetLongName("cumulative distribution function");
END_FECORE_CLASS();

