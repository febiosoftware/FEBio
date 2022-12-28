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
#include <limits>
#include "FEFiberNeoHookean.h"

//-----------------------------------------------------------------------------
// FEFiberNH
//-----------------------------------------------------------------------------

// define the material parameters
BEGIN_FECORE_CLASS(FEFiberNH, FEFiberMaterial)
	ADD_PARAMETER(m_mu, FE_RANGE_GREATER_OR_EQUAL(0.0), "mu")->setUnits(UNIT_PRESSURE);
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEFiberNH::FEFiberNH(FEModel* pfem) : FEFiberMaterial(pfem)
{ 
	m_mu = 0; 
	m_epsf = 0.0;
}

//-----------------------------------------------------------------------------
mat3ds FEFiberNH::FiberStress(FEMaterialPoint& mp, const vec3d& n0)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// deformation gradient
	mat3d &F = pt.m_F;
	double J = pt.m_J;
	
	// loop over all integration points
	mat3ds C = pt.RightCauchyGreen();
	
	// Calculate In = n0*C*n0
	double In_1 = n0*(C*n0) - 1.0;
	
	// only take fibers in tension into consideration
	mat3ds s;
	if (In_1 > 0.0)
	{
		// get the global spatial fiber direction in current configuration
		vec3d nt = F*n0;
		
		// calculate the outer product of nt
		mat3ds N = dyad(nt);
		
		// calculate the fiber stress
		s = N*(m_mu*In_1/J);
	}
	else
	{
		s.zero();
	}
	
	return s;
}

//-----------------------------------------------------------------------------
tens4ds FEFiberNH::FiberTangent(FEMaterialPoint& mp, const vec3d& n0)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// deformation gradient
	mat3d &F = pt.m_F;
	double J = pt.m_J;
	
	// loop over all integration points
	mat3ds C = pt.RightCauchyGreen();
	
	// Calculate In = n0*C*n0
	double In_1 = n0*(C*n0) - 1.0;
	
	// only take fibers in tension into consideration
	tens4ds c;
	const double eps = m_epsf*std::numeric_limits<double>::epsilon();
	if (In_1 > eps)
	{
		// get the global spatial fiber direction in current configuration
		vec3d nt = F*n0;
		
		// calculate the outer product of nt
		mat3ds N = dyad(nt);
		tens4ds NxN = dyad1s(N);
		
		// calculate the fiber tangent
		c = NxN*(2*m_mu/J);
	}
	else
	{
		c.zero();
	}
	
	return c;
}

//-----------------------------------------------------------------------------
//! Strain energy density
double FEFiberNH::FiberStrainEnergyDensity(FEMaterialPoint& mp, const vec3d& n0)
{
    
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// loop over all integration points
	mat3ds C = pt.RightCauchyGreen();
	
	// Calculate In = n0*C*n0
	double In_1 = n0*(C*n0) - 1.0;
	
	// only take fibers in tension into consideration
	double sed = 0.0;
	if (In_1 > 0.0)
        sed = 0.25*m_mu*In_1*In_1;
    
    return sed;
}

// define the material parameters
BEGIN_FECORE_CLASS(FEElasticFiberNH, FEElasticFiberMaterial)
	ADD_PARAMETER(m_fib.m_mu, FE_RANGE_GREATER_OR_EQUAL(0.0), "mu")->setUnits(UNIT_PRESSURE);
END_FECORE_CLASS();
