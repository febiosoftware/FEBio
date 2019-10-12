/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, The Trustees of Columbia University in 
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
BEGIN_FECORE_CLASS(FEFiberPowLinearUncoupled, FEElasticFiberMaterialUC)
	ADD_PARAMETER(m_E    , FE_RANGE_GREATER(0.0), "E"    );
	ADD_PARAMETER(m_lam0 , FE_RANGE_GREATER(1.0), "lam0" );
	ADD_PARAMETER(m_beta , FE_RANGE_GREATER_OR_EQUAL(2.0), "beta" );
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEFiberPowLinearUncoupled::FEFiberPowLinearUncoupled(FEModel* pfem) : FEElasticFiberMaterialUC(pfem)
{ 

}

//-----------------------------------------------------------------------------
bool FEFiberPowLinearUncoupled::Validate()
{
	if (FEElasticFiberMaterialUC::Validate() == false) return false;

	// initialize material constants
	m_I0 = m_lam0*m_lam0;
	m_ksi = m_E / 4 / (m_beta - 1)*pow(m_I0, -3. / 2.)*pow(m_I0 - 1, 2 - m_beta);
	m_b = m_ksi*pow(m_I0 - 1, m_beta - 1) + m_E / 2 / sqrt(m_I0);

	return true;
}

//-----------------------------------------------------------------------------
mat3ds FEFiberPowLinearUncoupled::DevFiberStress(FEMaterialPoint& mp, const vec3d& n0)
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
    
    // only take fibers in tension into consideration
	const double eps = 0;
	if (In - 1 >= eps)
    {
        // get the global spatial fiber direction in current configuration
        vec3d nt = F*n0/sqrt(In);
        
        // calculate the outer product of nt
        mat3ds N = dyad(nt);
        
        // calculate the fiber stress magnitude
        double sn = (In < m_I0) ?
        2*In*m_ksi*pow(In-1, m_beta-1) :
        2*m_b*In - m_E/2;
        
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
tens4ds FEFiberPowLinearUncoupled::DevFiberTangent(FEMaterialPoint& mp, const vec3d& n0)
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
        
        // calculate the fiber stress magnitude
        double sn = (In < m_I0) ?
        2*In*m_ksi*pow(In-1, m_beta-1) :
        2*m_b*In - m_E/2;
        
        // calculate the fiber stress
        s = N*(sn/J);
        
        // calculate modulus
        double cn = (In < m_I0) ?
        4*In*In*m_ksi*(m_beta-1)*pow(In-1, m_beta-2) :
        m_E;
        
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
double FEFiberPowLinearUncoupled::DevFiberStrainEnergyDensity(FEMaterialPoint& mp, const vec3d& n0)
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
        // calculate strain energy density
        sed = (In < m_I0) ?
        m_ksi/m_beta*pow(In-1, m_beta) :
        m_b*(In-m_I0) - m_E/4*log(In/m_I0) + m_ksi/m_beta*pow(m_I0-1, m_beta);
    }
    
    return sed;
}
