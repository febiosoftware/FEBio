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
#include "FEFiberIntegrationTrapezoidal.h"

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

class FEFiberIntegrationTrapezoidal::Iterator : public FEFiberIntegrationSchemeIterator
{
public:
	Iterator(FEMaterialPoint* mp, int nth)
	{
		m_nth = nth;

		a0 = vec3d(1,0,0);
		a1 = vec3d(0,1,0);

		double pi = 4 * atan(1.0);
		dth = pi / m_nth;  // integrate from 0 to pi

		i = -1;
		Next();
	}

	bool IsValid()
	{
		return (i < m_nth);
	}

	// move to the next integration point
	bool Next()
	{
		++i;
		if (i < m_nth)
		{
			double theta = i*dth;
			m_fiber = a0*cos(theta) + a1*sin(theta);

			// Multiply by 2 since fibers along theta+pi have same stress as along theta
			m_weight = dth*2.0;

			return true;
		}
		else return false;
	}

public:
	int	m_nth;
	double dth;
	vec3d a0, a1;

	int i;
};

//-----------------------------------------------------------------------------
// FEFiberIntegrationTrapezoidal
//-----------------------------------------------------------------------------

// register the material with the framework

// define the material parameters
BEGIN_FECORE_CLASS(FEFiberIntegrationTrapezoidal, FEFiberIntegrationScheme)
	ADD_PARAMETER(m_nth, FE_RANGE_GREATER(0), "nth");
END_FECORE_CLASS();

FEFiberIntegrationTrapezoidal::FEFiberIntegrationTrapezoidal(FEModel* pfem) : FEFiberIntegrationScheme(pfem)
{ 
	m_nth = 12; 
}

FEFiberIntegrationTrapezoidal::~FEFiberIntegrationTrapezoidal()
{
}

//-----------------------------------------------------------------------------
FEFiberIntegrationSchemeIterator* FEFiberIntegrationTrapezoidal::GetIterator(FEMaterialPoint* mp)
{
	return new Iterator(mp, m_nth);
}

/*
//-----------------------------------------------------------------------------
mat3ds FEFiberIntegrationTrapezoidal::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
    // initialize stress tensor
	mat3ds s;
	s.zero();
    
    double theta;
    double pi = 4*atan(1.0);
    double dth = pi/m_nth;  // integrate from 0 to pi
    
	// get the local coordinate systems
	mat3d Q = GetLocalCS(mp);
	vec3d a0(Q(0,0),Q(1,0),Q(2,0)); // local x-direction unit vector
    vec3d a1(Q(0,1),Q(1,1),Q(2,1)); // local y-direction unit vector
    
    vec3d n0e, n0a;
    
    // loop over all integration points
    for (int i=0; i<m_nth; ++i) {
        theta = i*dth;
        
        // set fiber direction in x-y plane of local coordinate system
        n0a = a0*cos(theta) + a1*sin(theta);
        // evaluate local fiber distribution
        double R = m_pFDD->FiberDensity(n0a);
        
        // rotate to global configuration to set fiber direction
        n0e = Q*n0a;
        m_pFmat->SetFiberDirection(mp, n0e);
        
        // calculate the stress
        s += m_pFmat->Stress(mp)*(R*dth);
    }
    
    // Multiply by 2 since fibers along theta+pi have same stress as along theta
	return s*2;
}

//-----------------------------------------------------------------------------
tens4ds FEFiberIntegrationTrapezoidal::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
    // initialize stress tensor
	tens4ds c;
	c.zero();
    
    double theta;
    double pi = 4*atan(1.0);
    double dth = pi/m_nth;  // integrate from 0 to pi
    
	// get the local coordinate systems
	mat3d Q = GetLocalCS(mp);
	vec3d a0(Q(0,0),Q(1,0),Q(2,0)); // local x-direction unit vector
    vec3d a1(Q(0,1),Q(1,1),Q(2,1)); // local y-direction unit vector
    
    vec3d n0e, n0a;
    
    // loop over all integration points
    for (int i=0; i<m_nth; ++i) {
        theta = i*dth;
        
        // set fiber direction in x-y plane of local coordinate system
        n0a = a0*cos(theta) + a1*sin(theta);
        // evaluate local fiber distribution
        double R = m_pFDD->FiberDensity(n0a);
        
        // rotate to global configuration to set fiber direction
        n0e = Q*n0a;
        m_pFmat->SetFiberDirection(mp, n0e);
        
        // calculate the stress
        c += m_pFmat->Tangent(mp)*(R*dth);
    }

    // Multiply by 2 since fibers along theta+pi have same stress as along theta
    return c*2;
}

//-----------------------------------------------------------------------------
double FEFiberIntegrationTrapezoidal::StrainEnergyDensity(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
    // initialize strain energy density
	double sed = 0.0;
    
    double theta;
    double pi = 4*atan(1.0);
    double dth = pi/m_nth;  // integrate from 0 to pi
    
	// get the local coordinate systems
	mat3d Q = GetLocalCS(mp);
	vec3d a0(Q(0,0),Q(1,0),Q(2,0)); // local x-direction unit vector
    vec3d a1(Q(0,1),Q(1,1),Q(2,1)); // local y-direction unit vector
    
    vec3d n0e, n0a;
    
    // loop over all integration points
    for (int i=0; i<m_nth; ++i) {
        theta = i*dth;
        
        // set fiber direction in x-y plane of local coordinate system
        n0a = a0*cos(theta) + a1*sin(theta);
        // evaluate local fiber distribution
        double R = m_pFDD->FiberDensity(n0a);
        
        // rotate to global configuration to set fiber direction
        n0e = Q*n0a;
        m_pFmat->SetFiberDirection(mp, n0e);
        
        // calculate the stress
        sed += m_pFmat->StrainEnergyDensity(mp)*(R*dth);
    }
    
    // Multiply by 2 since fibers along theta+pi have same sed as along theta
	return sed*2;
}

//-----------------------------------------------------------------------------
double FEFiberIntegrationTrapezoidal::IntegratedFiberDensity()
{
    // initialize integrated fiber density distribution
    double IFD = 1;
    
	double C = 0;
    
    double theta;
    double pi = 4*atan(1.0);
    double dth = pi/m_nth;  // integrate from 0 to pi
    
    vec3d a0(1,0,0); // local x-direction unit vector
    vec3d a1(0,1,0); // local y-direction unit vector
    
    vec3d n0a;
    
    // loop over all integration points
    for (int i=0; i<m_nth; ++i) {
        theta = i*dth;
        
        // set fiber direction in x-y plane of local coordinate system
        n0a = a0*cos(theta) + a1*sin(theta);
        // evaluate local fiber distribution
        double R = m_pFDD->FiberDensity(n0a);
        
        // integrate the fiber distribution
        C += R*dth;
    }
    
    // Multiply by 2 to take advantange of symmetry
	IFD = C*2;
    return IFD;
}
*/
