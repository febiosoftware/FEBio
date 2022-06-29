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
#include "FECustomFiberDistribution.h"
#include <FECore/sphericalHarmonics.h>
#include <FECore/spherePoints.h>

BEGIN_FECORE_CLASS(FECustomFiberDistribution, FEElasticMaterial)

	// material properties
	ADD_PROPERTY(m_pFmat, "fibers");
    ADD_PARAMETER(m_shpHar, "shp_harmonics");

	ADD_PROPERTY(m_Q, "mat_axis")->SetFlags(FEProperty::Optional);

END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FECustomFiberDistribution::FECustomFiberDistribution(FEModel* pfem) : FEElasticMaterial(pfem)
{
	m_pFmat = 0;
    m_ODF.reserve(NPTS);
}

//-----------------------------------------------------------------------------
FECustomFiberDistribution::~FECustomFiberDistribution() {}

//-----------------------------------------------------------------------------
FEMaterialPoint* FECustomFiberDistribution::CreateMaterialPointData()
{
	FEMaterialPoint* mp = FEElasticMaterial::CreateMaterialPointData();
	mp->SetNext(m_pFmat->CreateMaterialPointData());
    return mp;
}

//-----------------------------------------------------------------------------
bool FECustomFiberDistribution::Init()
{
    // initialize base class
	if (FEElasticMaterial::Init() == false) return false;

    // Get harmonic order by looking at number of coefficients
    int order = (sqrt(8*m_shpHar.size() + 1) - 3)/2;

    double* theta = new double[NPTS] {};
    double* phi = new double[NPTS] {};

    getSphereCoords(NPTS, XCOORDS, YCOORDS, ZCOORDS, theta, phi);

    auto T = compSH(order, NPTS, theta, phi);

    delete[] theta;
    delete[] phi;

    (*T).mult(m_shpHar, m_ODF);

    // Many of the ODF values end up being 0. This allows us to skip them
    // during the later evaluations.

    std::vector<int> indices;

    double mean = 0;
    for(int index = 0; index < NPTS; index++)
    {
        mean += m_ODF[index];
    }
    mean /= NPTS*100;

    for(int index = 0; index < NPTS; index++)
    {
        if(m_ODF[index] > mean)
        {
            indices.push_back(index);
        }
    }

    for(auto index : indices)
    {
        if(ZCOORDS[index] > 0)
        {
            m_posIndices.push_back(index);
        }
        else if(ZCOORDS[index] == 0)
        {
            m_zeroIndices.push_back(index);
        }
    }

	return true;
}

//-----------------------------------------------------------------------------
//! Serialization
void FECustomFiberDistribution::Serialize(DumpStream& ar)
{	
	FEElasticMaterial::Serialize(ar);
	if (ar.IsShallow()) return;
}

//-----------------------------------------------------------------------------
//! calculate stress at material point
mat3ds FECustomFiberDistribution::Stress(FEMaterialPoint& mp)
{ 
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    FEFiberMaterialPoint& fp = *mp.ExtractData<FEFiberMaterialPoint>();

	// calculate stress
	mat3ds s; s.zero();

	// get the local coordinate system
	mat3d Q = GetLocalCS(mp);
    
    for(int index : m_posIndices)
    {
        // get the fiber direction for that fiber distribution
        vec3d N(XCOORDS[index], YCOORDS[index], ZCOORDS[index]);

        // evaluate ellipsoidally distributed material coefficients
        double R = m_ODF[index]*2;
        
        // convert fiber to global coordinates
        vec3d n0 = Q*N;
        
        // calculate the stress
        s += m_pFmat->FiberStress(pt, fp.FiberPreStretch(n0))*(R);
	}

    for(int index : m_zeroIndices)
    {
        // get the fiber direction for that fiber distribution
        vec3d N(XCOORDS[index], YCOORDS[index], ZCOORDS[index]);

        // evaluate ellipsoidally distributed material coefficients
        double R = m_ODF[index];
        
        // convert fiber to global coordinates
        vec3d n0 = Q*N;
        
        // calculate the stress
        s += m_pFmat->FiberStress(pt, fp.FiberPreStretch(n0))*(R);
	}

	return s;
}

//-----------------------------------------------------------------------------
//! calculate tangent stiffness at material point
tens4ds FECustomFiberDistribution::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    FEFiberMaterialPoint& fp = *mp.ExtractData<FEFiberMaterialPoint>();

	// get the local coordinate system
	mat3d Q = GetLocalCS(mp);
    
	// initialize stress tensor
	tens4ds c;
	c.zero();
    for(int index : m_posIndices)
    {
        // get the fiber direction for that fiber distribution
        vec3d N(XCOORDS[index], YCOORDS[index], ZCOORDS[index]);
        
        // evaluate ellipsoidally distributed material coefficients
        double R = m_ODF[index]*2;
        
        // convert fiber to global coordinates
        vec3d n0 = Q*N;

        // calculate the tangent
        c += m_pFmat->FiberTangent(mp, fp.FiberPreStretch(n0))*(R);
	}

    for(int index : m_zeroIndices)
    {
        // get the fiber direction for that fiber distribution
        vec3d N(XCOORDS[index], YCOORDS[index], ZCOORDS[index]);
        
        // evaluate ellipsoidally distributed material coefficients
        double R = m_ODF[index];
        
        // convert fiber to global coordinates
        vec3d n0 = Q*N;

        // calculate the tangent
        c += m_pFmat->FiberTangent(mp, fp.FiberPreStretch(n0))*(R);
	}

	return c;
}

//-----------------------------------------------------------------------------
//! calculate strain energy density at material point
//double FEContinuousFiberDistribution::StrainEnergyDensity(FEMaterialPoint& pt) { return m_pFint->StrainEnergyDensity(pt); }
double FECustomFiberDistribution::StrainEnergyDensity(FEMaterialPoint& mp)
{ 
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    FEFiberMaterialPoint& fp = *mp.ExtractData<FEFiberMaterialPoint>();

	// get the local coordinate system
	mat3d Q = GetLocalCS(mp);
    
	double sed = 0.0;
    for(int index : m_posIndices)
    {
        // get the fiber direction for that fiber distribution
        vec3d N(XCOORDS[index], YCOORDS[index], ZCOORDS[index]);
        
        // evaluate ellipsoidally distributed material coefficients
        double R = m_ODF[index]*2;
        
        // convert fiber to global coordinates
        vec3d n0 = Q*N;

        // calculate the stress
        sed += m_pFmat->FiberStrainEnergyDensity(mp, fp.FiberPreStretch(n0))*(R);
    }

    for(int index : m_zeroIndices)
    {
        // get the fiber direction for that fiber distribution
        vec3d N(XCOORDS[index], YCOORDS[index], ZCOORDS[index]);
        
        // evaluate ellipsoidally distributed material coefficients
        double R = m_ODF[index];
        
        // convert fiber to global coordinates
        vec3d n0 = Q*N;

        // calculate the stress
        sed += m_pFmat->FiberStrainEnergyDensity(mp, fp.FiberPreStretch(n0))*(R);
    }

	return sed;
}