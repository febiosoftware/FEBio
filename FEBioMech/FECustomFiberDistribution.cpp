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
#include <FEAMR/sphericalHarmonics.h>
#include <FEAMR/spherePoints.h>

BEGIN_FECORE_CLASS(FECustomFiberDistribution, FEElasticMaterial)

	// material properties
	ADD_PROPERTY(m_pFmat, "fibers");
    ADD_PARAMETER(m_shpHar, "shp_harmonics");

	ADD_PROPERTY(m_Q, "mat_axis")->SetFlags(FEProperty::Optional);

END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FECustomFiberDistribution::FECustomFiberDistribution(FEModel* pfem) 
    : FEElasticMaterial(pfem), m_lengthScale(10), m_hausd(0.05), m_grad(1.3)
{

	m_pFmat = 0;
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

    std::vector<double> gradient;
    altGradient(order, m_shpHar, gradient);
    
    vector<vec3i> elems;
    remesh(gradient, m_lengthScale, m_hausd, m_grad, m_nodePos, elems);

    int NN = m_nodePos.size();
    int NE = elems.size();

    double* xCoords = new double[NN] {};
    double* yCoords = new double[NN] {};
    double* zCoords = new double[NN] {};
    for(int index = 0; index < NN; index++)
    {
        vec3d vec = m_nodePos[index];

        xCoords[index] = vec.x;
        yCoords[index] = vec.y;
        zCoords[index] = vec.z;
    }

    double* theta = new double[NN] {};
    double* phi = new double[NN] {};

    getSphereCoords(NN, xCoords, yCoords, zCoords, theta, phi);

    auto T = compSH(order, NN, theta, phi);

    delete[] xCoords;
    delete[] yCoords;
    delete[] zCoords;
    delete[] theta;
    delete[] phi;

    m_ODF.resize(NN);

    (*T).mult(m_shpHar, m_ODF);

    vector<double> area(NN, 0);

    for(int index = 0; index < NE; index++)
    {
        int n0 = elems[index].x;
        int n1 = elems[index].y;
        int n2 = elems[index].z;

        vec3d n01 = m_nodePos[n0] - m_nodePos[n1];
        vec3d n02 = m_nodePos[n0] - m_nodePos[n2];

        double temp = (n01^n02).Length()/2;

        area[n0] += temp;
        area[n1] += temp;
        area[n2] += temp;
    }

    double sum = 0;
    for(int index = 0; index < NN; index++)
    {
        double val = m_ODF[index];

        if(m_nodePos[index].z != 0)
        {
            val *= 2;
        }

        val *= area[index];

        m_ODF[index] = val;

        sum += val;
    }

    for(int index = 0; index < NN; index++)
    {
        m_ODF[index] /= sum;
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

    for(int index = 0; index < m_ODF.size(); index++)
    {
        // evaluate ellipsoidally distributed material coefficients
        double R = m_ODF[index];
        
        // convert fiber to global coordinates
        vec3d n0 = Q*m_nodePos[index];
        
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
    for(int index = 0; index < m_ODF.size(); index++)
    {
        //evaluate ellipsoidally distributed material coefficients
        double R = m_ODF[index];
        
        // convert fiber to global coordinates
        vec3d n0 = Q*m_nodePos[index];

        // calculate the tangent
        c += m_pFmat->FiberTangent(mp, fp.FiberPreStretch(n0))*(R);
	}

	return c;
}

//-----------------------------------------------------------------------------
//! calculate strain energy density at material point
double FECustomFiberDistribution::StrainEnergyDensity(FEMaterialPoint& mp)
{ 
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    FEFiberMaterialPoint& fp = *mp.ExtractData<FEFiberMaterialPoint>();

	// get the local coordinate system
	mat3d Q = GetLocalCS(mp);

    double sed = 0.0;
    for(int index = 0; index < m_ODF.size(); index++)
    {
        // evaluate ellipsoidally distributed material coefficients
        double R = m_ODF[index];
        
        // convert fiber to global coordinates
        vec3d n0 = Q*m_nodePos[index];

        // calculate the stress
        sed += m_pFmat->FiberStrainEnergyDensity(mp, fp.FiberPreStretch(n0))*(R);
    }

	return sed;
}