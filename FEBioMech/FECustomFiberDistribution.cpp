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
#include <cmath>
#include <algorithm>
#include <FECore/FEModel.h>
#include <FECore/FEDomain.h>
#include <FECore/FENode.h>
#include <FEAMR/sphericalHarmonics.h>
#include <FEAMR/spherePoints.h>

//=============================================================================
BEGIN_FECORE_CLASS(FEFiberODF, FECoreClass)
	ADD_PARAMETER(m_shpHar, "odf");
	ADD_PARAMETER(m_pos, "position");
END_FECORE_CLASS();

//=============================================================================

void FEElementODF::calcODF(std::vector<std::vector<double>>& ODFs)
{
    double maxWeight = std::distance(m_weights.begin(), std::max_element(m_weights.begin(), m_weights.end()));
    m_ODF = ODFs[maxWeight];

    // define gradient descent step size
    double eps = 0.25;

    // threshold for convegence 
    double thresh = 5e-8; 

    double prevPhi = 3;
    double nPhi = 2;

    std::vector<std::vector<double>> logmaps(ODFs.size(), std::vector<double>(NPTS, 0));

    while(nPhi < prevPhi)
    {
        // compute logmaps from current meanPODF to every other ODF
        for(int index = 0; index < m_ODF.size(); index++)
        {
            auto& currentLogmap = logmaps[index];
            auto& currentODF = ODFs[index];
            double currentWeight = m_weights[index];

            double dot = 0;
            for(int index2 = 0; index2 < NPTS; index2++)
            {
                dot += m_ODF[index2]*currentODF[index2];
            }

            // if the two vectors are the same, tangent is the 0 vector
            if(abs(1-dot) < 10e-12)
            {
                logmaps[index] = std::vector<double>(NPTS, 0);
            }
            else
            {
                double denom = sqrt(1-dot*dot)*acos(dot);

                for(int index2 = 0; index2 < NPTS; index2++)
                {
                    currentLogmap[index2] = (currentODF[index2] - dot*m_ODF[index2])/denom;
                }
            }
        }

        // compute tangent vector of each 
        std::vector<double> phi(m_ODF.size(), 0);
        for(int index = 0; index < logmaps.size(); index++)
        {
            for(int index2 = 0; index2 < NPTS; index2++)
            {
                phi[index2] += logmaps[index][index2]*m_weights[index];
            }

        }

        prevPhi = nPhi;

        // take norm of phi
        nPhi = 0;
        for(double val : phi)
        {
            nPhi += val*val;
        }
        nPhi = sqrt(nPhi);

        if(nPhi < thresh) break;

        // apply exponential map to mean ODF in direction of tanget vector
        for(int index = 0; index < phi.size(); index++)
        {
            phi[index] = phi[index]*eps;
        }

        double normPhi = 0;
        for(int index = 0; index < NPTS; index++)
        {
            normPhi += phi[index]*phi[index];
        }
        normPhi = sqrt(normPhi);

        if(normPhi >= 10e-12)
        {
            for(int index = 0; index < NPTS; index++)
            {
                m_ODF[index] = cos(normPhi)*m_ODF[index] + sin(normPhi)*phi[index]/normPhi;
            } 
        }
    }

    //  undo sqaure root transform to obtain mean ODF
    for(int index = 0; index < NPTS; index++)
    {
        m_ODF[index] = m_ODF[index]*m_ODF[index];
    }
}

//=============================================================================
BEGIN_FECORE_CLASS(FECustomFiberDistribution, FEElasticMaterial)

	// material properties
	ADD_PROPERTY(m_pFmat, "fibers");
    ADD_PROPERTY(m_ODF, "shp_harmonics");

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

    // Calculate the ODFs for each of the image subregions
    for(auto& ODF : m_ODF)
    {
        reconstructODF(ODF->m_shpHar, ODF->m_ODF);
    }

    FEMesh& mesh = GetFEModel()->GetMesh();
    
    auto& domainList = GetDomainList();
    for(int index = 0; index < domainList.Domains(); index++)
    {
        FEDomain* domain = domainList.GetDomain(index);

        for(int el = 0; el < domain->Elements(); el++)
        {
            FEElement& element = domain->ElementRef(el);

            FEElementODF* odf = new FEElementODF(NPTS, m_ODF.size());

            // Calculate the centroid of the element
            for(int node = 0; node < element.Nodes(); node++)
            {
                odf->m_pos += mesh.Node(element.m_node[node]).m_r0;
            }
            odf->m_pos/element.Nodes();

            // Calculate weight to each ODF
            double sum = 0;
            for(int index = 0; index < m_ODF.size(); index++)
            {
                auto current = m_ODF[index];

                double weight = 1/(odf->m_pos - current->m_pos).Length();
                odf->m_weights[index] = weight;
                sum += weight;
            }
            // Normalize weights
            for(int index = 0; index < m_ODF.size(); index++)
            {
                odf->m_weights[index] /= sum;
            }

            // Add element odf object to map
            m_ElODF[element.GetID()] = odf;
        }
    }

    // Apply the square root transform to each ODF
    std::vector<std::vector<double>> pODF;
    pODF.resize(m_ODF.size());

    for(int index = 0; index < m_ODF.size(); index++)
    {
        auto& currentODF = m_ODF[index]->m_ODF;
        auto& currentPODF = pODF[index];

        currentPODF.resize(NPTS);

        for(int index2 = 0; index2 < NPTS; index2++)
        {
            currentPODF[index2] = sqrt(currentODF[index2]);
        }
    }

    // Compute logmaps from each ODF to every other ODF
    std::vector<std::vector<std::vector<double>>> logmaps;
    logmaps.resize(m_ODF.size());
    
    for(int index = 0; index < m_ODF.size(); index++)
    {
        auto& currentLogmapSet = logmaps[index];
        currentLogmapSet.resize(m_ODF.size());
        auto& meanODF = pODF[index];

        for(int index2 = 0; index2 < m_ODF.size(); index2++)
        {
            auto& currentLogmap = currentLogmapSet[index2];
            currentLogmap.resize(NPTS);

            auto& currentODF = pODF[index2];

            // if the two vectors are the same, tangent is the 0 vector
            if(index == index2)
            {
                currentLogmapSet[index2] = std::vector<double>(NPTS,0);
            }
            else
            {
                double dot = 0;
                for(int index3 = 0; index3 < NPTS; index3++)
                {
                    dot += meanODF[index3]*currentODF[index3];
                }

                double denom = sqrt(1-dot*dot)*acos(dot);

                for(int index3 = 0; index3 < NPTS; index3++)
                {
                    currentLogmap[index3] = (currentODF[index3] - dot*meanODF[index3])/denom;
                }
            }

        }       
    }

    //     // Get harmonic order by looking at number of coefficients
    //     int order = (sqrt(8*ODF->m_shpHar.size() + 1) - 3)/2;

    //     std::vector<double> gradient;
    //     altGradient(order, ODF->m_shpHar, gradient);
        
    //     vector<vec3i> elems;
    //     remesh(gradient, m_lengthScale, m_hausd, m_grad, ODF->m_nodePos, elems);

    //     int NN = ODF->m_nodePos.size();
    //     int NE = elems.size();

    //     double* xCoords = new double[NN] {};
    //     double* yCoords = new double[NN] {};
    //     double* zCoords = new double[NN] {};
    //     for(int index = 0; index < NN; index++)
    //     {
    //         vec3d vec = ODF->m_nodePos[index];

    //         xCoords[index] = vec.x;
    //         yCoords[index] = vec.y;
    //         zCoords[index] = vec.z;
    //     }

    //     double* theta = new double[NN] {};
    //     double* phi = new double[NN] {};

    //     getSphereCoords(NN, xCoords, yCoords, zCoords, theta, phi);

    //     auto T = compSH(order, NN, theta, phi);

    //     delete[] xCoords;
    //     delete[] yCoords;
    //     delete[] zCoords;
    //     delete[] theta;
    //     delete[] phi;

    //     ODF->m_ODF.resize(NN);

    //     (*T).mult(ODF->m_shpHar, ODF->m_ODF);

    //     vector<double> area(NN, 0);

    //     for(int index = 0; index < NE; index++)
    //     {
    //         int n0 = elems[index].x;
    //         int n1 = elems[index].y;
    //         int n2 = elems[index].z;

    //         vec3d n01 = ODF->m_nodePos[n0] - ODF->m_nodePos[n1];
    //         vec3d n02 = ODF->m_nodePos[n0] - ODF->m_nodePos[n2];

    //         double temp = (n01^n02).Length()/2;

    //         area[n0] += temp;
    //         area[n1] += temp;
    //         area[n2] += temp;
    //     }

    //     double sum = 0;
    //     for(int index = 0; index < NN; index++)
    //     {
    //         double val = ODF->m_ODF[index];

    //         if(ODF->m_nodePos[index].z != 0)
    //         {
    //             val *= 2;
    //         }

    //         val *= area[index];

    //         ODF->m_ODF[index] = val;

    //         sum += val;
    //     }

    //     for(int index = 0; index < NN; index++)
    //     {
    //         ODF->m_ODF[index] /= sum;
    //     }
    // }

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
        // // evaluate ellipsoidally distributed material coefficients
        // double R = m_ODF[index];
        
        // // convert fiber to global coordinates
        // vec3d n0 = Q*m_nodePos[index];
        
        // // calculate the stress
        // s += m_pFmat->FiberStress(pt, fp.FiberPreStretch(n0))*(R);
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
        // //evaluate ellipsoidally distributed material coefficients
        // double R = m_ODF[index];
        
        // // convert fiber to global coordinates
        // vec3d n0 = Q*m_nodePos[index];

        // // calculate the tangent
        // c += m_pFmat->FiberTangent(mp, fp.FiberPreStretch(n0))*(R);
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
        // // evaluate ellipsoidally distributed material coefficients
        // double R = m_ODF[index];
        
        // // convert fiber to global coordinates
        // vec3d n0 = Q*m_nodePos[index];

        // // calculate the stress
        // sed += m_pFmat->FiberStrainEnergyDensity(mp, fp.FiberPreStretch(n0))*(R);
    }

	return sed;
}