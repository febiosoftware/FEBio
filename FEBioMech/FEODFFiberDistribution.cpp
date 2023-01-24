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
#include "FEODFFiberDistribution.h"
#include <cmath>
#include <algorithm>
#include <FECore/FEModel.h>
#include <FECore/FEDomain.h>
#include <FECore/FENode.h>
#include <FEAMR/sphericalHarmonics.h>
#include <FEAMR/spherePoints.h>

#include <FECore/Timer.h>
#include <iostream>

//=============================================================================
BEGIN_FECORE_CLASS(FEFiberODF, FECoreClass)
	ADD_PARAMETER(m_shpHar, "shp_harmonics");
	ADD_PARAMETER(m_pos, "position");
END_FECORE_CLASS();

//=============================================================================

void FEElementODF::calcODF(std::vector<std::vector<double>>& ODFs)
{
    m_ODF.resize(NPTS);

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
        for(int index = 0; index < ODFs.size(); index++)
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
BEGIN_FECORE_CLASS(FEODFFiberDistribution, FEElasticMaterial)

	// material properties
	ADD_PROPERTY(m_pFmat, "fibers");
    ADD_PROPERTY(m_ODF, "fiber-odf");

	ADD_PROPERTY(m_Q, "mat_axis")->SetFlags(FEProperty::Optional);

END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEODFFiberDistribution::FEODFFiberDistribution(FEModel* pfem) 
    : FEElasticMaterial(pfem), m_lengthScale(10), m_hausd(0.05), m_grad(1.3), m_interpolate(false),
        m_theta(nullptr), m_phi(nullptr)
{
	m_pFmat = 0;
}

//-----------------------------------------------------------------------------
FEODFFiberDistribution::~FEODFFiberDistribution()
{
    if(m_theta) delete[] m_theta;
    if(m_phi) delete[] m_phi;
}

//-----------------------------------------------------------------------------
FEMaterialPointData* FEODFFiberDistribution::CreateMaterialPointData()
{
	FEMaterialPointData* mp = FEElasticMaterial::CreateMaterialPointData();
	mp->SetNext(m_pFmat->CreateMaterialPointData());
    return mp;
}

//-----------------------------------------------------------------------------
bool FEODFFiberDistribution::Init()
{
    Timer totalTime, initTime, interpTime, reduceTime;

    double remeshTime = 0;

    totalTime.start();
    initTime.start();

    // initialize base class
	if (FEElasticMaterial::Init() == false) return false;

    m_interpolate = m_ODF.size() > 1;

    // Get harmonic order by looking at number of coefficients of first ODF
    m_order = (sqrt(8*m_ODF[0]->m_shpHar.size() + 1) - 3)/2;

    // Initialize spherical coordinates and T matrix
    m_theta = new double[NPTS];
    m_phi = new double[NPTS];
    getSphereCoords(NPTS, XCOORDS, YCOORDS, ZCOORDS, m_theta, m_phi);
    m_T = compSH(m_order, NPTS, m_theta, m_phi);
    matrix transposeT = m_T->transpose();
    matrix B = (transposeT*(*m_T)).inverse()*transposeT;

    // Calculate the ODFs for each of the image subregions
    #pragma omp parallel for
    for(int index = 0; index < m_ODF.size(); index++)
    {
        FEFiberODF* ODF = m_ODF[index];
        reconstructODF(ODF->m_shpHar, ODF->m_ODF, NPTS, m_theta, m_phi);
    }

    initTime.stop();

    if(m_interpolate)
    {
        // Apply the square root transform to each ODF
        std::vector<std::vector<double>> pODF;
        pODF.resize(m_ODF.size());

        #pragma omp parallel for
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

        FEMesh& mesh = GetFEModel()->GetMesh();

        std::vector<int> ids;
        
        for(int index = 0; index < mesh.Domains(); index++)
        {
            FEDomain* domain = &mesh.Domain(index);
			FEMaterial* mat = dynamic_cast<FEMaterial*>(domain->GetMaterial());
			if (mat != this->GetParent())
			{
                continue;
            }

            #pragma omp parallel for
            for(int el = 0; el < domain->Elements(); el++)
            {
                FEElement& element = domain->ElementRef(el);

                FEElementODF* odf = new FEElementODF(m_ODF.size());

                // Calculate the centroid of the element
                for(int node = 0; node < element.Nodes(); node++)
                {
                    odf->m_pos += mesh.Node(element.m_node[node]).m_r0;
                }
				odf->m_pos /= (double)element.Nodes();

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

                // interpTime.start();
                // odf->calcODF(pODF);
                // interpTime.stop();

                // Reduce the number of points in the ODF
                // reduceTime.start();
                // reduceODF(odf, B);
                
                

                // reduceTime.stop();

                // Add element odf object to map
                #pragma omp critical
                m_ElODF[element.GetID()] = odf;

                #pragma omp critical
                ids.push_back(element.GetID());

            }
        }

        interpTime.start();
        #pragma omp parallel for
        for(int index = 0; index < ids.size(); index++)
        {
            m_ElODF[ids[index]]->calcODF(pODF);
        }
        interpTime.stop();

        reduceTime.start();
        #pragma omp parallel for
        for(int index = 0; index < ids.size(); index++)
        {
            reduceODF(m_ElODF[ids[index]], B, &remeshTime);
        }
        reduceTime.stop();

    }
    else
    {
        reduceODF(m_ODF[0], B, new double);
    }

    delete[] m_theta; m_theta = nullptr;
    delete[] m_phi; m_phi = nullptr;

    totalTime.stop();

    // std::cout << "Total Time: " << totalTime.peek() << std::endl;
    // std::cout << "Init Time: " << initTime.peek() << std::endl;
    // std::cout << "Interp Time: " << interpTime.peek() << std::endl;
    // std::cout << "Reduce Time: " << reduceTime.peek() << std::endl;
    // std::cout << "Remesh Time: " << remeshTime/36 << std::endl;

	return true;
}

void FEODFFiberDistribution::reduceODF(FEBaseODF* ODF, matrix& B, double* time)
{
    Timer remeshTime;
    vector<double>* sphHar;

    if(dynamic_cast<FEElementODF*>(ODF))
    {
        // Calculate spherical harmonics
        int sphrHarmSize = (m_order+1)*(m_order+2)/2;
        sphHar = new vector<double>(sphrHarmSize);

        B.mult(ODF->m_ODF, *sphHar);
    }
    else
    {
        sphHar = &dynamic_cast<FEFiberODF*>(ODF)->m_shpHar;
    }

    // Calculate the graident
    std::vector<double> gradient;
    altGradient(m_order, ODF->m_ODF, gradient);
    
    // Remesh the sphere
    vector<vec3i> elems;
    remeshTime.start();
    remesh(gradient, m_lengthScale, m_hausd, m_grad, ODF->m_nodePos, elems);
    remeshTime.stop();

    #pragma omp critical
    *time += remeshTime.peek();

    int NN = ODF->m_nodePos.size();
    int NE = elems.size();

    // Store the new coordinates
    double* xCoords = new double[NN] {};
    double* yCoords = new double[NN] {};
    double* zCoords = new double[NN] {};
    for(int index = 0; index < NN; index++)
    {
        vec3d vec = ODF->m_nodePos[index];

        xCoords[index] = vec.x;
        yCoords[index] = vec.y;
        zCoords[index] = vec.z;
    }

    // Convert the new coordinates to spherical coordinates
    double* theta = new double[NN] {};
    double* phi = new double[NN] {};
    getSphereCoords(NN, xCoords, yCoords, zCoords, theta, phi);

    // Compute the new ODF values
    auto T = compSH(m_order, NN, theta, phi); 
    ODF->m_ODF.resize(NN);
    (*T).mult(*sphHar, ODF->m_ODF);

    // Properly scale the ODF values based on the sizes of the surrounding elements
    vector<double> area(NN, 0);
    for(int index = 0; index < NE; index++)
    {
        int n0 = elems[index].x;
        int n1 = elems[index].y;
        int n2 = elems[index].z;

        vec3d n01 = ODF->m_nodePos[n0] - ODF->m_nodePos[n1];
        vec3d n02 = ODF->m_nodePos[n0] - ODF->m_nodePos[n2];

        double temp = (n01^n02).Length()/2;

        area[n0] += temp;
        area[n1] += temp;
        area[n2] += temp;
    }

    // Normalize the ODF
    double sum = 0;
    for(int index = 0; index < NN; index++)
    {
        double val = ODF->m_ODF[index];

        if(ODF->m_nodePos[index].z != 0)
        {
            val *= 2;
        }

        val *= area[index];

        ODF->m_ODF[index] = val;

        sum += val;
    }

    for(int index = 0; index < NN; index++)
    {
        ODF->m_ODF[index] /= sum;
    }

    // Cleanup
    delete[] xCoords;
    delete[] yCoords;
    delete[] zCoords;
    delete[] theta;
    delete[] phi;

    if(dynamic_cast<FEElementODF*>(ODF))
    {
        delete sphHar;
    }
}

//-----------------------------------------------------------------------------
//! Serialization
void FEODFFiberDistribution::Serialize(DumpStream& ar)
{	
	FEElasticMaterial::Serialize(ar);
	if (ar.IsShallow()) return;
}

//-----------------------------------------------------------------------------
//! calculate stress at material point
mat3ds FEODFFiberDistribution::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    FEFiberMaterialPoint& fp = *mp.ExtractData<FEFiberMaterialPoint>();

    FEBaseODF* ODF;
    if(m_interpolate)
    {
        ODF = m_ElODF[mp.m_elem->GetID()];
    }
    else
    {
        ODF = m_ODF[0];
    }

	// calculate stress
	mat3ds s; s.zero();

	// get the local coordinate system
	mat3d Q = GetLocalCS(mp);

    for(int index = 0; index < ODF->m_ODF.size(); index++)
    {
        // evaluate ellipsoidally distributed material coefficients
        double R = ODF->m_ODF[index];
        
        // convert fiber to global coordinates
        vec3d n0 = Q*ODF->m_nodePos[index];
        
        // calculate the stress
        s += m_pFmat->FiberStress(mp, fp.FiberPreStretch(n0))*(R);
    }

	return s;
}

//-----------------------------------------------------------------------------
//! calculate tangent stiffness at material point
tens4ds FEODFFiberDistribution::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    FEFiberMaterialPoint& fp = *mp.ExtractData<FEFiberMaterialPoint>();

    FEBaseODF* ODF;
    if(m_interpolate)
    {
        ODF = m_ElODF[mp.m_elem->GetID()];
    }
    else
    {
        ODF = m_ODF[0];
    }

	// get the local coordinate system
	mat3d Q = GetLocalCS(mp);

    // initialize stress tensor
	tens4ds c;
	c.zero();
    for(int index = 0; index < ODF->m_ODF.size(); index++)
    {
        //evaluate ellipsoidally distributed material coefficients
        double R = ODF->m_ODF[index];
        
        // convert fiber to global coordinates
        vec3d n0 = Q*ODF->m_nodePos[index];

        // calculate the tangent
        c += m_pFmat->FiberTangent(mp, fp.FiberPreStretch(n0))*(R);
	}

	return c;
}

//-----------------------------------------------------------------------------
//! calculate strain energy density at material point
double FEODFFiberDistribution::StrainEnergyDensity(FEMaterialPoint& mp)
{ 
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    FEFiberMaterialPoint& fp = *mp.ExtractData<FEFiberMaterialPoint>();

    FEBaseODF* ODF;
    if(m_interpolate)
    {
        ODF = m_ElODF[mp.m_elem->GetID()];
    }
    else
    {
        ODF = m_ODF[0];
    }

	// get the local coordinate system
	mat3d Q = GetLocalCS(mp);

    double sed = 0.0;
    for(int index = 0; index < ODF->m_ODF.size(); index++)
    {
        // evaluate ellipsoidally distributed material coefficients
        double R = ODF->m_ODF[index];
        
        // convert fiber to global coordinates
        vec3d n0 = Q*ODF->m_nodePos[index];

        // calculate the stress
        sed += m_pFmat->FiberStrainEnergyDensity(mp, fp.FiberPreStretch(n0))*(R);
    }

	return sed;
}