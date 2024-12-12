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



#pragma once
#include <FECore/FECoreClass.h>
#include "FEElasticMaterial.h"
#include "FEFiberMaterial.h"
#include "FEFiberDensityDistribution.h"
#include "FEFiberIntegrationScheme.h"
#include "FEFiberMaterialPoint.h"
#include <unordered_map>

class FEODFFiberDistribution;

class FEBaseODF
{
public:
    FEBaseODF() : m_pos(0,0,0) {}
    virtual ~FEBaseODF() {}

public:
    vec3d	m_pos;	//!< position
    std::vector<double> m_ODF; //!< ODF values
    std::vector<vec3d> m_nodePos; //!< node positions
};

class FEFiberODF : public FECoreClass, public FEBaseODF
{
public:
	FEFiberODF(FEModel* fem) : FECoreClass(fem) {}

public:
    std::vector<double> m_shpHar;	//!< spherical harmonics Values

	DECLARE_FECORE_CLASS();
	FECORE_BASE_CLASS(FEFiberODF);
};

class FEElementODF : public FEBaseODF
{
public:
    FEElementODF(int weightsSize) : m_weights(weightsSize) {}

    void calcODF(std::vector<std::vector<double>>& ODFs);

public:
    std::vector<double> m_weights; //!< weights for interpolation
};

//  This material is a container for a fiber material, a fiber density
//  distribution, and an integration scheme.
class FEODFFiberDistribution : public FEElasticMaterial
{
public:
    FEODFFiberDistribution(FEModel* pfem);
    ~FEODFFiberDistribution();
    
    // returns a pointer to a new material point object
    FEMaterialPointData* CreateMaterialPointData() override;
    
    // Initialization
    bool Init() override;

public:
	//! calculate stress at material point
	mat3ds Stress(FEMaterialPoint& pt) override;
    
	//! calculate tangent stiffness at material point
	tens4ds Tangent(FEMaterialPoint& pt) override;
    
	//! calculate strain energy density at material point
	double StrainEnergyDensity(FEMaterialPoint& pt) override;

	//! Serialization
	void Serialize(DumpStream& ar) override;

protected:
	FEFiberMaterial*			m_pFmat;    // pointer to fiber material
    std::vector<FEFiberODF*>    m_ODF;     // ODF objects
    std::unordered_map<int, FEElementODF*> m_ElODF; // element ODF objects

	DECLARE_FECORE_CLASS();

private:
    void reduceODF(FEBaseODF* ODF, matrix& B, double* time);

private:
    bool m_interpolate; // whether or not we're interpolating at each element
    int m_order; // Spherical Harmonic Order

    std::unique_ptr<matrix> m_T;

    double* m_theta;
    double* m_phi;

    double m_lengthScale;
    double m_hausd;
    double m_grad;
};