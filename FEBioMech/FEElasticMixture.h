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
#include "FEElasticMaterial.h"

//-----------------------------------------------------------------------------
//! Material point data for mixtures
//!
class FEElasticMixtureMaterialPoint : public FEMaterialPointArray
{
public:
	//! constructor
	FEElasticMixtureMaterialPoint();

	//! Copy material point data
	FEMaterialPointData* Copy() override;

	//! material point initialization
	void Init() override;

	//! data serialization
	void Serialize(DumpStream& ar) override;

public:
	vector<double>				m_w;	//!< material weights
};

//-----------------------------------------------------------------------------
//! Elastic mixtures

//! This class describes a mixture of elastic solids.  The user must declare
//! elastic solids that can be combined within this class.  The stress and
//! tangent tensors evaluated in this class represent the sum of the respective
//! tensors of all the solids forming the mixture.

//! \todo This class defines two accessor interfaces. Modify to use the FEMaterial interface only.
class FEElasticMixture : public FEElasticMaterial
{
public:
	FEElasticMixture(FEModel* pfem);

	// returns a pointer to a new material point object
	FEMaterialPointData* CreateMaterialPointData() override;

	// return number of materials
	int Materials() { return (int)m_pMat.size(); }

	// return a material component
	FEElasticMaterial* GetMaterial(int i) { return m_pMat[i]; }

	// Add a material component
	void AddMaterial(FEElasticMaterial* pm);

    //! specialized material points
    void UpdateSpecializedMaterialPoints(FEMaterialPoint& mp, const FETimeInfo& tp) override;
    
public:
    //! data initialization
    bool Init() override;
    
	//! calculate stress at material point
	mat3ds Stress(FEMaterialPoint& pt) override;
		
	//! calculate tangent stiffness at material point
	tens4ds Tangent(FEMaterialPoint& pt) override;
		
	//! calculate strain energy density at material point
	double StrainEnergyDensity(FEMaterialPoint& pt) override;

public:
    double StrongBondSED(FEMaterialPoint& pt) override;
    double WeakBondSED(FEMaterialPoint& pt) override;
    
private:
	std::vector<FEElasticMaterial*>	m_pMat;	//!< pointers to elastic materials

	DECLARE_FECORE_CLASS();
};
