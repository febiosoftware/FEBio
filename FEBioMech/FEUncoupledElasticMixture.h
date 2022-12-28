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
#include "FEUncoupledMaterial.h"
#include "FEElasticMixture.h"

//-----------------------------------------------------------------------------
//! Uncoupled elastic mixtures

//! This class describes a mixture of uncoupled elastic solids.  The user must declare
//! uncoupled elastic solids that can be combined within this class.  The stress and
//! tangent tensors evaluated in this class represent the sum of the respective
//! tensors of all the solids forming the mixture.
//! \todo This class defines two accessor interfaces. Modify to use the FEMaterial interface only.

class FEUncoupledElasticMixture : public FEUncoupledMaterial
{
public:
	FEUncoupledElasticMixture(FEModel* pfem);

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
	mat3ds DevStress(FEMaterialPoint& pt) override;
	
	//! calculate tangent stiffness at material point
	tens4ds DevTangent(FEMaterialPoint& pt) override;
	
	//! calculate strain energy density at material point
	double DevStrainEnergyDensity(FEMaterialPoint& pt) override;
    
	//! the density is the sum of the constituent densities
	double Density(FEMaterialPoint& mp) override;

public:
    double StrongBondDevSED(FEMaterialPoint& pt) override;
    double WeakBondDevSED(FEMaterialPoint& pt) override;

private:
	// TODO: temporarily reverted back to uncoupled materials. This was needed to make sure that 
	//       FEBio Studio displays the uncoupled materials as options. 
	//       Need to figure out a way to allow elastic materials again.
	std::vector<FEUncoupledMaterial*>	m_pMat;	//!< pointers to elastic materials
//	std::vector<FEElasticMaterial*>	m_pMat;	//!< pointers to elastic materials

	DECLARE_FECORE_CLASS();
};
