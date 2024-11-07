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
#include "FEBioMech/FEElasticMaterial.h"
#include "FEBioMech/FERemodelingElasticMaterial.h"
#include "FEMultiphasic.h"
#include "febiomix_api.h"

//-----------------------------------------------------------------------------
//! This is a container material for an elastic material associated with a SBM. Its strain energy density, stress
//! and elasticity are scaled by the mass fraction of the SBM, relative to the SBM's true density.

class FEBIOMIX_API FERemodelingSolid : public FEElasticMaterial, public FERemodelingInterface
{
public:
    FERemodelingSolid(FEModel* pfem) : FEElasticMaterial(pfem) { m_sbm = -1; m_lsbm = -1; m_pMat = nullptr; m_pMP = nullptr; }
	
protected:
	int		m_lsbm;	//!< local id of solid-bound molecule

public:
	//! data initialization and checking
	bool Init() override;

	//! serialization
	void Serialize(DumpStream& ar) override;

	//! calculate stress at material point
	mat3ds Stress(FEMaterialPoint& pt) override;
	
	//! calculate tangent stiffness at material point
	tens4ds Tangent(FEMaterialPoint& pt) override;

    //! evaluate referential mass density
    double Density(FEMaterialPoint& pt) override;
    
	//! Create material point data
	FEMaterialPointData* CreateMaterialPointData() override;
	
    //! calculate strain energy density at material point
    double StrainEnergyDensity(FEMaterialPoint& pt) override { return StrainEnergy(pt); }
    
    //! update specialize material point data
    void UpdateSpecializedMaterialPoints(FEMaterialPoint& mp, const FETimeInfo& tp) override;
    
public: // --- remodeling interface ---

	//! calculate strain energy density at material point
	double StrainEnergy(FEMaterialPoint& pt) override;

	//! calculate tangent of strain energy density with mass density
	double Tangent_SE_Density(FEMaterialPoint& pt) override;
	
	//! calculate tangent of stress with mass density
	mat3ds Tangent_Stress_Density(FEMaterialPoint& pt) override;

public:
    FEElasticMaterial*  m_pMat; //!< elastic material which obeys simple remodeling rule

private:
    FEMultiphasic*      m_pMP;  //!< multiphasic domain containing this material

    // declare the parameter list
	DECLARE_FECORE_CLASS();
};
