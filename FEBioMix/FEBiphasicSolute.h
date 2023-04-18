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
#include "FEBiphasic.h"
#include "FESolutesMaterialPoint.h"
#include "FESolute.h"
#include "FEOsmoticCoefficient.h"
#include "FESoluteInterface.h"
#include <FECore/FEModelParam.h>

//-----------------------------------------------------------------------------
//! Base class for solute diffusion in biphasic materials.

class FEBIOMIX_API FEBiphasicSolute : public FEMaterial, public FEBiphasicInterface, public FESoluteInterface_T<FESolutesMaterialPoint>
{
public:
	FEBiphasicSolute(FEModel* pfem);
	
	// returns a pointer to a new material point object
	FEMaterialPointData* CreateMaterialPointData() override;

	// Get the elastic component (overridden from FEMaterial)
	FEElasticMaterial* GetElasticMaterial() { return m_pSolid; }

	//! Get the solid
	FEElasticMaterial* GetSolid() { return m_pSolid; }

	//! Get the permeability
	FEHydraulicPermeability* GetPermeability() { return m_pPerm; }

// solute interface
public:
	// number of solutes 
	int Solutes() override { return 1; }

	//! Get the solute
	FESolute* GetSolute(int i=0) override { return (i==0 ? (FESolute*)m_pSolute : 0); }
	double GetReferentialFixedChargeDensity(const FEMaterialPoint& mp) override;

	double GetFixedChargeDensity(const FEMaterialPoint& mp) override {
		const FESolutesMaterialPoint* spt = (mp.ExtractData<FESolutesMaterialPoint>());
		return spt->m_cF;
	}

	//! Get the osmotic coefficient
	FEOsmoticCoefficient* GetOsmoticCoefficient() override { return m_pOsmC; }

public: // overridden from FEBiphasicInterface
    //! Evaluate and return solid referential volume fraction
    double SolidReferentialVolumeFraction(FEMaterialPoint& mp) override {
        double phisr = m_phi0(mp);
        FEBiphasicMaterialPoint* pt = (mp.ExtractData<FEBiphasicMaterialPoint>());
        pt->m_phi0 = pt->m_phi0t = phisr;
        return phisr;
    };
    // Return solid referential volume fraction
	double GetReferentialSolidVolumeFraction(const FEMaterialPoint& mp) override {
		const FEBiphasicMaterialPoint* pt = (mp.ExtractData<FEBiphasicMaterialPoint>());
		return pt->m_phi0t;
	}
    // Return actual fluid pressure
	double GetActualFluidPressure(const FEMaterialPoint& mp) override {
		const FEBiphasicMaterialPoint* pt = (mp.ExtractData<FEBiphasicMaterialPoint>());
		return pt->m_pa;
	}

public:
	bool Init() override;
    
    //! specialized material points
    void UpdateSpecializedMaterialPoints(FEMaterialPoint& mp, const FETimeInfo& tp) override;

	//! serialization
	void Serialize(DumpStream& ar) override;
	
	//! calculate stress at material point
	mat3ds Stress(FEMaterialPoint& pt);
	
	//! calculate tangent stiffness at material point
	tens4ds Tangent(FEMaterialPoint& pt);
	
	//! calculate fluid (solvent) flux
	vec3d FluidFlux(FEMaterialPoint& pt);
	
	//! calculate solute molar flux
	vec3d SoluteFlux(FEMaterialPoint& pt);
	
	//! actual fluid pressure (as opposed to effective pressure)
	double Pressure(FEMaterialPoint& pt);
	
	//! actual concentration (as opposed to effective concentration)
	double Concentration(FEMaterialPoint& pt);

	//! referential concentration (normalized to mixture volume in reference state)
	double ReferentialConcentration(FEMaterialPoint& pt);

	//! porosity
	double Porosity(FEMaterialPoint& pt);
	
    //! partition coefficient derivatives
    void PartitionCoefficientFunctions(FEMaterialPoint& mp, double& kappa,
                                       double& dkdJ, double& dkdc);
	//! fluid density
	double FluidDensity() { return m_rhoTw; }
	
public: // material parameters
	double						m_rhoTw;		//!< true fluid density
	FEParamDouble				m_phi0;			//!< solid volume fraction in reference configuration

public:
	double						m_Mu;			//!< solute molecular weight
	double						m_rhoTu;		//!< true solute density
	double						m_Rgas;			//!< universal gas constant
	double						m_Tabs;			//!< absolute temperature

private: // material properties
	FEElasticMaterial*			m_pSolid;		//!< pointer to elastic solid material
	FEHydraulicPermeability*	m_pPerm;		//!< pointer to permeability material
	FEOsmoticCoefficient*		m_pOsmC;		//!< pointer to osmotic coefficient material
	FESolute*					m_pSolute;		//!< pointer to solute material

	DECLARE_FECORE_CLASS();
};
