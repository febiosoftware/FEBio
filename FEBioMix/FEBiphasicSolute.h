#pragma once
#include "FEBiphasic.h"
#include "FESolutesMaterialPoint.h"
#include "FESolute.h"
#include "FEOsmoticCoefficient.h"

//-----------------------------------------------------------------------------
//! Base class for solute diffusion in biphasic materials.

class FEBiphasicSolute : public FEMaterial
{
public:
	FEBiphasicSolute(FEModel* pfem);
	
	// returns a pointer to a new material point object
	FEMaterialPoint* CreateMaterialPointData();

	// Get the elastic component (overridden from FEMaterial)
	FEElasticMaterial* GetElasticMaterial() { return m_pSolid->GetElasticMaterial(); }

	//! Get the solid
	FEElasticMaterial* GetSolid() { return m_pSolid; }

	//! Get the permeability
	FEHydraulicPermeability* GetPermeability() { return m_pPerm; }

	//! Get the osmotic coefficient
	FEOsmoticCoefficient* GetOsmoticCoefficient() { return m_pOsmC; }

	//! Get the solute
	FESolute* GetSolute() { return m_pSolute; }

public:
	//! return number of material properties
	int MaterialProperties();

	//! return a material property
	FEProperty* GetMaterialProperty(int nid);

public:
	void Init();
	
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
	
	//! fluid density
	double FluidDensity() { return m_rhoTw; }
	
public: // material parameters
	double						m_rhoTw;		//!< true fluid density
	double						m_rhoTu;		//!< true solute density
	double						m_Mu;			//!< solute molecular weight
	double						m_phi0;			//!< solid volume fraction in reference configuration
	double						m_Rgas;			//!< universal gas constant
	double						m_Tabs;			//!< absolute temperature

private: // material properties
	FEPropertyT<FEElasticMaterial>			m_pSolid;		//!< pointer to elastic solid material
	FEPropertyT<FEHydraulicPermeability>	m_pPerm;		//!< pointer to permeability material
	FEPropertyT<FEOsmoticCoefficient>		m_pOsmC;		//!< pointer to osmotic coefficient material
	FEPropertyT<FESolute>					m_pSolute;		//!< pointer to solute material

	DECLARE_PARAMETER_LIST();
};
