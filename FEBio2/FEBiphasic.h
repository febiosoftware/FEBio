#pragma once
#include "FECore/FEMaterial.h"

//-----------------------------------------------------------------------------
//! Base class for hydraulic permeability of porous materials.
//! These materials need to define the permeability and tangent permeability functions.
//!
class FEHydraulicPermeability : public FEMaterial
	{
	public:
		FEHydraulicPermeability() {}
		virtual ~FEHydraulicPermeability(){}
		
		//! hydraulic permeability
		virtual mat3ds Permeability(FEMaterialPoint& pt) = 0;
		
		//! tangent of hydraulic permeability with respect to strain
		virtual tens4ds Tangent_Permeability_Strain(FEMaterialPoint& mp) = 0;
		
		//! tangent of hydraulic permeability with respect to concentration
		mat3ds Tangent_Permeability_Concentration(FEMaterialPoint& mp);
		
		void Init();
		
	};

//-----------------------------------------------------------------------------
//! Base class for biphasic materials.

class FEBiphasic : public FEMaterial
{
public:
	FEBiphasic();
	
	// returns a pointer to a new material point object
	FEMaterialPoint* CreateMaterialPointData() 
	{ 
		return new FEPoroElasticMaterialPoint(m_pSolid->CreateMaterialPointData());
	}
	
public:
	void Init();
	
	//! calculate stress at material point
	mat3ds Stress(FEMaterialPoint& pt);
	
	//! calculate tangent stiffness at material point
	tens4ds Tangent(FEMaterialPoint& pt);

	//! return the permeability tensor as a matrix
	void Permeability(double k[3][3], FEMaterialPoint& pt);
	
	//! calculate fluid flux
	vec3d Flux(FEMaterialPoint& pt);
	
	//! calculate actual fluid pressure
	double Pressure(FEMaterialPoint& pt);

	//! porosity
	double Porosity(FEMaterialPoint& pt);
	
	//! fluid density
	double FluidDensity() { return m_rhoTw; } 

	//! Serialization
	void Serialize(DumpFile& ar);
	
public:
	double						m_rhoTw;	//!< true fluid density
	double						m_phi0;		//!< solid volume fraction in reference configuration
	FEElasticMaterial*			m_pSolid;	//!< pointer to elastic solid material
	FEHydraulicPermeability*	m_pPerm;	//!< pointer to permeability material
	
	// declare as registered
	DECLARE_REGISTERED(FEBiphasic);
	
	DECLARE_PARAMETER_LIST();
	
};
