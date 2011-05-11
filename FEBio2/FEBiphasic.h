#pragma once
#include "FECore/FEMaterial.h"

//-----------------------------------------------------------------------------
//! Base class for hydraulic permeability of porous materials.
//! These materials need to define the permeability and tangent permeability functions.
//!
class FEHydraulicPermeability : public FEMaterial
	{
	public:
		FEHydraulicPermeability() {m_phi0 = -1; }
		virtual ~FEHydraulicPermeability(){}
		
		//! hydraulic permeability
		virtual mat3ds Permeability(FEMaterialPoint& pt) = 0;
		
		//! tangent of hydraulic permeability with respect to strain
		virtual tens4ds Tangent_Permeability_Strain(FEMaterialPoint& mp) = 0;
		
		//! tangent of hydraulic permeability with respect to concentration
		mat3ds Tangent_Permeability_Concentration(FEMaterialPoint& mp);
		
		void Init();
		
	public:
		double	m_phi0;			//!< solid volume fraction in reference state
		
		// declare parameter list
		DECLARE_PARAMETER_LIST();
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
	
	//! calculate fluid flux
	vec3d Flux(FEMaterialPoint& pt);
	
	//! calculate actual fluid pressure
	double Pressure(FEMaterialPoint& pt);
	
	//! permeability
	void Permeability(double k[3][3], FEMaterialPoint& pt);
	
	//! tangent of permeability
	tens4ds Tangent_Permeability_Strain(FEMaterialPoint& pt);
	
	//! porosity
	double Porosity(FEMaterialPoint& pt);
	
	//! fluid density
	double FluidDensity() { return m_rhoTw; } 

	//! Serialization
	void Serialize(DumpFile& ar);
	
public:
	double						m_rhoTw;	//!< true fluid density
	FEElasticMaterial*			m_pSolid;	//!< pointer to elastic solid material
	FEHydraulicPermeability*	m_pPerm;	//!< pointer to permeability material
	
	// declare as registered
	DECLARE_REGISTERED(FEBiphasic);
	
	DECLARE_PARAMETER_LIST();
	
};

//-----------------------------------------------------------------------------
//! Base class for solute diffusivity.
//! These materials need to define the diffusivity and tangent diffusivity functions.
//!
class FESoluteDiffusivity : public FEMaterial
	{
	public:
		//! solute diffusivity
		virtual mat3ds Diffusivity(FEMaterialPoint& pt) = 0;
		
		//! tangent of diffusivity with respect to strain
		virtual tens4ds Tangent_Diffusivity_Strain(FEMaterialPoint& mp) = 0;
		
		//! tangent of diffusivity with respect to solute concentration
		virtual mat3ds Tangent_Diffusivity_Concentration(FEMaterialPoint& mp) = 0;
		
		//! solute diffusivity in free solution
		virtual double Free_Diffusivity(FEMaterialPoint& pt) = 0;
		
	};

//-----------------------------------------------------------------------------
//! Base class for solute solubility.
//! These materials need to define the solubility and tangent solubility functions.
//!
class FESoluteSolubility : public FEMaterial
	{
	public:
		//! solute solubility
		virtual double Solubility(FEMaterialPoint& pt) = 0;
		
		//! tangent of solubility with respect to strain
		virtual double Tangent_Solubility_Strain(FEMaterialPoint& mp) = 0;
		
		//! tangent of solubility with respect to concentration
		virtual double Tangent_Solubility_Concentration(FEMaterialPoint& mp) = 0;
		
		//! cross derivative of solubility with respect to strain and concentration
		virtual double Tangent_Solubility_Strain_Concentration(FEMaterialPoint& mp) = 0;
		
		//! second derivative of solubility with respect to strain
		virtual double Tangent_Solubility_Strain_Strain(FEMaterialPoint& mp) = 0;
		
	};

//-----------------------------------------------------------------------------
//! Base class for osmotic coefficient.
//! These materials need to define the osmotic coefficient and tangent functions.
//!
class FEOsmoticCoefficient : public FEMaterial
	{
	public:
		//! osmotic coefficient
		virtual double OsmoticCoefficient(FEMaterialPoint& pt) = 0;
		
		//! tangent of osmotic coefficient with respect to strain
		virtual double Tangent_OsmoticCoefficient_Strain(FEMaterialPoint& mp) = 0;
		
		//! tangent of osmotic coefficient with respect to concentration
		virtual double Tangent_OsmoticCoefficient_Concentration(FEMaterialPoint& mp) = 0;
		
	};

//-----------------------------------------------------------------------------
//! Base class for solute diffusion in biphasic materials.

class FEBiphasicSolute : public FEMaterial
{
public:
	FEBiphasicSolute();
	
	// returns a pointer to a new material point object
	FEMaterialPoint* CreateMaterialPointData() 
	{ 
		return new FESolutePoroElasticMaterialPoint(m_pSolid->CreateMaterialPointData());
	}
	
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
	
	//! porosity
	double Porosity(FEMaterialPoint& pt);
	
	//! fluid density
	double FluidDensity() { return m_rhoTw; }
	
	//! solute density
	double SoluteDensity() { return m_rhoTu; }
	
	//! solute molecular weight
	double SoluteMolecularWeight() { return m_Mu; }

	//! Serialization
	void Serialize(DumpFile& ar);
	
public:
	double						m_rhoTw;		//!< true fluid density
	double						m_rhoTu;		//!< true solute density
	double						m_Mu;			//!< solute molecular weight
	FEElasticMaterial*			m_pSolid;		//!< pointer to elastic solid material
	FEHydraulicPermeability*	m_pPerm;		//!< pointer to permeability material
	FESoluteDiffusivity*		m_pDiff;		//!< pointer to diffusivity material
	FESoluteSolubility*			m_pSolub;		//!< pointer to solubility material
	FEOsmoticCoefficient*		m_pOsmC;		//!< pointer to osmotic coefficient material
	double						m_Rgas;			//!< universal gas constant
	double						m_Tabs;			//!< absolute temperature
	
	// declare as registered
	DECLARE_REGISTERED(FEBiphasicSolute);
	
	DECLARE_PARAMETER_LIST();
	
};
