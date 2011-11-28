#pragma once
#include "FEBiphasic.h"

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
//! Base class for solute supply.
//! These materials need to define the solute supply and tangent supply functions.
//! The solute supply has units of moles/(referential mixture volume)/time
//!
class FESoluteSupply : public FEMaterial
{
public:
	virtual void Init() {}
	
	//! solute supply
	virtual double Supply(FEMaterialPoint& pt) = 0;
	
	//! solute supply under steady-state conditions
	virtual double SupplySS(FEMaterialPoint& pt) = 0;
	
	//! tangent of solute supply with respect to strain
	virtual double Tangent_Supply_Strain(FEMaterialPoint& mp) = 0;
	
	//! tangent of solute supply with respect to solute concentration
	virtual double Tangent_Supply_Concentration(FEMaterialPoint& mp) = 0;
	
	//! receptor-ligand complex supply
	virtual double ReceptorLigandSupply(FEMaterialPoint& pt) = 0;
	
	//! receptor-ligand concentration under steady-state conditions
	virtual double ReceptorLigandConcentrationSS(FEMaterialPoint& pt) = 0;
	
	//! referential solid supply
	virtual double SolidSupply(FEMaterialPoint& pt) = 0;
	
	//! referential solid concentration under steady-state conditions
	virtual double SolidConcentrationSS(FEMaterialPoint& pt) = 0;
};

//-----------------------------------------------------------------------------

class FESoluteMaterialPoint : public FEMaterialPoint
{
public:
	FESoluteMaterialPoint(FEMaterialPoint* ppt) : FEMaterialPoint(ppt) {}
	
	FEMaterialPoint* Copy()
	{
		FESoluteMaterialPoint* pt = new FESoluteMaterialPoint(*this);
		if (m_pt) pt->m_pt = m_pt->Copy();
		return pt;
	}
	
	void Serialize(DumpFile& ar)
	{
		if (ar.IsSaving())
		{
			ar << m_c << m_gradc << m_j << m_ca << m_crc << m_crcp;
		}
		else
		{
			ar >> m_c >> m_gradc >> m_j >> m_ca >> m_crc >> m_crcp;
		}
		
		if (m_pt) m_pt->Serialize(ar);
	}
	
	void Init(bool bflag)
	{
		if (bflag)
		{
			m_c = m_ca = 0;
			m_gradc = vec3d(0,0,0);
			m_j = vec3d(0,0,0);
			m_crc = m_crcp = 0;
		}
		
		if (m_pt) m_pt->Init(bflag);
	}
	
public:
	// solute material data
	double		m_c;		//!< effective solute concentration
	vec3d		m_gradc;	//!< spatial gradient of c
	vec3d		m_j;		//!< solute molar flux
	double		m_ca;		//!< actual solute concentration
	double		m_crc;		//!< referential concentration of receptor-ligand complex
	double		m_crcp;		//!< m_crc at previous time point
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
		return new FESoluteMaterialPoint
		(new FEBiphasicMaterialPoint(m_pSolid->CreateMaterialPointData()));
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

	//! referential concentration (normalized to mixture volume in reference state)
	double ReferentialConcentration(FEMaterialPoint& pt);

	//! porosity
	double Porosity(FEMaterialPoint& pt);
	
	//! fluid density
	double FluidDensity() { return m_rhoTw; }
	
	//! solute density
	double SoluteDensity() { return m_rhoTu; }
	
	//! solute molecular weight
	double SoluteMolarMass() { return m_Mu; }

	//! Serialization
	void Serialize(DumpFile& ar);
	
public:
	double						m_rhoTw;		//!< true fluid density
	double						m_rhoTu;		//!< true solute density
	double						m_Mu;			//!< solute molecular weight
	double						m_phi0;			//!< solid volume fraction in reference configuration
	FEElasticMaterial*			m_pSolid;		//!< pointer to elastic solid material
	FEHydraulicPermeability*	m_pPerm;		//!< pointer to permeability material
	FESoluteDiffusivity*		m_pDiff;		//!< pointer to diffusivity material
	FESoluteSolubility*			m_pSolub;		//!< pointer to solubility material
	FEOsmoticCoefficient*		m_pOsmC;		//!< pointer to osmotic coefficient material
	FESoluteSupply*				m_pSupp;		//!< pointer to solute supply material
	double						m_Rgas;			//!< universal gas constant
	double						m_Tabs;			//!< absolute temperature
	
	// declare as registered
	DECLARE_REGISTERED(FEBiphasicSolute);
	
	DECLARE_PARAMETER_LIST();
};
