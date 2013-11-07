#pragma once
#include "FEBiphasic.h"
#include "FESolute.h"

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
	virtual double Tangent_OsmoticCoefficient_Concentration(FEMaterialPoint& mp, const int isol) = 0;
};

//-----------------------------------------------------------------------------
//! material point data for solute in biphasic-solute material
class FESoluteMaterialPoint : public FEMaterialPoint
{
public:
	//! constructor
	FESoluteMaterialPoint(FEMaterialPoint* ppt);
	
	//! shallow copy 
	FEMaterialPoint* Copy();
	
	//! serialize data
	void Serialize(DumpFile& ar);

	//! Data initialization
	void Init(bool bflag);
	
public:
	// solute material data
	double		m_c;		//!< effective solute concentration
	vec3d		m_gradc;	//!< spatial gradient of c
	vec3d		m_j;		//!< solute molar flux
	double		m_ca;		//!< actual solute concentration
	double		m_crc;		//!< referential concentration of receptor-ligand complex
	double		m_crcp;		//!< m_crc at previous time point
	double		m_crchat;	//!< referential receptor-ligand complex supply
	double		m_crchatp;	//!< m_crchat at previous time point
};

//-----------------------------------------------------------------------------
//! Base class for solute diffusion in biphasic materials.

class FEBiphasicSolute : public FEMultiMaterial
{
public:
	FEBiphasicSolute();
	
	// returns a pointer to a new material point object
	FEMaterialPoint* CreateMaterialPointData() 
	{ 
		return new FESoluteMaterialPoint
		(new FEBiphasicMaterialPoint(m_pSolid->CreateMaterialPointData()));
	}

	// Get the elastic component (overridden from FEMaterial)
	FEElasticMaterial* GetElasticMaterial() { return m_pSolid->GetElasticMaterial(); }

	// find a material parameter
	FEParam* GetParameter(const ParamString& s);

public:
	//! return number of material properties
	int Properties();

	//! return a material property
	FEMaterial* GetProperty(int i);

	//! find a material property index ( returns <0 for error)
	int FindPropertyIndex(const char* szname);

	//! set a material property (returns false on error)
	bool SetProperty(int i, FEMaterial* pm);
	
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
	
	//! Serialization
	void Serialize(DumpFile& ar);
	
public: // material parameters
	double						m_rhoTw;		//!< true fluid density
	double						m_rhoTu;		//!< true solute density
	double						m_Mu;			//!< solute molecular weight
	double						m_phi0;			//!< solid volume fraction in reference configuration
	double						m_Rgas;			//!< universal gas constant
	double						m_Tabs;			//!< absolute temperature

public: // material properties
	FEElasticMaterial*			m_pSolid;		//!< pointer to elastic solid material
	FEHydraulicPermeability*	m_pPerm;		//!< pointer to permeability material
	FEOsmoticCoefficient*		m_pOsmC;		//!< pointer to osmotic coefficient material
	FESolute*					m_pSolute;		//!< pointer to solute material

	DECLARE_PARAMETER_LIST();
};
