#pragma once
#include "FEBiphasicSolute.h"

//-----------------------------------------------------------------------------
// Monovalent salt (cation + anion)

class FESaltMaterialPoint : public FEMaterialPoint
{
public:
	FESaltMaterialPoint(FEMaterialPoint* ppt) : FEMaterialPoint(ppt) {}
	
	FEMaterialPoint* Copy()
	{
		FESaltMaterialPoint* pt = new FESaltMaterialPoint(*this);
		if (m_pt) pt->m_pt = m_pt->Copy();
		return pt;
	}
	
	void Serialize(DumpFile& ar)
	{
		if (ar.IsSaving())
		{
			ar << m_c[0] << m_gradc[0] << m_j[0] << m_ca[0];
			ar << m_c[1] << m_gradc[1] << m_j[1] << m_ca[1];
		}
		else
		{
			ar >> m_c[0] >> m_gradc[0] >> m_j[0] >> m_ca[0];
			ar >> m_c[1] >> m_gradc[1] >> m_j[1] >> m_ca[1];
		}
		
		if (m_pt) m_pt->Serialize(ar);
	}
	
	void Init(bool bflag)
	{
		if (bflag)
		{
			m_c[0] = m_c[1] = 0;
			m_ca[0] = m_ca[1] = 0;
			m_psi = m_cF = 0;
			m_gradc[0] = m_gradc[1] = vec3d(0,0,0);
			m_j[0] = m_j[1] = m_Ie = vec3d(0,0,0);
		}
		
		if (m_pt) m_pt->Init(bflag);
	}
	
public:
	// salt material data
	double		m_c[2];		//!< effective concentration (0=cation, 1=anion)
	vec3d		m_gradc[2];	//!< spatial gradient of concentration
	vec3d		m_j[2];		//!< solute molar flux
	double		m_ca[2];	//!< actual solute concentration
	double		m_psi;		//!< electric potential
	vec3d		m_Ie;		//!< current density
	double		m_cF;		//!< fixed charge density in current configuration
};

//-----------------------------------------------------------------------------
//! Base class for triphasic materials.

class FETriphasic : public FEMultiMaterial
{
public:
	FETriphasic(FEModel* pfem);

	void Init();
	
	// returns a pointer to a new material point object
	FEMaterialPoint* CreateMaterialPointData() 
	{ 
		return new FESaltMaterialPoint(new FEBiphasicMaterialPoint(m_pSolid->CreateMaterialPointData()));
	}

	// Get the elastic component (overridden from FEMaterial)
	FEElasticMaterial* GetElasticMaterial() { return m_pSolid->GetElasticMaterial(); }

	//! find a material parameter
	FEParam* GetParameter(const ParamString& s);

public:

	//! calculate stress at material point
	mat3ds Stress(FEMaterialPoint& pt);
	
	//! calculate tangent stiffness at material point
	tens4ds Tangent(FEMaterialPoint& pt);
	
	//! calculate fluid (solvent) flux
	vec3d FluidFlux(FEMaterialPoint& pt);
	
	//! calculate solute molar flux
	vec3d SoluteFlux(FEMaterialPoint& pt, const int ion);
	
	//! actual fluid pressure (as opposed to effective pressure)
	double Pressure(FEMaterialPoint& pt);
	
	//! actual concentration (as opposed to effective concentration)
	double Concentration(FEMaterialPoint& pt, const int ion);
	
	//! porosity
	double Porosity(FEMaterialPoint& pt);
	
	//! fixed charge density
	double FixedChargeDensity(FEMaterialPoint& pt);
	
	//! electric potential
	double ElectricPotential(FEMaterialPoint& pt, const bool eform=false);
	
	//! current density
	vec3d CurrentDensity(FEMaterialPoint& pt);
	
	//! fluid density
	double FluidDensity() { return m_rhoTw; }
	
	//! solute density
	double SoluteDensity(const int ion) { return m_pSolute[ion]->Density(); }
	
	//! solute molar mass
	double SoluteMolarMass(const int ion) { return m_pSolute[ion]->MolarMass(); }
	
	//! solute charge number
	int SoluteChargeNumber(const int ion) { return m_pSolute[ion]->ChargeNumber(); }
	
	//! Serialization
	void Serialize(DumpFile& ar);

	//! Add a solute component
	void AddSolute(FESolute* ps);

public:
	//! return number of material properties
	int Properties();

	//! return a material property
	FEMaterial* GetProperty(int i);

	//! find a material property index ( returns <0 for error)
	int FindPropertyIndex(const char* szname);

	//! set a material property (returns false on error)
	bool SetProperty(int i, FEMaterial* pm);
	
public: // material parameters
	double						m_phi0;			//!< solid volume fraction in reference configuration
	double						m_rhoTw;		//!< true fluid density
	double						m_cFr;			//!< fixed charge density in reference configurations
	double						m_Rgas;			//!< universal gas constant
	double						m_Tabs;			//!< absolute temperature
	double						m_Fc;			//!< Faraday's constant
	double						m_penalty;		//!< penalty for enforcing electroneutrality

public: // material properties
	FEElasticMaterial*			m_pSolid;		//!< pointer to elastic solid material
	FEHydraulicPermeability*	m_pPerm;		//!< pointer to permeability material
	FEOsmoticCoefficient*		m_pOsmC;		//!< pointer to osmotic coefficient material
	vector<FESolute*>			m_pSolute;		//!< pointer to solute materials

	DECLARE_PARAMETER_LIST();
};
