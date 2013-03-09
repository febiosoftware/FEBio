#pragma once
#include "FEBiphasicSolute.h"
#include "FESolute.h"

//-----------------------------------------------------------------------------
// Multiple solutes

class FESolutesMaterialPoint : public FEMaterialPoint
{
public:
	FESolutesMaterialPoint(FEMaterialPoint* ppt) : FEMaterialPoint(ppt) {}
	
	FEMaterialPoint* Copy()
	{
		FESolutesMaterialPoint* pt = new FESolutesMaterialPoint(*this);
		if (m_pt) pt->m_pt = m_pt->Copy();
		return pt;
	}
	
	void Serialize(DumpFile& ar)
	{
		if (ar.IsSaving())
		{
			for (int i=0; i<m_nsol; ++i)
				ar << m_c[i] << m_gradc[i] << m_j[i] << m_ca[i];
		}
		else
		{
			for (int i=0; i<m_nsol; ++i)
				ar >> m_c[i] >> m_gradc[i] >> m_j[i] >> m_ca[i];
		}
		
		if (m_pt) m_pt->Serialize(ar);
	}
	
	void Init(bool bflag)
	{
		if (bflag)
		{
			m_nsol = 0;
			m_psi = m_cF = 0;
			m_Ie = vec3d(0,0,0);
		}
		
		if (m_pt) m_pt->Init(bflag);
	}
	
public:
	// solutes material data
	int				m_nsol;		//!< number of solutes
	vector<double>	m_c;		//!< effective concentration
	vector<vec3d>	m_gradc;	//!< spatial gradient of concentration
	vector<vec3d>	m_j;		//!< solute molar flux
	vector<double>	m_ca;		//!< actual solute concentration
	double			m_psi;		//!< electric potential
	vec3d			m_Ie;		//!< current density
	double			m_cF;		//!< fixed charge density in current configuration
};

//-----------------------------------------------------------------------------
//! Base class for multiphasic materials.

class FEMultiphasic : public FEMultiMaterial
{
public:
	FEMultiphasic();
	
	// returns a pointer to a new material point object
	FEMaterialPoint* CreateMaterialPointData() 
	{ 
		return new FESolutesMaterialPoint(new FEBiphasicMaterialPoint(m_pSolid->CreateMaterialPointData()));
	}
	
	// return elastic material component
	FEElasticMaterial* GetElasticMaterial() { return m_pSolid->GetElasticMaterial(); }

	// find a material parameter
	FEParam* GetParameter(const ParamString& s);
	
public:
	void Init();
	
	//! calculate stress at material point
	mat3ds Stress(FEMaterialPoint& pt);
	
	//! calculate tangent stiffness at material point
	tens4ds Tangent(FEMaterialPoint& pt);
	
	//! calculate fluid (solvent) flux
	vec3d FluidFlux(FEMaterialPoint& pt);
	
	//! calculate solute molar flux
	vec3d SoluteFlux(FEMaterialPoint& pt, const int sol);
	
	//! actual fluid pressure (as opposed to effective pressure)
	double Pressure(FEMaterialPoint& pt);
	
	//! actual concentration (as opposed to effective concentration)
	double Concentration(FEMaterialPoint& pt, const int sol);
	
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
	double SoluteDensity(const int sol) { return m_pSolute[sol]->Density(); }
	
	//! solute molar mass
	double SoluteMolarMass(const int sol) { return m_pSolute[sol]->MolarMass(); }
	
	//! solute charge number
	int SoluteChargeNumber(const int sol) { return m_pSolute[sol]->ChargeNumber(); }
	
	//! Serialization
	void Serialize(DumpFile& ar);
		
public:
	double						m_phi0;			//!< solid volume fraction in reference configuration
	double						m_rhoTw;		//!< true fluid density
	double						m_cFr;			//!< fixed charge density in reference configurations
	double						m_Rgas;			//!< universal gas constant
	double						m_Tabs;			//!< absolute temperature
	double						m_Fc;			//!< Faraday's constant
	double						m_penalty;		//!< penalty for enforcing electroneutrality
	int							m_zmin;			//!< minimum charge number in mixture
	int							m_ndeg;			//!< polynomial degree of zeta in electroneutrality

public:
	FEElasticMaterial*			m_pSolid;		//!< pointer to elastic solid material
	FEHydraulicPermeability*	m_pPerm;		//!< pointer to permeability material
	FEOsmoticCoefficient*		m_pOsmC;		//!< pointer to osmotic coefficient material
	FESolventSupply*			m_pSupp;		//!< pointer to solvent supply material
	vector<FESolute*>			m_pSolute;		//!< pointer to solute materials

	// declare as registered
	DECLARE_REGISTERED(FEMultiphasic);
	
	DECLARE_PARAMETER_LIST();
};
