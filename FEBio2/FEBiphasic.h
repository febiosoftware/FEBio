#pragma once
#include "FECore/FEElasticMaterial.h"

//-----------------------------------------------------------------------------
// Biphasic material point class.
//
class FEBiphasicMaterialPoint : public FEMaterialPoint
{
public:
	FEBiphasicMaterialPoint(FEMaterialPoint* ppt) : FEMaterialPoint(ppt) {}

	FEMaterialPoint* Copy()
	{
		FEBiphasicMaterialPoint* pt = new FEBiphasicMaterialPoint(*this);
		if (m_pt) pt->m_pt = m_pt->Copy();
		return pt;
	}

	void Serialize(DumpFile& ar)
	{
		if (ar.IsSaving())
		{
			ar << m_p << m_gradp << m_w << m_pa << m_phiw;
		}
		else
		{
			ar >> m_p >> m_gradp >> m_w >> m_pa >> m_phiw;
		}

		if (m_pt) m_pt->Serialize(ar);
	}

	void Init(bool bflag)
	{
		if (bflag)
		{
			m_p = m_pa = 0;
			m_gradp = vec3d(0,0,0);
			m_w = vec3d(0,0,0);
			m_phiw = 1;
		}

		if (m_pt) m_pt->Init(bflag);
	}

public:
	// poro-elastic material data
	// The actual fluid pressure is the same as the effective fluid pressure
	// in a poroelastic material without solute(s).  The actual fluid pressure
	// is included here so that models that include both poroelastic and
	// solute-poroelastic domains produce plotfiles with consistent fluid
	// pressure fields.
	double		m_p;		//!< fluid pressure
	vec3d		m_gradp;	//!< spatial gradient of p
	vec3d		m_w;		//!< fluid flux
	double		m_pa;		//!< actual fluid pressure
	double		m_phiw;		//!< porosity in current configuration
};

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
		return new FEBiphasicMaterialPoint(m_pSolid->CreateMaterialPointData());
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
