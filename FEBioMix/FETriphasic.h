#pragma once
#include "FEMultiphasic.h"

//-----------------------------------------------------------------------------
//! Base class for triphasic materials.

class FETriphasic : public FEMaterial
{
public:
	FETriphasic(FEModel* pfem);

	// initialization
	bool Init();

	// serialization
	void Serialize(DumpStream& ar);
	
	// returns a pointer to a new material point object
	FEMaterialPoint* CreateMaterialPointData() 
	{ 
		return new FESolutesMaterialPoint(new FEBiphasicMaterialPoint(m_pSolid->CreateMaterialPointData()));
	}

	// Get the elastic component (overridden from FEMaterial)
	FEElasticMaterial* GetElasticMaterial() { return m_pSolid->GetElasticMaterial(); }

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
	
    //! partition coefficient
    double PartitionCoefficient(FEMaterialPoint& pt, const int sol);
    
    //! partition coefficient derivatives
    void PartitionCoefficientFunctions(FEMaterialPoint& mp, vector<double>& kappa,
                                       vector<double>& dkdJ,
                                       vector< vector<double> >& dkdc);
    //! fluid density
	double FluidDensity() { return m_rhoTw; }
	
	//! solute density
	double SoluteDensity(const int ion) { return m_pSolute[ion]->Density(); }
	
	//! solute molar mass
	double SoluteMolarMass(const int ion) { return m_pSolute[ion]->MolarMass(); }
	
	//! solute charge number
	int SoluteChargeNumber(const int ion) { return m_pSolute[ion]->ChargeNumber(); }
	
	//! Add a solute component
	void AddSolute(FESolute* ps);

public:
    FEElasticMaterial*			GetSolid()				{ return m_pSolid; }
    FEHydraulicPermeability*	GetPermeability()		{ return m_pPerm;  }
    FEOsmoticCoefficient*		GetOsmoticCoefficient() { return m_pOsmC;  }
    FESolute*					GetSolute(int i)        { return m_pSolute[i]; }
    
    int Solutes		() { return (int) m_pSolute.size(); }

public: // material parameters
	double						m_phi0;			//!< solid volume fraction in reference configuration
	double						m_rhoTw;		//!< true fluid density
	double						m_cFr;			//!< fixed charge density in reference configurations
	double						m_penalty;		//!< penalty for enforcing electroneutrality

public:
	double						m_Rgas;			//!< universal gas constant
	double						m_Tabs;			//!< absolute temperature
	double						m_Fc;			//!< Faraday's constant

public: // material properties
	FEPropertyT<FEElasticMaterial>			m_pSolid;		//!< pointer to elastic solid material
	FEPropertyT<FEHydraulicPermeability>	m_pPerm;		//!< pointer to permeability material
	FEPropertyT<FEOsmoticCoefficient>		m_pOsmC;		//!< pointer to osmotic coefficient material
	FEVecPropertyT<FESolute>				m_pSolute;		//!< pointer to solute materials

	DECLARE_PARAMETER_LIST();
};
