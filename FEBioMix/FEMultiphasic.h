#pragma once
#include "FEBiphasic.h"
#include "FESolutesMaterialPoint.h"
#include "FESolute.h"
#include "FEOsmoticCoefficient.h"
#include "FEChemicalReaction.h"
#include "FEMembraneReaction.h"
#include "FESoluteInterface.h"

//-----------------------------------------------------------------------------
//! Base class for multiphasic materials.

class FEMultiphasic : public FEMaterial, public FESoluteInterface
{
public:
	//! constructor
	FEMultiphasic(FEModel* pfem);

	//! initialization
	bool Init() override;

	//! Serialization
	void Serialize(DumpStream& ar) override;

	// returns a pointer to a new material point object
	virtual FEMaterialPoint* CreateMaterialPointData() override = 0;
	
	// return elastic material component
	FEElasticMaterial* GetElasticMaterial() override { return m_pSolid->GetElasticMaterial(); }

    //! Update solid bound molecules
    virtual void UpdateSolidBoundMolecules(FEMaterialPoint& mp) = 0;

public:
	
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

	//! partition coefficient
	double PartitionCoefficient(FEMaterialPoint& pt, const int sol);
	
	//! partition coefficients and their derivatives
	void PartitionCoefficientFunctions(FEMaterialPoint& mp, vector<double>& kappa,
									   vector<double>& dkdJ,
									   vector< vector<double> >& dkdc,
                                       vector< vector<double> >& dkdr,
                                       vector< vector<double> >& dkdJr,
                                       vector< vector< vector<double> > >& dkdrc);
	
	//! solid referential apparent density
	double SolidReferentialApparentDensity(FEMaterialPoint& pt);
	
	//! solid referential volume fraction
	double SolidReferentialVolumeFraction(FEMaterialPoint& pt);

	//! actual concentration (as opposed to effective concentration)
	double Concentration(FEMaterialPoint& pt, const int sol);
	
	//! porosity
	double Porosity(FEMaterialPoint& pt);
	
	//! fixed charge density
    virtual double FixedChargeDensity(FEMaterialPoint& pt);
	
	//! electric potential
	double ElectricPotential(FEMaterialPoint& pt, const bool eform=false);
	
	//! current density
	vec3d CurrentDensity(FEMaterialPoint& pt);

	//! fluid true density
	double FluidDensity() { return m_rhoTw; }
	
	//! solute density
	double SoluteDensity(const int sol) { return m_pSolute[sol]->Density(); }
	
	//! solute molar mass
	double SoluteMolarMass(const int sol) { return m_pSolute[sol]->MolarMass(); }
	
	//! solute charge number
	int SoluteChargeNumber(const int sol) { return m_pSolute[sol]->ChargeNumber(); }
	
	//! SBM density
	double SBMDensity(const int sbm) { return m_pSBM[sbm]->Density(); }
	
	//! SBM molar mass
	double SBMMolarMass(const int sbm) { return m_pSBM[sbm]->MolarMass(); }
	
	//! SBM charge number
	int SBMChargeNumber(const int sbm) { return m_pSBM[sbm]->ChargeNumber(); }
	
	//! SBM actual concentration (molar concentration in current configuration)
	double SBMConcentration(FEMaterialPoint& pt, const int sbm) {
		FEElasticMaterialPoint& ept = *pt.ExtractData<FEElasticMaterialPoint>();
		FEBiphasicMaterialPoint& bpt = *pt.ExtractData<FEBiphasicMaterialPoint>();
		FESolutesMaterialPoint& spt = *pt.ExtractData<FESolutesMaterialPoint>();
		return spt.m_sbmr[sbm]/(ept.m_J-bpt.m_phi0)/SBMMolarMass(sbm);
	}

	//! SBM referential volume fraction
	double SBMReferentialVolumeFraction(FEMaterialPoint& pt, const int sbm) {
		FESolutesMaterialPoint& spt = *pt.ExtractData<FESolutesMaterialPoint>();
		return spt.m_sbmr[sbm]/SBMDensity(sbm);
	}

	//! find local SBM ID from global one
	int FindLocalSBMID(int nid);

	//! Add a solute
	void AddSolute(FESolute* psol);

	//! Add a solid bound molecule
	void AddSolidBoundMolecule(FESolidBoundMolecule* psbm);

	//! Add a chemical reaction
	void AddChemicalReaction(FEChemicalReaction* pcr);
    
    //! Add a membrane reaction
    void AddMembraneReaction(FEMembraneReaction* pcr);

// solute interface
public:
	int Solutes() override { return (int)m_pSolute.size(); }
	FESolute* GetSolute(int i) override { return m_pSolute[i]; }

public:
	FEElasticMaterial*			GetSolid()				{ return m_pSolid; }
	FEHydraulicPermeability*	GetPermeability()		{ return m_pPerm;  }
	FEOsmoticCoefficient*		GetOsmoticCoefficient() { return m_pOsmC;  }
	FESolventSupply*			GetSolventSupply()		{ return m_pSupp;  }
	FESolidBoundMolecule*		GetSBM				(int i) { return m_pSBM[i];    }
	FEChemicalReaction*			GetReaction			(int i) { return m_pReact[i];  }
    FEMembraneReaction*         GetMembraneReaction (int i) { return m_pMReact[i]; }

	int SBMs		     () { return (int) m_pSBM.size();	}
	int Reactions	     () { return (int) m_pReact.size();	}
    int MembraneReactions() { return (int) m_pMReact.size();}

public: // parameters
	double	m_phi0;			//!< solid volume fraction in reference configuration
	double	m_rhoTw;		//!< true fluid density
	double	m_penalty;		//!< penalty for enforcing electroneutrality
	double	m_cFr;			//!< fixed charge density in reference configurations

public:
	double	m_Rgas;			//!< universal gas constant
	double	m_Tabs;			//!< absolute temperature
	double	m_Fc;			//!< Faraday's constant
	int		m_zmin;			//!< minimum charge number in mixture
	int		m_ndeg;			//!< polynomial degree of zeta in electroneutrality

protected:
	// material properties
	FEPropertyT<FEElasticMaterial>			m_pSolid;		//!< pointer to elastic solid material
	FEPropertyT<FEHydraulicPermeability>	m_pPerm;		//!< pointer to permeability material
	FEPropertyT<FEOsmoticCoefficient>		m_pOsmC;		//!< pointer to osmotic coefficient material
	FEPropertyT<FESolventSupply>			m_pSupp;		//!< pointer to solvent supply material
	FEVecPropertyT<FESolute>				m_pSolute;		//!< pointer to solute materials
	FEVecPropertyT<FESolidBoundMolecule>	m_pSBM;			//!< pointer to solid-bound molecule materials
	FEVecPropertyT<FEChemicalReaction>		m_pReact;		//!< pointer to chemical reactions
    FEVecPropertyT<FEMembraneReaction>      m_pMReact;      //!< pointer to membrane reactions

	DECLARE_PARAMETER_LIST();
};
