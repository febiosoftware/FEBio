#pragma once
#include "FEBiphasic.h"
#include "FESolutesMaterialPoint.h"
#include "FESolute.h"
#include "FEOsmoticCoefficient.h"
#include "FEChemicalReaction.h"

//-----------------------------------------------------------------------------
//! Base class for multiphasic materials.

class FEMultiphasic : public FEMaterial
{
public:
	//! constructor
	FEMultiphasic(FEModel* pfem);

	//! initialization
	void Init();

	// returns a pointer to a new material point object
	FEMaterialPoint* CreateMaterialPointData();
	
	// return elastic material component
	FEElasticMaterial* GetElasticMaterial() { return m_pSolid->GetElasticMaterial(); }

	// find a material parameter
	FEParam* GetParameter(const ParamString& s);
	
public:
	//! return number of material properties
	int Properties();

	//! return a material property
	FECoreBase* GetProperty(int i);

	//! find a material property index ( returns <0 for error)
	int FindPropertyIndex(const char* szname);

	//! set a material property (returns false on error)
	bool SetProperty(int i, FECoreBase* pm);

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
									   vector<double>& dkdJJ,
									   vector< vector<double> >& dkdJc,
									   vector< vector< vector<double> > >& dkdcc,
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
	double FixedChargeDensity(FEMaterialPoint& pt);
	
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
	
	//! Serialization
	void Serialize(DumpFile& ar);

	// initialize chemical reaction
	void InitializeReaction(FEChemicalReaction* m_pReact);

	//! Add a solute
	void AddSolute(FESolute* psol);

	//! Add a solid bound molecule
	void AddSolidBoundMolecule(FESolidBoundMolecule* psbm);

	//! Add a chemical reaction
	void AddChemicalReaction(FEChemicalReaction* pcr);

public:
	FEElasticMaterial*			GetSolid()				{ return m_pSolid; }
	FEHydraulicPermeability*	GetPermeability()		{ return m_pPerm;  }
	FEOsmoticCoefficient*		GetOsmoticCoefficient() { return m_pOsmC;  }
	FESolventSupply*			GetSolventSupply()		{ return m_pSupp;  }
	FESolute*					GetSolute			(int i) { return m_pSolute[i]; }
	FESolidBoundMolecule*		GetSBM				(int i) { return m_pSBM[i];    }
	FEChemicalReaction*			GetReaction			(int i) { return m_pReact[i];  }

	int Solutes		() { return (int) m_pSolute.size(); }
	int SBMs		() { return (int) m_pSBM.size();	}
	int Reactions	() { return (int) m_pReact.size();	}

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

private:
	FEElasticMaterial*				m_pSolid;		//!< pointer to elastic solid material
	FEHydraulicPermeability*		m_pPerm;		//!< pointer to permeability material
	FEOsmoticCoefficient*			m_pOsmC;		//!< pointer to osmotic coefficient material
	FESolventSupply*				m_pSupp;		//!< pointer to solvent supply material
	vector<FESolute*>				m_pSolute;		//!< pointer to solute materials
	vector<FESolidBoundMolecule*>	m_pSBM;			//!< pointer to solid-bound molecule materials
	vector<FEChemicalReaction*>		m_pReact;		//!< pointer to chemical reactions

	DECLARE_PARAMETER_LIST();
};
