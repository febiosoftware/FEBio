#pragma once
#include "FEBiphasicSolute.h"
#include <map>

class FEMultiphasic;
class FEChemicalReaction;

//-----------------------------------------------------------------------------
//! Class for storing material point data for solute materials

class FESolutesMaterialPoint : public FEMaterialPoint
{
public:
	//! Constructor
	FESolutesMaterialPoint(FEMaterialPoint* ppt) : FEMaterialPoint(ppt) {}
	
	//! Create a shallow copy
	FEMaterialPoint* Copy();
	
	//! serialize data
	void Serialize(DumpFile& ar);
	
	//! Initialize material point data
	void Init(bool bflag);
	
public:
	// solutes material data
	int				m_nsol;		//!< number of solutes
	vector<double>	m_c;		//!< effective solute concentration
	vector<vec3d>	m_gradc;	//!< spatial gradient of solute concentration
	vector<vec3d>	m_j;		//!< solute molar flux
	vector<double>	m_ca;		//!< actual solute concentration
	double			m_psi;		//!< electric potential
	vec3d			m_Ie;		//!< current density
	double			m_cF;		//!< fixed charge density in current configuration
	int				m_nsbm;		//!< number of solid-bound molecules
	vector<double>	m_sbmr;		//!< referential mass concentration of solid-bound molecules
	vector<double>	m_sbmrp;	//!< m_sbmr at previoust time step
	vector<double>	m_sbmrhat;	//!< referential mass supply of solid-bound molecules
	vector<double>	m_sbmrmin;	//!< minimum value of m_sbmr
	vector<double>	m_sbmrmax;	//!< maximum value of m_sbmr
	vector<double>	m_k;		//!< solute partition coefficient
	vector<double>	m_dkdJ;		//!< 1st deriv of m_k with strain (J)
	vector<double>	m_dkdJJ;	//!< 2nd deriv of m_k with strain (J)
	vector< vector<double> >	m_dkdc;			//!< 1st deriv of m_k with effective concentration
	vector< vector<double> >	m_dkdJc;		//!< cross deriv of m_k with J and c
	vector< vector< vector<double> > > m_dkdcc;	// 2nd deriv of m_k with c
};

//-----------------------------------------------------------------------------
//! Base class for reaction rates.

class FEReactionRate : public FEMaterial
{
public:
	virtual void Init() {}
	
	//! reaction rate at material point
	virtual double ReactionRate(FEMaterialPoint& pt) = 0;
	
	//! tangent of reaction rate with strain at material point
	virtual mat3ds Tangent_ReactionRate_Strain(FEMaterialPoint& pt) = 0;
	
	//! tangent of reaction rate with effective fluid pressure at material point
	virtual double Tangent_ReactionRate_Pressure(FEMaterialPoint& pt) = 0;
	
public:
	FEChemicalReaction*	m_pReact;	//!< pointer to parent chemical reaction

};

//-----------------------------------------------------------------------------
//! Base class for chemical reactions.

typedef std::map<int,int> intmap;
typedef std::map<int,int>::iterator itrmap;

class FEChemicalReaction : public FEMaterial
{
public:
    FEChemicalReaction() {m_Vovr = false; m_pMP = 0; }
    
	virtual void Init() {
		if (m_pFwd) InitializeReactionRate(m_pFwd);
		if (m_pRev) InitializeReactionRate(m_pRev);
	}
	
	//! molar supply at material point
	virtual double ReactionSupply(FEMaterialPoint& pt) = 0;
	
	//! tangent of molar supply with strain at material point
	virtual mat3ds Tangent_ReactionSupply_Strain(FEMaterialPoint& pt) = 0;
	
	//! tangent of molar supply with effective pressure at material point
	virtual double Tangent_ReactionSupply_Pressure(FEMaterialPoint& pt) = 0;
	
	//! tangent of molar supply with effective concentration at material point
	virtual double Tangent_ReactionSupply_Concentration(FEMaterialPoint& pt, const int sol) = 0;
	
	// initialize chemical reaction rate
	void InitializeReactionRate(FEReactionRate* m_pRate) {m_pRate->m_pReact = this; }
	
public:
	intmap			m_solR;		//!< stoichiometric coefficients of solute reactants (input)
	intmap			m_solP;		//!< stoichiometric coefficients of solute products (input)
	intmap			m_sbmR;		//!< stoichiometric coefficients of solid-bound reactants (input)
	intmap			m_sbmP;		//!< stoichiometric coefficients of solid-bound products (input)
	
public:
	int				m_nsol;		//!< number of solutes in the mixture
	vector<int>		m_vR;		//!< stoichiometric coefficients of reactants
	vector<int>		m_vP;		//!< stoichiometric coefficients of products
	vector<int>		m_v;		//!< net stoichiometric coefficients of reactants and products
	FEReactionRate*	m_pFwd;		//!< pointer to forward reaction rate
	FEReactionRate*	m_pRev;		//!< pointer to reverse reaction rate
    double          m_Vbar;     //!< weighted molar volume of reactants and products
    bool            m_Vovr;     //!< override flag for m_Vbar
	
public:
	FEMultiphasic*	m_pMP;		//!< pointer to multiphasic material where reaction occurs
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

	//! partition coefficient
	double PartitionCoefficient(FEMaterialPoint& pt, const int sol);
	
	//! partition coefficients and their derivatives
	void PartitionCoefficientFunctions(FEMaterialPoint& mp, vector<double>& kappa,
									   vector<double>& dkdJ,
									   vector< vector<double> >& dkdc,
									   vector<double>& dkdJJ,
									   vector< vector<double> >& dkdJc,
									   vector< vector< vector<double> > >& dkdcc);
	
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
	double FluidTrueDensity() { return m_rhoTw; }
	
	//! solute density
	double SoluteTrueDensity(const int sol) { return m_pSolute[sol]->TrueDensity(); }
	
	//! solute molar mass
	double SoluteMolarMass(const int sol) { return m_pSolute[sol]->MolarMass(); }
	
	//! solute charge number
	int SoluteChargeNumber(const int sol) { return m_pSolute[sol]->ChargeNumber(); }
	
	//! SBM density
	double SBMTrueDensity(const int sbm) { return m_pSBM[sbm]->TrueDensity(); }
	
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
		return spt.m_sbmr[sbm]/SBMTrueDensity(sbm);
	}

	//! find local SBM ID from global one
	int FindLocalSBMID(int nid);
	
	//! Serialization
	void Serialize(DumpFile& ar);

	// initialize chemical reaction
	void InitializeReaction(FEChemicalReaction* m_pReact);

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
	FEElasticMaterial*				m_pSolid;		//!< pointer to elastic solid material
	FEHydraulicPermeability*		m_pPerm;		//!< pointer to permeability material
	FEOsmoticCoefficient*			m_pOsmC;		//!< pointer to osmotic coefficient material
	FESolventSupply*				m_pSupp;		//!< pointer to solvent supply material
	vector<FESolute*>				m_pSolute;		//!< pointer to solute materials
	vector<FESolidBoundMolecule*>	m_pSBM;			//!< pointer to solid-bound molecule materials
	vector<FEChemicalReaction*>		m_pReact;		//!< pointer to chemical reactions

	// declare as registered
	DECLARE_REGISTERED(FEMultiphasic);
	
	DECLARE_PARAMETER_LIST();
};
