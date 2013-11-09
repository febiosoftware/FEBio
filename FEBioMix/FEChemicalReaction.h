#pragma once
#include "FECore/FEMaterial.h"
#include <map>

//-----------------------------------------------------------------------------
class FEMultiphasic;
class FEChemicalReaction;

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
	//! constructor
    FEChemicalReaction();
    
	//! initialization
	virtual void Init();	

	// initialize chemical reaction rate
	void InitializeReactionRate(FEReactionRate* m_pRate);

public:
	//! return number of properties
	int Properties();

	//! return a material property
	FEMaterial* GetProperty(int i);

	//! find a material property index ( returns <0 for error)
	int FindPropertyIndex(const char* szname);

	//! set a material property (returns false on error)
	bool SetProperty(int i, FEMaterial* pm);

	void SetParameter(FEParam& p);

	bool SetParameterAttribute(FEParam& p, const char* szatt, const char* szval);

public:
	//! set the forward reaction rate
	void SetForwardReactionRate(FEReactionRate* pfwd) { m_pFwd = pfwd; }

	//! set the reverse reaction rate
	void SetReverseReactionRate(FEReactionRate* prev) { m_pRev = prev; }

	//! set the stoichiometric coefficients of solid-bound reactions
	void SetSolidReactantsCoefficients(int id, int vR) { m_sbmR.insert(std::pair<int, int>(id, vR)); }

	//! set the stoichiometric coefficients of solid-bound produts
	void SetSolidProductsCoefficients(int id, int vR) { m_sbmP.insert(std::pair<int, int>(id, vR)); }

	//! set the stoichiometric coefficients of solute reactions
	void SetSoluteReactantsCoefficients(int id, int vR) { m_solR.insert(std::pair<int, int>(id, vR)); }

	//! set the stoichiometric coefficients of solute products
	void SetSoluteProductsCoefficients(int id, int vR) { m_solP.insert(std::pair<int, int>(id, vR)); }

public:

	//! molar supply at material point
	virtual double ReactionSupply(FEMaterialPoint& pt) = 0;
	
	//! tangent of molar supply with strain at material point
	virtual mat3ds Tangent_ReactionSupply_Strain(FEMaterialPoint& pt) = 0;
	
	//! tangent of molar supply with effective pressure at material point
	virtual double Tangent_ReactionSupply_Pressure(FEMaterialPoint& pt) = 0;
	
	//! tangent of molar supply with effective concentration at material point
	virtual double Tangent_ReactionSupply_Concentration(FEMaterialPoint& pt, const int sol) = 0;

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
    double          m_Vbar;     //!< weighted molar volume of reactants and products
    bool            m_Vovr;     //!< override flag for m_Vbar
	int				m_vRtmp;	//!< helper variable for reading in stoichiometric coefficients for reactants
	int				m_vPtmp;	//!< helper variable for reading in stoichiometric coefficients for products

	FEReactionRate*	m_pFwd;		//!< pointer to forward reaction rate
	FEReactionRate*	m_pRev;		//!< pointer to reverse reaction rate

public:
	FEMultiphasic*	m_pMP;		//!< pointer to multiphasic material where reaction occurs

	DECLARE_PARAMETER_LIST();
};
