#pragma once
#include "FECore/FEMaterial.h"

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
	virtual mat3ds Tangent_Diffusivity_Concentration(FEMaterialPoint& mp, const int isol) = 0;
	
	//! solute diffusivity in free solution
	virtual double Free_Diffusivity(FEMaterialPoint& pt) = 0;
	
	//! tangent of free diffusivity with respect to solute concentration
	virtual double Tangent_Free_Diffusivity_Concentration(FEMaterialPoint& pt, const int isol) = 0;
	
	//! set solute ID
	void SetSoluteID(const int ID) {m_ID = ID;}
	
private:
	int	m_ID;		//!< solute ID
	
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
	virtual double Tangent_Solubility_Concentration(FEMaterialPoint& mp, const int isol) = 0;
	
	//! cross derivative of solubility with respect to strain and concentration
	virtual double Tangent_Solubility_Strain_Concentration(FEMaterialPoint& mp, const int isol) = 0;
	
	//! second derivative of solubility with respect to strain
	virtual double Tangent_Solubility_Strain_Strain(FEMaterialPoint& mp) = 0;
	
	//! second derivative of solubility with respect to concentration
	virtual double Tangent_Solubility_Concentration_Concentration(FEMaterialPoint& mp, 
																  const int isol, const int jsol) = 0;
	
	//! set solute ID
	void SetSoluteID(const int ID) {m_ID = ID;}
	
private:
	int	m_ID;		//!< solute ID
	
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
//! Base class for solute materials.

class FESolute : public FEMultiMaterial
{
public:
	FESolute();

public:
	void Init();
	
	//! solute density
	double Density() { return m_rhoT; }
	
	//! solute molecular weight
	double MolarMass() { return m_M; }
	
	//! solute charge number
	int ChargeNumber() { return m_z; }
	
	//! Serialization
	void Serialize(DumpFile& ar);

	//! set solute ID
	void SetSoluteID(const int ID) {m_ID = ID;}
	
	//! get solute ID
	int GetSoluteID() {return m_ID;}

	//! Find a material parameter
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

	//! set the material attribute
	bool SetAttribute(const char* szname, const char* szval);
	
private:
	int						m_ID;		//!< solute ID in global table
	
public: // material parameters
	double					m_rhoT;		//!< true solute density
	double					m_M;		//!< solute molecular weight
	int						m_z;		//!< charge number of solute

public: // material properties
	FESoluteDiffusivity*	m_pDiff;	//!< pointer to diffusivity material
	FESoluteSolubility*		m_pSolub;	//!< pointer to solubility material
	FESoluteSupply*			m_pSupp;	//!< pointer to solute supply material
	
	DECLARE_REGISTERED(FESolute);
};

//-----------------------------------------------------------------------------
//! Base class for solid-bound molecules.

class FESolidBoundMolecule : public FEMaterial
{
public:
	FESolidBoundMolecule();
	
public:
	void Init();
	
	//! solute density
	double Density() { return m_rhoT; }
	
	//! solute molecular weight
	double MolarMass() { return m_M; }
	
	//! solute charge number
	int ChargeNumber() { return m_z; }
	
	//! Serialization
	void Serialize(DumpFile& ar);
	
	//! set solute ID
	void SetSBMID(const int ID) {m_ID = ID;}
	
	//! get SBM ID
	int GetSBMID() {return m_ID;}
	
public:
	//! set the material attribute
	bool SetAttribute(const char* szname, const char* szval);

	//! Find a material parameter
	FEParam* GetParameter(const ParamString& s);
	
private:
	int						m_ID;		//!< SBM ID in global table
	
public:
	double					m_rhoT;		//!< true SBM density
	double					m_M;		//!< SBM molar mass
	int						m_z;		//!< charge number of SBM
	double					m_rho0;		//!< initial referential (apparent) density of SBM
	double					m_rhomin;	//!< minimum referential (apparent) density of SBM
	double					m_rhomax;	//!< maximum referential (apparent) density of SBM
	
	// declare as registered
	DECLARE_REGISTERED(FESolidBoundMolecule);
	
	DECLARE_PARAMETER_LIST();
};
