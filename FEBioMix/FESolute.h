/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, The Trustees of Columbia University in 
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/



#pragma once
#include <FECore/FEMaterial.h>
#include <FECore/FEGlobalData.h>
#include <FECore/tens4d.h>
#include "febiomix_api.h"

//-----------------------------------------------------------------------------
//! Base class for solute diffusivity.
//! These materials need to define the diffusivity and tangent diffusivity functions.
//!
class FEBIOMIX_API FESoluteDiffusivity : public FEMaterial
{
public:
	//! constructor
	FESoluteDiffusivity(FEModel* pfem) : FEMaterial(pfem) {}

	//! solute diffusivity
	virtual mat3ds Diffusivity(FEMaterialPoint& pt) = 0;
	
	//! tangent of diffusivity with respect to strain
	virtual tens4dmm Tangent_Diffusivity_Strain(FEMaterialPoint& mp) = 0;
	
	//! tangent of diffusivity with respect to solute concentration
	virtual mat3ds Tangent_Diffusivity_Concentration(FEMaterialPoint& mp, const int isol) = 0;
	
	//! solute diffusivity in free solution
	virtual double Free_Diffusivity(FEMaterialPoint& pt) = 0;
	
	//! tangent of free diffusivity with respect to solute concentration
	virtual double Tangent_Free_Diffusivity_Concentration(FEMaterialPoint& pt, const int isol) = 0;
	
	//! set solute ID
	void SetSoluteID(const int ID) {m_ID = ID;}
	
	//! set solute ID
	int GetSoluteID() { return m_ID;}
	
private:
	int	m_ID;		//!< solute ID
	
};

//-----------------------------------------------------------------------------
//! Base class for solute solubility.
//! These materials need to define the solubility and tangent solubility functions.
//!
class FEBIOMIX_API FESoluteSolubility : public FEMaterial
{
public:
	//! constructor
	FESoluteSolubility(FEModel* pfem) : FEMaterial(pfem) {}

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
	
	//! set solute ID
	int GetSoluteID() { return m_ID;}
	
private:
	int	m_ID;		//!< solute ID
	
};

//-----------------------------------------------------------------------------
//! Base class for solute supply.
//! These materials need to define the solute supply and tangent supply functions.
//! The solute supply has units of moles/(referential mixture volume)/time
//!
class FEBIOMIX_API FESoluteSupply : public FEMaterial
{
public:
	//! constructor
	FESoluteSupply(FEModel* pfem) : FEMaterial(pfem) {}

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
//! Global solute data
//! This structure uniquely identifies a solute in multiphasic problems
class FESoluteData : public FEGlobalData 
{
public:
	FESoluteData(FEModel* pfem);

	//! initialization
	bool Init() override;

public:
	double	m_rhoT;			//!< true solute density
	double	m_M;			//!< solute molecular weight
	int		m_z;			//!< solute charge number

	DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
//! Base class for solute materials.

class FESolute : public FEMaterial
{
public:
	FESolute(FEModel* pfem);

public:
	bool Init() override;
	
	//! solute density
	double Density() { return m_rhoT; }
	
	//! solute molecular weight
	double MolarMass() { return m_M; }
	
	//! solute charge number
	int ChargeNumber() { return m_z; }
	
	//! Serialization
	void Serialize(DumpStream& ar) override;

	//! set solute ID
	void SetSoluteID(const int ID) {m_ID = ID;}
	
	//! get solute ID
	int GetSoluteID() {return m_ID;}

	//! set solute local ID
	void SetSoluteLocalID(const int LID) {m_LID = LID;}
	
	//! get solute local ID
	int GetSoluteLocalID() {return m_LID;}
  
	//! return the solute's dof
	int GetSoluteDOF() const { return m_ID - 1; }

private:
	FESoluteData* FindSoluteData(int nid);
	
private:
	int						m_ID;		//!< solute ID in global table
    int                     m_LID;      //!< solute local ID in parent material
	
public: // material parameters
	double					m_rhoT;		//!< true solute density
	double					m_M;		//!< solute molecular weight
	int						m_z;		//!< charge number of solute

public: // material properties
	FESoluteDiffusivity*	m_pDiff;	//!< pointer to diffusivity material
	FESoluteSolubility*		m_pSolub;	//!< pointer to solubility material
	FESoluteSupply*			m_pSupp;	//!< pointer to solute supply material

	DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
//! Global solid-bound molecule (SBM) data.
class FESBMData : public FEGlobalData
{
public:
	FESBMData(FEModel* pfem);

public:
	double	m_rhoT;			//!< SBM true density
	double	m_M;			//!< SBM molar mass
	int		m_z;			//!< SBM charge number

	DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
//! Base class for solid-bound molecules.

class FESolidBoundMolecule : public FEMaterial
{
public:
	FESolidBoundMolecule(FEModel* pfem);
	
public:
	bool Init() override;
	
	//! solute density
	double Density() { return m_rhoT; }
	
	//! solute molecular weight
	double MolarMass() { return m_M; }
	
	//! solute charge number
	int ChargeNumber() { return m_z; }
	
	//! Serialization
	void Serialize(DumpStream& ar) override;
	
	//! set solute ID
	void SetSBMID(const int ID) {m_ID = ID;}
	
	//! get SBM ID
	int GetSBMID() {return m_ID;}
	
private:
	FESBMData* FindSBMData(int nid);

private:
	int						m_ID;		//!< SBM ID in global table
	
public:
	double					m_rhoT;		//!< true SBM density
	double					m_M;		//!< SBM molar mass
	int						m_z;		//!< charge number of SBM
	double					m_rho0;		//!< initial referential (apparent) density of SBM
	double					m_rhomin;	//!< minimum referential (apparent) density of SBM
	double					m_rhomax;	//!< maximum referential (apparent) density of SBM
	
	DECLARE_FECORE_CLASS();
};
