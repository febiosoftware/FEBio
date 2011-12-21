/*
 *  FESupplySynthesisBinding.h
 *
 *  Created by Gerard Ateshian on 8/18/11.
 *
 */

#pragma once
#include "FEBiphasicSolute.h"

//-----------------------------------------------------------------------------
// This class implements a material that has a solute supply based on
// receptor-ligand binding kinetics as described by the Langmuir or Hill equation
// and also includes a constant supply rate (synthesis when positive, degradation
// when negative).

class FESupplySynthesisBinding : public FESoluteSupply
{
public:
	//! constructor
	FESupplySynthesisBinding();
	
	//! Solute supply
	double Supply(FEMaterialPoint& pt);
	
	//! Tangent of supply with respect to strain
	double Tangent_Supply_Strain(FEMaterialPoint& mp);
	
	//! Tangent of supply with respect to concentration
	double Tangent_Supply_Concentration(FEMaterialPoint& mp);
	
	//! receptor-ligand complex supply
	double ReceptorLigandSupply(FEMaterialPoint& mp);
	
	//! Solute supply at steady-state
	double SupplySS(FEMaterialPoint& pt);
	
	//! receptor-ligand complex concentration at steady-state
	double ReceptorLigandConcentrationSS(FEMaterialPoint& mp);

	//! referential solid supply
	double SolidSupply(FEMaterialPoint& pt);
	
	//! referential solid volume fraction under steady-state conditions
	double SolidConcentrationSS(FEMaterialPoint& pt);
	
	//! data initialization and checking
	void Init();
	
public:
	double	m_supp;			//!< synthesis rate
	double	m_kf;			//!< forward reaction rate constant
	double	m_kr;			//!< reverse reaction rate constant
	double	m_crt;			//!< total receptor concentration (referential)
	
	// declare parameter list
	DECLARE_PARAMETER_LIST();
};
