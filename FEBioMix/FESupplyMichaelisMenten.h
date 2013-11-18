/*
 *  FESupplyMichaelisMenten.h
 *  FEBioXCode
 *
 *  Created by Gerard Ateshian on 6/27/12.
 *
 */

#pragma once
#include "FEBiphasicSolute.h"

//-----------------------------------------------------------------------------
// This class implements a material that has a solute supply based on
// Michaelis-Menten kinetics

class FESupplyMichaelisMenten :	public FESoluteSupply
{
public:
	//! constructor
	FESupplyMichaelisMenten(FEModel* pfem);
	
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
	double	m_Vmax;			//!< maximum uptake rate
	double	m_Km;			//!< concentration at which half-maximum rate occurs
	
	// declare as registered
	DECLARE_REGISTERED(FESupplyMichaelisMenten);
	
	// declare parameter list
	DECLARE_PARAMETER_LIST();
};
