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

class FEBIOMIX_API FESupplyMichaelisMenten :	public FESoluteSupply
{
public:
	//! constructor
	FESupplyMichaelisMenten(FEModel* pfem);
	
	//! Solute supply
	double Supply(FEMaterialPoint& pt) override;
	
	//! Tangent of supply with respect to strain
	double Tangent_Supply_Strain(FEMaterialPoint& mp) override;
	
	//! Tangent of supply with respect to concentration
	double Tangent_Supply_Concentration(FEMaterialPoint& mp) override;
	
	//! receptor-ligand complex supply
	double ReceptorLigandSupply(FEMaterialPoint& mp) override;
	
	//! Solute supply at steady-state
	double SupplySS(FEMaterialPoint& pt) override;
	
	//! receptor-ligand complex concentration at steady-state
	double ReceptorLigandConcentrationSS(FEMaterialPoint& mp) override;
	
	//! referential solid supply
	double SolidSupply(FEMaterialPoint& pt) override;
	
	//! referential solid volume fraction under steady-state conditions
	double SolidConcentrationSS(FEMaterialPoint& pt) override;
	
public:
	double	m_Vmax;			//!< maximum uptake rate
	double	m_Km;			//!< concentration at which half-maximum rate occurs
	
	// declare parameter list
	DECLARE_FECORE_CLASS();
};
