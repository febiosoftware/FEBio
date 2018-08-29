#pragma once
#include "FEBiphasicSolute.h"

//-----------------------------------------------------------------------------
// This class implements a material that has a constant solute supply

class FESupplyConst :	public FESoluteSupply
{
public:
	//! constructor
	FESupplyConst(FEModel* pfem);
	
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
	double	m_supp;			//!< solute supply
	
	// declare parameter list
	DECLARE_PARAMETER_LIST();
};
