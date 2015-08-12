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

public:
	double	m_supp;			//!< solute supply
	
	// declare parameter list
	DECLARE_PARAMETER_LIST();
};
