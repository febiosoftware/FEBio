#pragma once
#include "FEBiphasic.h"

//-----------------------------------------------------------------------------
// This class implements a material that has a constant solute supply

class FESupplyConst :	public FESoluteSupply
{
public:
	//! constructor
	FESupplyConst();
	
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
	
	//! data initialization and checking
	void Init();
	
public:
	double	m_supp;			//!< solute supply
	
	// declare as registered
	DECLARE_REGISTERED(FESupplyConst);
	
	// declare parameter list
	DECLARE_PARAMETER_LIST();
};
