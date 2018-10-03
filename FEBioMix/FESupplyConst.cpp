/*
 *  FESupplyConst.cpp
 *
 */

#include "FESupplyConst.h"

// define the material parameters
BEGIN_FECORE_CLASS(FESupplyConst, FESoluteSupply)
	ADD_PARAMETER(m_supp, "supp");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor. 
FESupplyConst::FESupplyConst(FEModel* pfem) : FESoluteSupply(pfem)
{
	m_supp = 0;
}

//-----------------------------------------------------------------------------
//! Solute supply
double FESupplyConst::Supply(FEMaterialPoint& mp)
{
	// --- constant solubility ---
	
	return m_supp;
}

//-----------------------------------------------------------------------------
//! Tangent of solute supply with respect to strain
double FESupplyConst::Tangent_Supply_Strain(FEMaterialPoint &mp)
{
	return 0;
}

//-----------------------------------------------------------------------------
//! Tangent of solute supply with respect to concentration
double FESupplyConst::Tangent_Supply_Concentration(FEMaterialPoint &mp)
{
	return 0;
}

//-----------------------------------------------------------------------------
//! Receptor-ligand complex supply
double FESupplyConst::ReceptorLigandSupply(FEMaterialPoint &mp)
{
	return 0;
}

//-----------------------------------------------------------------------------
//! Solute supply at steady-state
double FESupplyConst::SupplySS(FEMaterialPoint& mp)
{
	// --- constant solubility ---
	
	return m_supp;
}

//-----------------------------------------------------------------------------
//! Receptor-ligand concentration at steady-state
double FESupplyConst::ReceptorLigandConcentrationSS(FEMaterialPoint& mp)
{
	return 0;
}

//-----------------------------------------------------------------------------
//! Referential solid supply (moles of solid/referential volume/time)
double FESupplyConst::SolidSupply(FEMaterialPoint& mp)
{
	return ReceptorLigandSupply(mp);
}

//-----------------------------------------------------------------------------
//! Referential solid concentration (moles of solid/referential volume)
//! at steady-state
double FESupplyConst::SolidConcentrationSS(FEMaterialPoint& mp)
{
	return 0;
}


