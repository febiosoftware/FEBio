/*
 *  FESuppConst.cpp
 *
 */

#include "stdafx.h"
#include "FESuppConst.h"

// register the material with the framework
REGISTER_MATERIAL(FESuppConst, "supp-const");

// define the material parameters
BEGIN_PARAMETER_LIST(FESuppConst, FESoluteSupply)
	ADD_PARAMETER(m_supp, FE_PARAM_DOUBLE, "supp");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! Constructor. 
FESuppConst::FESuppConst()
{
	m_supp = 0;
}

//-----------------------------------------------------------------------------
//! Initialization. 
void FESuppConst::Init()
{
}

//-----------------------------------------------------------------------------
//! Solute supply
double FESuppConst::Supply(FEMaterialPoint& mp)
{
	// --- constant solubility ---
	
	return m_supp;
}

//-----------------------------------------------------------------------------
//! Tangent of solute supply with respect to strain
double FESuppConst::Tangent_Supply_Strain(FEMaterialPoint &mp)
{
	return 0;
}

//-----------------------------------------------------------------------------
//! Tangent of solute supply with respect to concentration
double FESuppConst::Tangent_Supply_Concentration(FEMaterialPoint &mp)
{
	return 0;
}

