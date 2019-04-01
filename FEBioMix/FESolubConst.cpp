#include "stdafx.h"
#include "FESolubConst.h"

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FESolubConst, FESoluteSolubility)
	ADD_PARAMETER(m_solub, FE_RANGE_GREATER_OR_EQUAL(0.0), "solub");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor. 
FESolubConst::FESolubConst(FEModel* pfem) : FESoluteSolubility(pfem)
{
	m_solub = 1;
}

//-----------------------------------------------------------------------------
//! Solubility
double FESolubConst::Solubility(FEMaterialPoint& mp)
{
	// --- constant solubility ---
	
	return m_solub;
}

//-----------------------------------------------------------------------------
//! Tangent of solubility with respect to strain
double FESolubConst::Tangent_Solubility_Strain(FEMaterialPoint &mp)
{
	return 0;
}

//-----------------------------------------------------------------------------
//! Tangent of solubility with respect to concentration
double FESolubConst::Tangent_Solubility_Concentration(FEMaterialPoint &mp, const int isol)
{
	return 0;
}

//-----------------------------------------------------------------------------
//! Cross derivative of solubility with respect to strain and concentration
double FESolubConst::Tangent_Solubility_Strain_Concentration(FEMaterialPoint &mp, const int isol)
{
	return 0;
}

//-----------------------------------------------------------------------------
//! Second derivative of solubility with respect to strain
double FESolubConst::Tangent_Solubility_Strain_Strain(FEMaterialPoint &mp)
{
	return 0;
}

//-----------------------------------------------------------------------------
//! Second derivative of solubility with respect to concentration
double FESolubConst::Tangent_Solubility_Concentration_Concentration(FEMaterialPoint &mp, const int isol, const int jsol)
{
	return 0;
}
