#include "FEOsmCoefConst.h"

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_PARAMETER_LIST(FEOsmCoefConst, FEOsmoticCoefficient)
	ADD_PARAMETER(m_osmcoef, FE_PARAM_DOUBLE, "osmcoef");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! Constructor. 
FEOsmCoefConst::FEOsmCoefConst(FEModel* pfem) : FEOsmoticCoefficient(pfem)
{
	m_osmcoef = 1;
}

//-----------------------------------------------------------------------------
//! Initialization. 
void FEOsmCoefConst::Init()
{
	FEOsmoticCoefficient::Init();
	if (m_osmcoef < 0) throw MaterialError("osmcoef must be >= 0");
}

//-----------------------------------------------------------------------------
//! Osmotic coefficient
double FEOsmCoefConst::OsmoticCoefficient(FEMaterialPoint& mp)
{
	// --- constant osmotic coefficient ---
	
	return m_osmcoef;
}

//-----------------------------------------------------------------------------
//! Tangent of osmotic coefficient with respect to strain
double FEOsmCoefConst::Tangent_OsmoticCoefficient_Strain(FEMaterialPoint &mp)
{
	return 0;
}

//-----------------------------------------------------------------------------
//! Tangent of osmotic coefficient with respect to concentration
double FEOsmCoefConst::Tangent_OsmoticCoefficient_Concentration(FEMaterialPoint &mp, const int isol)
{
	return 0;
}

