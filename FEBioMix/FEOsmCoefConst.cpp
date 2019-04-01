#include "stdafx.h"
#include "FEOsmCoefConst.h"

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FEOsmCoefConst, FEOsmoticCoefficient)
	ADD_PARAMETER(m_osmcoef, FE_RANGE_GREATER_OR_EQUAL(0.0), "osmcoef");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor. 
FEOsmCoefConst::FEOsmCoefConst(FEModel* pfem) : FEOsmoticCoefficient(pfem)
{
	m_osmcoef = 1;
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

