#pragma once
#include "FEBiphasicSolute.h"

//-----------------------------------------------------------------------------
// This class implements a material that has a constant osmotic coefficient

class FEOsmCoefConst :	public FEOsmoticCoefficient
{
public:
	//! constructor
	FEOsmCoefConst(FEModel* pfem);
	
	//! osmotic coefficient
	double OsmoticCoefficient(FEMaterialPoint& pt) override;
	
	//! Tangent of osmotic coefficient with respect to strain (J=detF)
	double Tangent_OsmoticCoefficient_Strain(FEMaterialPoint& mp) override;
	
	//! Tangent of osmotic coefficient with respect to concentration
	double Tangent_OsmoticCoefficient_Concentration(FEMaterialPoint& mp, const int isol) override;
	
public:
	double	m_osmcoef;			//!< osmotic coefficient
	
	// declare parameter list
	DECLARE_PARAMETER_LIST();
};
