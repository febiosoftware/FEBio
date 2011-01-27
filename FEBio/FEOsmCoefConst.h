#pragma once
#include "FEMaterial.h"

//-----------------------------------------------------------------------------
// This class implements a material that has a constant osmotic coefficient

class FEOsmCoefConst :	public FEOsmoticCoefficient
	{
	public:
		//! constructor
		FEOsmCoefConst();
		
		//! osmotic coefficient
		double OsmoticCoefficient(FEMaterialPoint& pt);
		
		//! Tangent of osmotic coefficient with respect to strain (J=detF)
		double Tangent_OsmoticCoefficient_Strain(FEMaterialPoint& mp);
		
		//! Tangent of osmotic coefficient with respect to concentration
		double Tangent_OsmoticCoefficient_Concentration(FEMaterialPoint& mp);
		
		//! data initialization and checking
		void Init();
		
	public:
		double	m_osmcoef;			//!< osmotic coefficient
		
		// declare as registered
		DECLARE_REGISTERED(FEOsmCoefConst);
		
		// declare parameter list
		DECLARE_PARAMETER_LIST();
	};
