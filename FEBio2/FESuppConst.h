#pragma once
#include "FEBiphasic.h"

//-----------------------------------------------------------------------------
// This class implements a material that has a constant solute supply

class FESuppConst :	public FESoluteSupply
{
public:
	//! constructor
	FESuppConst();
	
	//! diffusivity
	double Supply(FEMaterialPoint& pt);
	
	//! Tangent of supply with respect to strain
	double Tangent_Supply_Strain(FEMaterialPoint& mp);
	
	//! Tangent of supply with respect to concentration
	double Tangent_Supply_Concentration(FEMaterialPoint& mp);
	
	//! data initialization and checking
	void Init();
	
public:
	double	m_supp;			//!< diffusivity
	
	// declare as registered
	DECLARE_REGISTERED(FESuppConst);
	
	// declare parameter list
	DECLARE_PARAMETER_LIST();
};
