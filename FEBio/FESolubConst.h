#pragma once
#include "FEMaterial.h"

//-----------------------------------------------------------------------------
// This class implements a material that has a constant solute solubility

class FESolubConst : public FESoluteSolubility
	{
	public:
		//! constructor
		FESolubConst();
		
		//! solubility
		double Solubility(FEMaterialPoint& pt);
		
		//! Tangent of solubility with respect to strain
		double Tangent_Solubility_Strain(FEMaterialPoint& mp);
		
		//! Tangent of solubility with respect to concentration
		double Tangent_Solubility_Concentration(FEMaterialPoint& mp);
		
		//! Cross derivative of solubility with respect to strain and concentration
		double Tangent_Solubility_Strain_Concentration(FEMaterialPoint& mp);
		
		//! Second derivative of solubility with respect to strain
		double Tangent_Solubility_Strain_Strain(FEMaterialPoint& mp);
		
		//! data initialization and checking
		void Init();
		
	public:
		double	m_solub;			//!< solubility
		
		// declare as registered
		DECLARE_REGISTERED(FESolubConst);
		
		// declare parameter list
		DECLARE_PARAMETER_LIST();
	};
