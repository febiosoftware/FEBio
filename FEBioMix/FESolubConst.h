#pragma once
#include "FEBiphasicSolute.h"

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
	double Tangent_Solubility_Concentration(FEMaterialPoint& mp, const int isol);
	
	//! Cross derivative of solubility with respect to strain and concentration
	double Tangent_Solubility_Strain_Concentration(FEMaterialPoint& mp, const int isol);
	
	//! Second derivative of solubility with respect to strain
	double Tangent_Solubility_Strain_Strain(FEMaterialPoint& mp);

	//! Second derivative of solubility with respect to concentration
	double Tangent_Solubility_Concentration_Concentration(FEMaterialPoint& mp, const int isol, const int jsol);

	//! data initialization and checking
	void Init();
	
public:
	double	m_solub;			//!< solubility
	
	// declare parameter list
	DECLARE_PARAMETER_LIST();
};
