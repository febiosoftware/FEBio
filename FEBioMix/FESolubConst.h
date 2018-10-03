#pragma once
#include "FEBiphasicSolute.h"

//-----------------------------------------------------------------------------
// This class implements a material that has a constant solute solubility

class FESolubConst : public FESoluteSolubility
{
public:
	//! constructor
	FESolubConst(FEModel* pfem);
	
	//! solubility
	double Solubility(FEMaterialPoint& pt) override;
	
	//! Tangent of solubility with respect to strain
	double Tangent_Solubility_Strain(FEMaterialPoint& mp) override;
	
	//! Tangent of solubility with respect to concentration
	double Tangent_Solubility_Concentration(FEMaterialPoint& mp, const int isol) override;
	
	//! Cross derivative of solubility with respect to strain and concentration
	double Tangent_Solubility_Strain_Concentration(FEMaterialPoint& mp, const int isol) override;
	
	//! Second derivative of solubility with respect to strain
	double Tangent_Solubility_Strain_Strain(FEMaterialPoint& mp) override;

	//! Second derivative of solubility with respect to concentration
	double Tangent_Solubility_Concentration_Concentration(FEMaterialPoint& mp, const int isol, const int jsol) override;

public:
	double	m_solub;			//!< solubility
	
	// declare parameter list
	DECLARE_FECORE_CLASS();
};
