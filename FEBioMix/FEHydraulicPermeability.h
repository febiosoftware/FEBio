#pragma once
#include "FECore/FEMaterial.h"

//-----------------------------------------------------------------------------
//! Base class for hydraulic permeability of porous materials.
//! These materials need to define the permeability and tangent permeability functions.
//!
class FEHydraulicPermeability : public FEMaterial
{
public:
	FEHydraulicPermeability(FEModel* pfem) : FEMaterial(pfem) {}
	virtual ~FEHydraulicPermeability(){}
    
	//! hydraulic permeability
	virtual mat3ds Permeability(FEMaterialPoint& pt) = 0;
    
	//! tangent of hydraulic permeability with respect to strain
	virtual tens4ds Tangent_Permeability_Strain(FEMaterialPoint& mp) = 0;
    
	//! tangent of hydraulic permeability with respect to concentration
	mat3ds Tangent_Permeability_Concentration(FEMaterialPoint& mp, const int isol);
    
	void Init();
};

