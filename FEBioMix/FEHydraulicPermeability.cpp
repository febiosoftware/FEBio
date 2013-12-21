#include "FEHydraulicPermeability.h"

//-----------------------------------------------------------------------------
// Material parameters for FEHydraulicPermeability
void FEHydraulicPermeability::Init()
{
	FEMaterial::Init();
}

//-----------------------------------------------------------------------------
// Derivative of permeability w.r.t. solute concentration at material point
// Set this to zero by default because poroelasticity problems do not require it
mat3ds FEHydraulicPermeability::Tangent_Permeability_Concentration(FEMaterialPoint& pt, const int isol)
{
	return mat3ds(0,0,0,0,0,0);
}

