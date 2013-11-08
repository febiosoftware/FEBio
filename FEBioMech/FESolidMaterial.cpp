#include "stdafx.h"
#include "FESolidMaterial.h"

//-----------------------------------------------------------------------------
// Derivative of stress w.r.t. solute concentration at material point
// Set this to zero by default because elasticity problems do not require it
mat3ds FESolidMaterial::Tangent_Concentration(FEMaterialPoint& pt, const int isol)
{
	return mat3ds(0,0,0,0,0,0);
}
