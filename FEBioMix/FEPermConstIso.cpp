#include "FEPermConstIso.h"

// define the material parameters
BEGIN_PARAMETER_LIST(FEPermConstIso, FEHydraulicPermeability)
	ADD_PARAMETER(m_perm, FE_PARAM_DOUBLE, "perm");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! Constructor. 
FEPermConstIso::FEPermConstIso(FEModel* pfem) : FEHydraulicPermeability(pfem)
{
	m_perm = 1;
}

//-----------------------------------------------------------------------------
//! Initialization. 
void FEPermConstIso::Init()
{
	FEHydraulicPermeability::Init();
	if (m_perm < 0) throw MaterialError("perm must be >= 0");
}

//-----------------------------------------------------------------------------
//! Permeability tensor.
mat3ds FEPermConstIso::Permeability(FEMaterialPoint& mp)
{
	// --- constant isotropic permeability ---
	
	return mat3dd(m_perm);
}

//-----------------------------------------------------------------------------
//! Tangent of permeability
tens4ds FEPermConstIso::Tangent_Permeability_Strain(FEMaterialPoint &mp)
{
	tens4ds K;
	K.zero();
	return K;
}
