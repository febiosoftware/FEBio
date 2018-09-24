#include "FEPermConstIso.h"

// define the material parameters
BEGIN_PARAMETER_LIST(FEPermConstIso, FEHydraulicPermeability)
	ADD_PARAMETER(m_perm, FE_RANGE_GREATER_OR_EQUAL(0.0), "perm");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! Constructor. 
FEPermConstIso::FEPermConstIso(FEModel* pfem) : FEHydraulicPermeability(pfem)
{
	m_perm = 1;
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
