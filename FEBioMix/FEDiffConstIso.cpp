#include "stdafx.h"
#include "FEDiffConstIso.h"
#include <FECore/log.h>

// define the material parameters
BEGIN_FECORE_CLASS(FEDiffConstIso, FESoluteDiffusivity)
	ADD_PARAMETER(m_free_diff, FE_RANGE_GREATER         (0.0), "free_diff");
	ADD_PARAMETER(m_diff     , FE_RANGE_GREATER_OR_EQUAL(0.0), "diff"     );
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor. 
FEDiffConstIso::FEDiffConstIso(FEModel* pfem) : FESoluteDiffusivity(pfem)
{
	m_free_diff = m_diff = 1;
}

//-----------------------------------------------------------------------------
//! Validation
bool FEDiffConstIso::Validate()
{
	if (FESoluteDiffusivity::Validate() == false) return false;
	if (m_free_diff < m_diff) {
		feLogError("free_diff must be >= diff");
		return false;
	}
	return true;
}

//-----------------------------------------------------------------------------
//! Free diffusivity
double FEDiffConstIso::Free_Diffusivity(FEMaterialPoint& mp)
{
	return m_free_diff;
}

//-----------------------------------------------------------------------------
//! Tangent of free diffusivity with respect to concentration
double FEDiffConstIso::Tangent_Free_Diffusivity_Concentration(FEMaterialPoint& mp, const int isol)
{
	return 0;
}

//-----------------------------------------------------------------------------
//! Diffusivity tensor
mat3ds FEDiffConstIso::Diffusivity(FEMaterialPoint& mp)
{
	// --- constant isotropic diffusivity ---
	
	return mat3dd(m_diff);
}

//-----------------------------------------------------------------------------
//! Tangent of diffusivity with respect to strain
tens4ds FEDiffConstIso::Tangent_Diffusivity_Strain(FEMaterialPoint &mp)
{
	tens4ds D;
	D.zero();
	return D;
}

//-----------------------------------------------------------------------------
//! Tangent of diffusivity with respect to concentration
mat3ds FEDiffConstIso::Tangent_Diffusivity_Concentration(FEMaterialPoint &mp, const int isol)
{
	mat3ds d;
	d.zero();
	return d;
}
