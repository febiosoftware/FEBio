#include "FEDiffConstIso.h"

// define the material parameters
BEGIN_PARAMETER_LIST(FEDiffConstIso, FESoluteDiffusivity)
ADD_PARAMETER(m_free_diff, FE_PARAM_DOUBLE, "free_diff");
ADD_PARAMETER(m_diff, FE_PARAM_DOUBLE, "diff");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! Constructor. 
FEDiffConstIso::FEDiffConstIso(FEModel* pfem) : FESoluteDiffusivity(pfem)
{
	m_free_diff = m_diff = 1;
}

//-----------------------------------------------------------------------------
//! Initialization. 
void FEDiffConstIso::Init()
{
	if (m_free_diff <= 0) throw MaterialError("free_diff must be > 0");
	if (m_diff < 0) throw MaterialError("diff must be >= 0");
	if (m_free_diff < m_diff) throw MaterialError("free_diff must be >= diff");
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
