#include "FEDiffConstOrtho.h"

// define the material parameters
BEGIN_FECORE_CLASS(FEDiffConstOrtho, FESoluteDiffusivity)
	ADD_PARAMETER(m_free_diff, FE_RANGE_GREATER(0.0), "free_diff");
	ADD_PARAMETER(m_diff     , 3, FE_RANGE_GREATER_OR_EQUAL(0.0), "diff" );
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor. 
FEDiffConstOrtho::FEDiffConstOrtho(FEModel* pfem) : FESoluteDiffusivity(pfem)
{
	m_free_diff = m_diff[0] = m_diff[1] = m_diff[2] = 1;
}

//-----------------------------------------------------------------------------
//! Initialization. 
bool FEDiffConstOrtho::Validate()
{
	if (FESoluteDiffusivity::Validate() == false) return false;
	if (m_free_diff < m_diff[0]) return fecore_error("free_diff must be >= diff1");
	if (m_free_diff < m_diff[1]) return fecore_error("free_diff must be >= diff2");
	if (m_free_diff < m_diff[2]) return fecore_error("free_diff must be >= diff3");
	return true;
}

//-----------------------------------------------------------------------------
//! Free diffusivity
double FEDiffConstOrtho::Free_Diffusivity(FEMaterialPoint& mp)
{
	return m_free_diff;
}

//-----------------------------------------------------------------------------
//! Tangent of free diffusivity with respect to concentration
double FEDiffConstOrtho::Tangent_Free_Diffusivity_Concentration(FEMaterialPoint& mp, const int isol)
{
	return 0;
}

//-----------------------------------------------------------------------------
//! Diffusivity tensor
mat3ds FEDiffConstOrtho::Diffusivity(FEMaterialPoint& mp)
{
	vec3d a0;				// texture direction in reference configuration
	mat3ds d(0,0,0,0,0,0);	// diffusion tensor

	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// --- constant orthotropic diffusivity ---
	for (int i=0; i<3; i++) {	// Perform sum over all three texture directions
		
		// Copy the texture direction in the reference configuration to a0
		a0.x = pt.m_Q[0][i]; a0.y = pt.m_Q[1][i]; a0.z = pt.m_Q[2][i];
		
		// Evaluate the texture tensor in the current configuration
		d += dyad(a0)*m_diff[i];
	}
	
	return d;
}

//-----------------------------------------------------------------------------
//! Tangent of diffusivity with respect to strain
tens4ds FEDiffConstOrtho::Tangent_Diffusivity_Strain(FEMaterialPoint &mp)
{
	tens4ds D;
	D.zero();
	return D;
}

//-----------------------------------------------------------------------------
//! Tangent of diffusivity with respect to concentration
mat3ds FEDiffConstOrtho::Tangent_Diffusivity_Concentration(FEMaterialPoint &mp, const int isol)
{
	mat3ds d;
	d.zero();
	return d;
}
