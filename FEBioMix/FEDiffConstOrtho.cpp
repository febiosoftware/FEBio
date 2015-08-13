#include "FEDiffConstOrtho.h"

// define the material parameters
BEGIN_PARAMETER_LIST(FEDiffConstOrtho, FESoluteDiffusivity)
	ADD_PARAMETER2(m_free_diff, FE_PARAM_DOUBLE, FE_RANGE_GREATER(0.0), "free_diff");
	ADD_PARAMETERV(m_diff , FE_PARAM_DOUBLEV, 3, "diff" );
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! Constructor. 
FEDiffConstOrtho::FEDiffConstOrtho(FEModel* pfem) : FESoluteDiffusivity(pfem)
{
	m_free_diff = m_diff[0] = m_diff[1] = m_diff[2] = 1;
}

//-----------------------------------------------------------------------------
//! Initialization. 
void FEDiffConstOrtho::Init()
{
	FESoluteDiffusivity::Init();

	if (m_diff[0] < 0) throw MaterialError("diff1 must be >= 0");
	if (m_diff[1] < 0) throw MaterialError("diff2 must be >= 0");
	if (m_diff[2] < 0) throw MaterialError("diff3 must be >= 0");
	if (m_free_diff < m_diff[0]) throw MaterialError("free_diff must be >= diff1");
	if (m_free_diff < m_diff[1]) throw MaterialError("free_diff must be >= diff2");
	if (m_free_diff < m_diff[2]) throw MaterialError("free_diff must be >= diff3");
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
