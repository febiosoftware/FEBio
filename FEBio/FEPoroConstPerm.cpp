#include "stdafx.h"
#include "FEPoroConstPerm.h"


// register the material with the framework
REGISTER_MATERIAL(FEPoroConstPerm, "poroelastic");

// define the material parameters
BEGIN_PARAMETER_LIST(FEPoroConstPerm, FEPoroElastic)
	ADD_PARAMETER(m_perm, FE_PARAM_DOUBLE, "perm");
	ADD_PARAMETER(m_permv[0], FE_PARAM_DOUBLE, "permx");
	ADD_PARAMETER(m_permv[1], FE_PARAM_DOUBLE, "permy");
	ADD_PARAMETER(m_permv[2], FE_PARAM_DOUBLE, "permz");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! Constructor. 
FEPoroConstPerm::FEPoroConstPerm()
{
	m_perm = 1;
	m_permv[0] = m_permv[1] = m_permv[2] = 1;
}

//-----------------------------------------------------------------------------
//! Fluid flux.

vec3d FEPoroConstPerm::Flux(FEMaterialPoint& mp)
{
	FEPoroElasticMaterialPoint& pt = *mp.ExtractData<FEPoroElasticMaterialPoint>();

	// pressure gradient
	vec3d gradp = pt.m_gradp;

	// fluid flux w = -k*grad(p)
	vec3d w;
	w.x = -m_perm*m_permv[0]*gradp.x;
	w.y = -m_perm*m_permv[1]*gradp.y;
	w.z = -m_perm*m_permv[2]*gradp.z;

	return w;
}

//-----------------------------------------------------------------------------
//! Permeability tensor.
void FEPoroConstPerm::Permeability(double k[3][3], FEMaterialPoint& mp)
{
	FEPoroElasticMaterialPoint& pt = *mp.ExtractData<FEPoroElasticMaterialPoint>();

	// --- constant isotropic permeability ---

	k[0][0] = m_perm*m_permv[0];
	k[1][1] = m_perm*m_permv[1];
	k[2][2] = m_perm*m_permv[2];
	k[0][1] = k[0][2] = 0;
	k[1][0] = k[1][2] = 0;
	k[2][0] = k[2][1] = 0;
}

//-----------------------------------------------------------------------------
//! Tangent of permeability
tens4ds FEPoroConstPerm::Tangent_Permeability(FEMaterialPoint &mp)
{
	tens4ds K;
	K.zero();
	return K;
}
