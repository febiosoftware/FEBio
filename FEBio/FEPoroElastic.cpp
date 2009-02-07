// FEPoroElastic.cpp: implementation of the FEPoroElastic class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FEPoroElastic.h"

// register the material with the framework
REGISTER_MATERIAL(FEPoroElastic, "poroelastic");

// define the material parameters
BEGIN_PARAMETER_LIST(FEPoroElastic, FEMaterial)
	ADD_PARAMETER(m_perm, FE_PARAM_DOUBLE, "perm");
	ADD_PARAMETER(m_permv[0], FE_PARAM_DOUBLE, "permx");
	ADD_PARAMETER(m_permv[1], FE_PARAM_DOUBLE, "permy");
	ADD_PARAMETER(m_permv[2], FE_PARAM_DOUBLE, "permz");
END_PARAMETER_LIST();

//////////////////////////////////////////////////////////////////////
// FEPoroElastic
//////////////////////////////////////////////////////////////////////

FEPoroElastic::FEPoroElastic()
{
	m_perm = 1;
	m_permv[0] = m_permv[1] = m_permv[2] = 1;
}

mat3ds FEPoroElastic::Stress(FEMaterialPoint& mp)
{
	FEPoroElasticMaterialPoint& pt = *mp.ExtractData<FEPoroElasticMaterialPoint>();

	// calculate solid material stress
	mat3ds s = m_psmat->Stress(mp);

	// add fluid pressure
	s.xx() -= pt.m_p;
	s.yy() -= pt.m_p;
	s.zz() -= pt.m_p;

	return s;
}

void FEPoroElastic::Tangent(double D[6][6], FEMaterialPoint& mp)
{
	FEPoroElasticMaterialPoint& pt = *mp.ExtractData<FEPoroElasticMaterialPoint>();

	// call solid tangent routine
	m_psmat->Tangent(D, mp);

	// fluid pressure
	double p = pt.m_p;

	// adjust tangent for pressures
	D[0][0] -= -p;
	D[1][1] -= -p;
	D[2][2] -= -p;

	D[0][1] -= p; D[1][0] -= p;
	D[1][2] -= p; D[2][1] -= p;
	D[0][2] -= p; D[2][0] -= p;

	D[3][3] -= -p;
	D[4][4] -= -p;
	D[5][5] -= -p;
}

vec3d FEPoroElastic::Flux(FEMaterialPoint& mp)
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

void FEPoroElastic::Permeability(double k[3][3], FEMaterialPoint& mp)
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
