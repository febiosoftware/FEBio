#include "stdafx.h"
#include "FEViscoElasticMaterial.h"

// register the material with the framework
REGISTER_MATERIAL(FEViscoElasticMaterial, "viscoelastic");

// define the material parameters
BEGIN_PARAMETER_LIST(FEViscoElasticMaterial, FEMaterial)
	ADD_PARAMETER(m_t[0], FE_PARAM_DOUBLE, "t1");
	ADD_PARAMETER(m_t[1], FE_PARAM_DOUBLE, "t2");
	ADD_PARAMETER(m_t[2], FE_PARAM_DOUBLE, "t3");
	ADD_PARAMETER(m_t[3], FE_PARAM_DOUBLE, "t4");
	ADD_PARAMETER(m_t[4], FE_PARAM_DOUBLE, "t5");
	ADD_PARAMETER(m_t[5], FE_PARAM_DOUBLE, "t6");
	ADD_PARAMETER(m_g[0], FE_PARAM_DOUBLE, "g1");
	ADD_PARAMETER(m_g[1], FE_PARAM_DOUBLE, "g2");
	ADD_PARAMETER(m_g[2], FE_PARAM_DOUBLE, "g3");
	ADD_PARAMETER(m_g[3], FE_PARAM_DOUBLE, "g4");
	ADD_PARAMETER(m_g[4], FE_PARAM_DOUBLE, "g5");
	ADD_PARAMETER(m_g[5], FE_PARAM_DOUBLE, "g6");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! constructor
FEViscoElasticMaterial::FEViscoElasticMaterial() : FEMaterial(FE_VISCO_ELASTIC)
{
	m_pemat = 0;
	m_g0 = 1;
	for (int i=0; i<MAX_TERMS; ++i)
	{
		m_t[i] = 1;
		m_g[i] = 0;
	}
}

//-----------------------------------------------------------------------------
//! Stress function

mat3ds FEViscoElasticMaterial::Stress(FEMaterialPoint& mp)
{
	// get the elastic part
	FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();

	// get the viscoelsatic point data
	FEViscoElasticMaterialPoint& pt = *mp.ExtractData<FEViscoElasticMaterialPoint>();

	// get the new elastic stress
	pt.m_se = m_pemat->Stress(mp);

	double dt = mp.dt;

	double g, h;

	int i;
	m_g0 = 1;
//	for (i=0; i<MAX_TERMS; ++i) m_g0 -= m_g[i];

	// calculate the new visco-elastic stress
	mat3ds s = pt.m_se*m_g0;

	// "Cauchy" stress of previous timestep
	mat3ds sp = ep.push_forward(pt.m_Sep);

	for (i=0; i<MAX_TERMS; ++i)
	{
		g = exp(-dt/m_t[i]);
		h = (1 - g)/(dt/m_t[i]);

		pt.m_H[i] = pt.m_Hp[i]*g + (pt.m_se - sp)*h;
		s += pt.m_H[i]*m_g[i];
	}

	return s;
}

//-----------------------------------------------------------------------------

void FEViscoElasticMaterial::Tangent(double D[6][6], FEMaterialPoint& pt)
{
	// calculate the elastic tangent
	m_pemat->Tangent(D, pt);

	double dt = pt.dt;

	int i;
	m_g0 = 1;
//	for (i=0; i<MAX_TERMS; ++i) m_g0 -= m_g[i];

	// multiply with visco-factor
	double f = m_g0, g, h;
	for (i=0; i<MAX_TERMS; ++i)
	{
		g = exp(-dt/m_t[i]);
		h = ( 1 - exp(-dt/m_t[i]) )/( dt/m_t[i] );
		f += m_g[i]*h; 
	}

	D[0][0] *= f; D[0][1] *= f; D[0][2] *= f; D[0][3] *= f; D[0][4] *= f; D[0][5] *= f;
	D[1][0] *= f; D[1][1] *= f; D[1][2] *= f; D[1][3] *= f; D[1][4] *= f; D[1][5] *= f;
	D[2][0] *= f; D[2][1] *= f; D[2][2] *= f; D[2][3] *= f; D[2][4] *= f; D[2][5] *= f;
	D[3][0] *= f; D[3][1] *= f; D[3][2] *= f; D[3][3] *= f; D[3][4] *= f; D[3][5] *= f;
	D[4][0] *= f; D[4][1] *= f; D[4][2] *= f; D[4][3] *= f; D[4][4] *= f; D[4][5] *= f;
	D[5][0] *= f; D[5][1] *= f; D[5][2] *= f; D[5][3] *= f; D[5][4] *= f; D[5][5] *= f;
}
