// FEStVenantKirchhoff.cpp: implementation of the FEStVenantKirchhoff class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FEStVenantKirchhoff.h"

// register the material with the framework
REGISTER_MATERIAL(FEStVenantKirchhoff, "St.Venant-Kirchhoff");

// define the material parameters
BEGIN_PARAMETER_LIST(FEStVenantKirchhoff, FEElasticMaterial)
	ADD_PARAMETER(m_E, FE_PARAM_DOUBLE, "E");
	ADD_PARAMETER(m_v, FE_PARAM_DOUBLE, "v");
END_PARAMETER_LIST();

//////////////////////////////////////////////////////////////////////
// FEStVenantKirchhoff
//////////////////////////////////////////////////////////////////////

void FEStVenantKirchhoff::Init()
{
	// intialize base class
	FEElasticMaterial::Init();

	if (m_E <= 0) throw MaterialError("Invalid value for E");
	if (!INRANGE(m_v, -1.0, 0.5)) throw MaterialError("Invalid value for v");
}

mat3ds FEStVenantKirchhoff::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	mat3d &F = pt.F;
	double detF = pt.J;
	double detFi = 1.0 / detF;

	double b[3][3], b2[3][3];

	double trE;

	// lame parameters
	double lam = m_v*m_E/((1+m_v)*(1-2*m_v));
	double mu  = 0.5*m_E/(1+m_v);

	// calculate left Cauchy-Green tensor (ie. b-matrix)
	b[0][0] = F[0][0]*F[0][0]+F[0][1]*F[0][1]+F[0][2]*F[0][2];
	b[0][1] = F[0][0]*F[1][0]+F[0][1]*F[1][1]+F[0][2]*F[1][2];
	b[0][2] = F[0][0]*F[2][0]+F[0][1]*F[2][1]+F[0][2]*F[2][2];

	b[1][0] = F[1][0]*F[0][0]+F[1][1]*F[0][1]+F[1][2]*F[0][2];
	b[1][1] = F[1][0]*F[1][0]+F[1][1]*F[1][1]+F[1][2]*F[1][2];
	b[1][2] = F[1][0]*F[2][0]+F[1][1]*F[2][1]+F[1][2]*F[2][2];

	b[2][0] = F[2][0]*F[0][0]+F[2][1]*F[0][1]+F[2][2]*F[0][2];
	b[2][1] = F[2][0]*F[1][0]+F[2][1]*F[1][1]+F[2][2]*F[1][2];
	b[2][2] = F[2][0]*F[2][0]+F[2][1]*F[2][1]+F[2][2]*F[2][2];

	// calculate trace of Green-Lagrance strain tensor
	trE = 0.5*(b[0][0]+b[1][1]+b[2][2]-3);

	// calculate square of b-matrix
	// (we commented out the matrix components we do not need)
	b2[0][0] = b[0][0]*b[0][0]+b[0][1]*b[1][0]+b[0][2]*b[2][0];
	b2[0][1] = b[0][0]*b[0][1]+b[0][1]*b[1][1]+b[0][2]*b[2][1];
	b2[0][2] = b[0][0]*b[0][2]+b[0][1]*b[1][2]+b[0][2]*b[2][2];

//	b2[1][0] = b[1][0]*b[0][0]+b[1][1]*b[1][0]+b[1][2]*b[2][0];
	b2[1][1] = b[1][0]*b[0][1]+b[1][1]*b[1][1]+b[1][2]*b[2][1];
	b2[1][2] = b[1][0]*b[0][2]+b[1][1]*b[1][2]+b[1][2]*b[2][2];

//	b2[2][0] = b[2][0]*b[0][0]+b[2][1]*b[1][0]+b[2][2]*b[2][0];
//	b2[2][1] = b[2][0]*b[0][1]+b[2][1]*b[1][1]+b[2][2]*b[2][1];
	b2[2][2] = b[2][0]*b[0][2]+b[2][1]*b[1][2]+b[2][2]*b[2][2];

	// calculate stress
	mat3ds s;

	s.xx() = (lam*trE*b[0][0] + mu*(b2[0][0] - b[0][0]))*detFi;
	s.yy() = (lam*trE*b[1][1] + mu*(b2[1][1] - b[1][1]))*detFi;
	s.zz() = (lam*trE*b[2][2] + mu*(b2[2][2] - b[2][2]))*detFi;
	s.xy() = (lam*trE*b[0][1] + mu*(b2[0][1] - b[0][1]))*detFi;
	s.yz() = (lam*trE*b[1][2] + mu*(b2[1][2] - b[1][2]))*detFi;
	s.xz() = (lam*trE*b[0][2] + mu*(b2[0][2] - b[0][2]))*detFi;

	return s;
}

void FEStVenantKirchhoff::Tangent(double D[6][6], FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// deformation gradient
	mat3d& F = pt.F;
	double detF = pt.J;

	// left cauchy-green matrix (i.e. the 'b' matrix)
	// (we commented out the matrix components we do not need)
	double b[3][3];
	b[0][0] = F[0][0]*F[0][0]+F[0][1]*F[0][1]+F[0][2]*F[0][2];
	b[0][1] = F[0][0]*F[1][0]+F[0][1]*F[1][1]+F[0][2]*F[1][2];
	b[0][2] = F[0][0]*F[2][0]+F[0][1]*F[2][1]+F[0][2]*F[2][2];

//	b[1][0] = F[1][0]*F[0][0]+F[1][1]*F[0][1]+F[1][2]*F[0][2];
	b[1][1] = F[1][0]*F[1][0]+F[1][1]*F[1][1]+F[1][2]*F[1][2];
	b[1][2] = F[1][0]*F[2][0]+F[1][1]*F[2][1]+F[1][2]*F[2][2];

//	b[2][0] = F[2][0]*F[0][0]+F[2][1]*F[0][1]+F[2][2]*F[0][2];
//	b[2][1] = F[2][0]*F[1][0]+F[2][1]*F[1][1]+F[2][2]*F[1][2];
	b[2][2] = F[2][0]*F[2][0]+F[2][1]*F[2][1]+F[2][2]*F[2][2];

	// lame parameters
	double lam = m_v*m_E/((1+m_v)*(1-2*m_v));
	double mu  = 0.5*m_E/(1+m_v);

	double lam1 = lam / detF;
	double mu1  = 2.0*mu/detF;

	D[0][0] = lam1*b[0][0]*b[0][0] + mu1*(b[0][0]*b[0][0]);
	D[1][1] = lam1*b[1][1]*b[1][1] + mu1*(b[1][1]*b[1][1]);
	D[2][2] = lam1*b[2][2]*b[2][2] + mu1*(b[2][2]*b[2][2]);

	D[3][3] = lam1*b[0][1]*b[0][1] + mu1*0.5*(b[0][0]*b[1][1]+b[0][1]*b[0][1]);
	D[4][4] = lam1*b[1][2]*b[1][2] + mu1*0.5*(b[1][1]*b[2][2]+b[1][2]*b[1][2]);
	D[5][5] = lam1*b[0][2]*b[0][2] + mu1*0.5*(b[0][0]*b[2][2]+b[0][2]*b[0][2]);

	D[1][0] = D[0][1] = lam1*b[0][0]*b[1][1] + mu1*b[0][1]*b[0][1];
	D[2][0] = D[0][2] = lam1*b[0][0]*b[2][2] + mu1*b[0][2]*b[0][2];
	D[2][1] = D[1][2] = lam1*b[1][1]*b[2][2] + mu1*b[1][2]*b[1][2];

	D[3][0] = D[0][3] = lam1*b[0][0]*b[0][1] + mu1*b[0][0]*b[0][1];
	D[4][0] = D[0][4] = lam1*b[0][0]*b[1][2] + mu1*b[0][1]*b[0][2];
	D[5][0] = D[0][5] = lam1*b[0][0]*b[0][2] + mu1*b[0][0]*b[0][2];

	D[3][1] = D[1][3] = lam1*b[1][1]*b[0][1] + mu1*b[0][1]*b[1][1];
	D[4][1] = D[1][4] = lam1*b[1][1]*b[1][2] + mu1*b[1][1]*b[1][2];
	D[5][1] = D[1][5] = lam1*b[1][1]*b[0][2] + mu1*b[0][1]*b[1][2];

	D[3][2] = D[2][3] = lam1*b[2][2]*b[0][1] + mu1*b[0][2]*b[1][2];
	D[4][2] = D[2][4] = lam1*b[2][2]*b[1][2] + mu1*b[1][2]*b[2][2];
	D[5][2] = D[2][5] = lam1*b[2][2]*b[0][2] + mu1*b[0][2]*b[2][2];

	D[4][3] = D[3][4] = lam1*b[0][1]*b[1][2] + mu1*0.5*(b[0][1]*b[1][2] + b[0][2]*b[1][1]);
	D[5][3] = D[3][5] = lam1*b[0][1]*b[0][2] + mu1*0.5*(b[0][0]*b[1][2] + b[0][2]*b[0][1]);
	D[5][4] = D[4][5] = lam1*b[1][2]*b[0][2] + mu1*0.5*(b[0][1]*b[2][2] + b[1][2]*b[0][2]);

}
