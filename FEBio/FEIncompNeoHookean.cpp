// FEIncompNeoHookean.cpp: implementation of the FEIncompNeoHookean class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FEIncompNeoHookean.h"

// register the material with the framework
REGISTER_MATERIAL(FEIncompNeoHookean, "incomp neo-Hookean");

// define the material parameters
BEGIN_PARAMETER_LIST(FEIncompNeoHookean, FEIncompressibleMaterial)
	ADD_PARAMETER(m_G, FE_PARAM_DOUBLE, "G");
END_PARAMETER_LIST();

//////////////////////////////////////////////////////////////////////
// IncompNeoHookean
//////////////////////////////////////////////////////////////////////

mat3ds FEIncompNeoHookean::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	double muJ, Ib;

	mat3d &F = pt.F;
	double detF = pt.J;

	double p = pt.avgp;

	// calculate left Cauchy-Green tensor
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

	// first invariant of b = trace of b
	Ib = b[0][0]+b[1][1]+b[2][2];

	// calculate stress
	mat3ds s;

	muJ = m_G/pow(detF, 5.0/3.0);

	s.xx() = muJ*(b[0][0] - Ib/3.) + p;
	s.yy() = muJ*(b[1][1] - Ib/3.) + p;
	s.zz() = muJ*(b[2][2] - Ib/3.) + p;
	s.xy() = muJ*b[0][1];
	s.yz() = muJ*b[1][2];
	s.xz() = muJ*b[0][2];

	return s;
}

void FEIncompNeoHookean::Tangent(double D[6][6], FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// deformation gradient
	mat3d &F = pt.F;
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

	// trace of b
	double Ib = b[0][0]+b[1][1]+b[2][2];

	double muJ = m_G/pow(detF, 5.0/3.0);

	double p = pt.avgp; // average element pressure

	D[0][0] = 2.*muJ*(4.*Ib/9. - 2.*b[0][0]/3.) - p;
	D[1][1] = 2.*muJ*(4.*Ib/9. - 2.*b[1][1]/3.) - p;
	D[2][2] = 2.*muJ*(4.*Ib/9. - 2.*b[2][2]/3.) - p;

	D[1][0] = D[0][1] = 2.*muJ*(Ib/9. - (1./3.)*b[0][0] - (1./3.)*b[1][1]) + p;
	D[2][0] = D[0][2] = 2.*muJ*(Ib/9. - (1./3.)*b[0][0] - (1./3.)*b[2][2]) + p;
	D[2][1] = D[1][2] = 2.*muJ*(Ib/9. - (1./3.)*b[1][1] - (1./3.)*b[2][2]) + p;

	D[3][0] = D[0][3] = D[3][1] = D[1][3] = D[3][2] = D[2][3] = -2.*muJ/3.*b[0][1];
	D[4][0] = D[0][4] = D[4][1] = D[1][4] = D[4][2] = D[2][4] = -2.*muJ/3.*b[1][2];
	D[5][0] = D[0][5] = D[5][1] = D[1][5] = D[5][2] = D[2][5] = -2.*muJ/3.*b[0][2];

	D[3][3] = D[4][4] = D[5][5] = muJ*Ib/3. - p;
}
