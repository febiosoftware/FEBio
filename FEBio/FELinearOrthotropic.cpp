#include "stdafx.h"
#include "FELinearOrthotropic.h"

// register the material with the framework
REGISTER_MATERIAL(FELinearOrthotropic, "linear orthotropic");

// define the material parameters
BEGIN_PARAMETER_LIST(FELinearOrthotropic, FEElasticMaterial)
	ADD_PARAMETER(E1, FE_PARAM_DOUBLE, "E1");
	ADD_PARAMETER(E2, FE_PARAM_DOUBLE, "E2");
	ADD_PARAMETER(E3, FE_PARAM_DOUBLE, "E3");
	ADD_PARAMETER(G12, FE_PARAM_DOUBLE, "G12");
	ADD_PARAMETER(G23, FE_PARAM_DOUBLE, "G23");
	ADD_PARAMETER(G31, FE_PARAM_DOUBLE, "G31");
	ADD_PARAMETER(v12, FE_PARAM_DOUBLE, "v12");
	ADD_PARAMETER(v23, FE_PARAM_DOUBLE, "v23");
	ADD_PARAMETER(v31, FE_PARAM_DOUBLE, "v31");
END_PARAMETER_LIST();

void FELinearOrthotropic::Init()
{
	FEElasticMaterial::Init();

	if (E1 <= 0) throw MaterialError("Invalid value for E1");
	if (E2 <= 0) throw MaterialError("Invalid value for E2");
	if (E3 <= 0) throw MaterialError("Invalid value for E3");

	if (v12 > sqrt(E1/E2)) throw MaterialError("Invalid value for v12");
	if (v23 > sqrt(E2/E3)) throw MaterialError("Invalid value for v23");
	if (v31 > sqrt(E3/E1)) throw MaterialError("Invalid value for v31");
}

//-----------------------------------------------------------------------------
//! Calculates the stress for a linear orthotropic material
mat3ds FELinearOrthotropic::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// deformation gradient
	mat3d &F = pt.F;

	// caculate small strain tensor
	mat3ds e;
	e.xx() = F[0][0] - 1.0;
	e.yy() = F[1][1] - 1.0;
	e.zz() = F[2][2] - 1.0;
	e.xy() = 0.5*(F[0][1] + F[1][0]);
	e.xz() = 0.5*(F[0][2] + F[2][0]);
	e.yz() = 0.5*(F[1][2] + F[2][1]);

	// calculate stiffness matrix
	tens4ds C = Tangent(pt);

	double D[6][6] = {0};
	C.extract(D);

	// calculate stress
	mat3ds s;

	s.xx() = D[0][0]*e.xx() + D[0][1]*e.yy() + D[0][2]*e.zz();
	s.yy() = D[1][0]*e.xx() + D[1][1]*e.yy() + D[1][2]*e.zz();
	s.zz() = D[2][0]*e.xx() + D[2][1]*e.yy() + D[2][2]*e.zz();
	s.yz() = D[3][3]*e.yz();
	s.xz() = D[4][4]*e.xz();
	s.xy() = D[5][5]*e.xy();

	return s;
}

tens4ds FELinearOrthotropic::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	double v21 = E2/E1*v12;
	double v32 = E3/E2*v23;
	double v13 = E1/E3*v31;

	double G = 1.0/(1.0 - v12*v21 - v23*v32 - v31*v13 - 2.0*v21*v32*v13);

	double D[6][6] = {0};
	D[0][0] = E1*(1.0 - v23*v32)*G; 
	D[1][1] = E2*(1.0 - v13*v31)*G;
	D[2][2] = E3*(1.0 - v12*v21)*G;
	D[3][3] = G23;
	D[4][4] = G31;
	D[5][5] = G12;

	D[0][1] = D[1][0] = E1*(v21 + v31*v23)*G;
	D[0][2] = D[2][0] = E1*(v31 + v21*v32)*G;
	D[1][2] = D[2][1] = E2*(v32 + v12*v31)*G;

	return tens4ds(D);
}
