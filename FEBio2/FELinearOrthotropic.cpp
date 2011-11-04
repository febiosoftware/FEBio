#include "stdafx.h"
#include "FELinearOrthotropic.h"

//-----------------------------------------------------------------------------
// register the material with the framework
REGISTER_MATERIAL(FELinearOrthotropic, "linear orthotropic");

//-----------------------------------------------------------------------------
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

//-----------------------------------------------------------------------------
void FELinearOrthotropic::Init()
{
	FEElasticMaterial::Init();

	if (E1 <= 0) throw MaterialError("E1 should be positive");
	if (E2 <= 0) throw MaterialError("E2 should be positive");
	if (E3 <= 0) throw MaterialError("E3 should be positive");
	
	if (G12 < 0) throw MaterialError("G12 should be positive");
	if (G23 < 0) throw MaterialError("G23 should be positive");
	if (G31 < 0) throw MaterialError("G31 should be positive");
	
	if (v12 > sqrt(E1/E2)) throw MaterialError("Invalid value for v12. Let v12 <= sqrt(E1/E2)");
	if (v23 > sqrt(E2/E3)) throw MaterialError("Invalid value for v23. Let v23 <= sqrt(E2/E3)");
	if (v31 > sqrt(E3/E1)) throw MaterialError("Invalid value for v31. Let v31 <= sqrt(E3/E1)");
}

//-----------------------------------------------------------------------------
//! Calculates the stress for a linear orthotropic material
//! \todo implement material orientation
mat3ds FELinearOrthotropic::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// Evaluate the small-strain tensor
	mat3ds e = pt.SmallStrain();

	// get the tangent
	tens4ds C = Tangent(mp);
	
	// stress = C:e
	return C.dot(e);
}

//-----------------------------------------------------------------------------
//! \todo implement material orientation
tens4ds FELinearOrthotropic::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	double v21 = E1*v12/E2;
	double v32 = E2*v23/E3;
	double v13 = E3*v31/E1;

	double d = (1.0 - v12*v21 - v23*v32 - v13*v31 - 2.0*v12*v23*v31)/(E1*E2*E3);

	double C[6][6] = {0};
	C[0][0] = (1 - v23*v32)/(E2*E3*d  ); C[0][1] = (v12 + v13*v32)/(E2*E3*d); C[0][2] = (v13 + v12*v23)/(E2*E3*d);
	C[1][0] = C[0][1]                  ; C[1][1] = (1.0 - v13*v31)/(E1*E3*d); C[1][2] = (v23 + v21*v13)/(E1*E3*d);
	C[2][0] = C[0][2]                  ; C[2][1] = C[1][2]                  ; C[2][2] = (1.0 - v12*v21)/(E1*E2*d);
	C[3][3] = G12; C[4][4] = G23; C[5][5] = G31;

	return tens4ds(C);
}
