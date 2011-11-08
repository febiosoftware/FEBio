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
//! Check material parameters.
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
//! Calculates the stress for a linear orthotropic material. It calls the 
//! FElinearOrthotropic::Tangent function and contracts it with the small
//! strain tensor.
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
//! Calculates the elasticity tensor for an orthotropic material. The elasticity
//! tangent is first calculated in the local material coordinate system and then
//! rotated to the global coordinate system.
//! \todo come up with some verification problems for this material model
tens4ds FELinearOrthotropic::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// get the transformation tensor
	mat3d& Q = pt.Q;

	// set-up the elasticity tensor in the local coordinate system
	double v21 = E1*v12/E2;
	double v32 = E2*v23/E3;
	double v13 = E3*v31/E1;

	double d = (1.0 - v12*v21 - v23*v32 - v13*v31 - 2.0*v12*v23*v31)/(E1*E2*E3);

	double C[6][6] = {0};
	C[0][0] = (1 - v23*v32)/(E2*E3*d  ); C[0][1] = (v12 + v13*v32)/(E2*E3*d); C[0][2] = (v13 + v12*v23)/(E2*E3*d);
	C[1][0] = C[0][1]                  ; C[1][1] = (1.0 - v13*v31)/(E1*E3*d); C[1][2] = (v23 + v21*v13)/(E1*E3*d);
	C[2][0] = C[0][2]                  ; C[2][1] = C[1][2]                  ; C[2][2] = (1.0 - v12*v21)/(E1*E2*d);
	C[3][3] = G12; C[4][4] = G23; C[5][5] = G31;

	// index look-up tables (Voigt-to-Matrix and Matrix-to-Voigt)
	const int VTM[6][2] = {{0,0},{1,1},{2,2},{0,1},{1,2},{2,0}};
	const int MTV[3][3] = {{0,3,5},{3,1,4},{5,4,2}};

	// transform to the global coordinate system
	double D[6][6] = {0};
	for (int i=0; i<6; ++i)
		for (int j=0; j<6; ++j)
		{
			int r = VTM[i][0];
			int s = VTM[i][1];
			int t = VTM[j][0];
			int u = VTM[j][1];
			for (int k=0; k<3; ++k)
				for (int l=0; l<3; ++l)
					for (int m=0; m<3; ++m)
						for (int n=0; n<3; ++n)
						{
							int I = MTV[k][l];
							int J = MTV[m][n];
							D[i][j] += Q[r][k]*Q[s][l]*Q[t][m]*Q[u][n]*C[I][J];
						}
		}

	// return result
	return tens4ds(D);
}
