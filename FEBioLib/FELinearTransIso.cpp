#include "stdafx.h"
#include "FELinearTransIso.h"

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_PARAMETER_LIST(FELinearTransIso, FEElasticMaterial)
	ADD_PARAMETER(E1, FE_PARAM_DOUBLE, "E1");
	ADD_PARAMETER(E3, FE_PARAM_DOUBLE, "E3");
	ADD_PARAMETER(G23, FE_PARAM_DOUBLE, "G23");
	ADD_PARAMETER(v12, FE_PARAM_DOUBLE, "v12");
	ADD_PARAMETER(v31, FE_PARAM_DOUBLE, "v31");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! Check material parameters.
//! \todo check material parameters
void FELinearTransIso::Init()
{
	FEElasticMaterial::Init();
}

//-----------------------------------------------------------------------------
//! Calculates the stress for a linear orthotropic material. It calls the 
//! FElinearOrthotropic::Tangent function and contracts it with the small
//! strain tensor.
mat3ds FELinearTransIso::Stress(FEMaterialPoint& mp)
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
tens4ds FELinearTransIso::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// get the transformation tensor
	mat3d& Q = pt.m_Q;

	// set-up the elasticity tensor in the local coordinate system
	double v13 = E3*v31/E1;
	double G12 = 0.5*E1/(1.0 + v12);

	double d = (1.0 + v12)*(1.0 - v12 - 2.0*v13*v31)/(E1*E1*E3);

	double C[6][6] = {0};
	C[0][0] = (1 - v13*v31)/(E1*E3*d  ); C[0][1] = (v12 + v13*v31)/(E1*E3*d); C[0][2] = v13*(1.0 + v12)/(E1*E3*d);
	C[1][0] = C[0][1]                  ; C[1][1] = C[0][0]                  ; C[1][2] = C[0][2];
	C[2][0] = C[0][2]                  ; C[2][1] = C[1][2]                  ; C[2][2] = (1.0 - v12*v12)/(E1*E1*d);
	C[3][3] = G12; C[4][4] = G23; C[5][5] = G23;

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
