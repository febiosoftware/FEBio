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

	if (E1 <= 0) throw MaterialError("E1 should be positive");
	if (E2 <= 0) throw MaterialError("E2 should be positive");
	if (E3 <= 0) throw MaterialError("E3 should be positive");
	
	if (G12 < 0) throw MaterialError("G12 should be positive");
	if (G23 < 0) throw MaterialError("G23 should be positive");
	if (G31 < 0) throw MaterialError("G31 should be positive");
	
	if (v12 > sqrt(E1/E2)) throw MaterialError("Invalid value for v12. Let v12 <= sqrt(E1/E2)");
	if (v23 > sqrt(E2/E3)) throw MaterialError("Invalid value for v23. Let v23 <= sqrt(E2/E3)");
	if (v31 > sqrt(E3/E1)) throw MaterialError("Invalid value for v31. Let v31 <= sqrt(E3/E1)");
	
	// Evaluate Lame coefficients
	mu[0] = G12 + G31 - G23;
	mu[1] = G12 - G31 + G23;
	mu[2] =-G12 + G31 + G23;
	lam[0][0] = 1.0/E1; lam[0][1] = -v12/E1; lam[0][2] = -v31/E3;
	lam[1][0] = -v12/E1; lam[1][1] = 1.0/E2; lam[1][2] = -v23/E2;
	lam[2][0] = -v31/E3; lam[2][1] = -v23/E2; lam[2][2] = 1.0/E3;
	mat3d c(lam);
	c = c.inverse();
	lam[0][0] = c[0][0] - 2*mu[0];
	lam[1][1] = c[1][1] - 2*mu[1];
	lam[2][2] = c[2][2] - 2*mu[2];
	lam[1][2] = c[1][2]; lam[2][1] = c[2][1];
	lam[2][0] = c[2][0]; lam[0][2] = c[0][2];
	lam[0][1] = c[0][1]; lam[1][0] = c[1][0];
}

//-----------------------------------------------------------------------------
//! Calculates the stress for a linear orthotropic material
mat3ds FELinearOrthotropic::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	int i,j;
	mat3ds A0[3];		// structural tensor in reference configuration
	double K[3];		// Ka

	// Evaluate the strain
	mat3ds E=pt.Strain();
	
	for (i=0; i<3; i++) {	// Perform sum over all three texture directions
		// Copy the texture direction in the reference configuration to a0
		vec3d a0(pt.Q[0][i],pt.Q[1][i],pt.Q[2][i]);
		K[i] = a0*(E*a0);
		A0[i] = dyad(a0);			// Evaluate the texture tensor in the reference configuration
	}
	
	// Evaluate the stress
	mat3ds s;
	s.zero();		// Initialize for summation
	for (i=0; i<3; i++) {
		s += mu[i]*(A0[i]*E + E*A0[i]);
		for (j=0; j<3; j++)
			s += (A0[j]*K[i]+A0[i]*K[j])*(lam[i][j]/2.);
	}

	return s;
}

tens4ds FELinearOrthotropic::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	int i,j;
	mat3ds A0[3];		// texture tensor in current configuration
	
	for (i=0; i<3; i++) {	// Perform sum over all three texture directions
		// Copy the texture direction in the reference configuration to a0
		vec3d a0(pt.Q[0][i],pt.Q[1][i],pt.Q[2][i]);
		A0[i] = dyad(a0);			// Evaluate the texture tensor in the reference configuration
	}
	
	// Evaluate the elasticity tensor
	tens4ds C(0.0);
	mat3dd I(1.0);
	for (i=0; i<3; i++) {
		C += dyad4s(A0[i],I)*mu[i];
		for (j=0; j<3; j++)
			C += dyad1s(A0[i],A0[j])*(lam[i][j]/2.);
	}
	
	return C;
}
