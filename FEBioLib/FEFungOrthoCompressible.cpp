#include "stdafx.h"
#include "FEFungOrthoCompressible.h"

// define the material parameters
BEGIN_PARAMETER_LIST(FEFungOrthoCompressible, FEElasticMaterial)
ADD_PARAMETER(E1, FE_PARAM_DOUBLE, "E1");
ADD_PARAMETER(E2, FE_PARAM_DOUBLE, "E2");
ADD_PARAMETER(E3, FE_PARAM_DOUBLE, "E3");
ADD_PARAMETER(G12, FE_PARAM_DOUBLE, "G12");
ADD_PARAMETER(G23, FE_PARAM_DOUBLE, "G23");
ADD_PARAMETER(G31, FE_PARAM_DOUBLE, "G31");
ADD_PARAMETER(v12, FE_PARAM_DOUBLE, "v12");
ADD_PARAMETER(v23, FE_PARAM_DOUBLE, "v23");
ADD_PARAMETER(v31, FE_PARAM_DOUBLE, "v31");
ADD_PARAMETER(m_c, FE_PARAM_DOUBLE, "c");
ADD_PARAMETER(m_k, FE_PARAM_DOUBLE, "k");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! Data initialization
void FEFungOrthoCompressible::Init()
{
	FEElasticMaterial::Init();
	
	if (E1 <= 0) throw MaterialError("E1 should be positive");
	if (E2 <= 0) throw MaterialError("E2 should be positive");
	if (E3 <= 0) throw MaterialError("E3 should be positive");
	
	if (v12 > sqrt(E1/E2)) throw MaterialError("Invalid value for v12. Let v12 > sqrt(E1/E2)");
	if (v23 > sqrt(E2/E3)) throw MaterialError("Invalid value for v23. Let v23 > sqrt(E2/E3)");
	if (v31 > sqrt(E3/E1)) throw MaterialError("Invalid value for v31. Let v31 > sqrt(E3/E1)");
	
	if (m_c <= 0) throw MaterialError("c should be positive");

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
//! Calculates the stress
mat3ds FEFungOrthoCompressible::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	int i,j;
	vec3d a0[3];		// texture direction in reference configuration
	vec3d a[3];			// texture direction in current configuration
	mat3ds A[3];		// texture tensor in current configuration
	double K[3];		// Ka
	double L[3];		// La
	double eQ;			// exp(Q)
	mat3ds bmi;			// B - I
	// Evaluate the deformation gradient
	mat3d &F = pt.F;
	double J = pt.J;
	
	// calculate left and right Cauchy-Green tensor
	mat3ds b = pt.LeftCauchyGreen();
	mat3ds c = pt.RightCauchyGreen();
	mat3ds c2 = c*c;
	mat3dd I(1.);
	
	for (i=0; i<3; i++) {	// Perform sum over all three texture directions
		// Copy the texture direction in the reference configuration to a0
		a0[i].x = pt.Q[0][i]; a0[i].y = pt.Q[1][i]; a0[i].z = pt.Q[2][i];
		K[i] = a0[i]*(c*a0[i]);
		L[i] = a0[i]*(c2*a0[i]);
		a[i] = F*a0[i]/sqrt(K[i]);	// Evaluate the texture direction in the current configuration
		A[i] = dyad(a[i]);			// Evaluate the texture tensor in the current configuration
	}
	
	// Evaluate exp(Q)
	eQ = 0;
	for (i=0; i<3; i++) {
		eQ += 2*mu[i]*(L[i]-2*K[i]+1);
		for (j=0; j<3; j++)
			eQ += lam[i][j]*(K[i]-1)*(K[j]-1);
	}
	eQ = exp(eQ/(4.*m_c));
	
	// Evaluate the stress
	mat3ds s;
	s.zero();		// Initialize for summation
	bmi = b - I;
	for (i=0; i<3; i++) {
		s += mu[i]*K[i]*(A[i]*bmi + bmi*A[i]);
		for (j=0; j<3; j++)
			s += lam[i][j]*((K[i]-1)*K[j]*A[j]+(K[j]-1)*K[i]*A[i])/2.;
	}
	s *= eQ/(2.0*J);
	
	s += I*(m_k*log(J)/J);
	
	return s;
}

//-----------------------------------------------------------------------------
//! Calculates the tangent
tens4ds FEFungOrthoCompressible::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	int i,j;
	vec3d a0[3];		// texture direction in reference configuration
	vec3d a[3];			// texture direction in current configuration
	mat3ds A[3];		// texture tensor in current configuration
	double K[3];		// Ka
	double L[3];		// La
	double eQ;			// exp(Q)
	mat3ds bmi;			// B - I
	// Evaluate the strain and texture
	mat3d &F = pt.F;
	double J = pt.J;
	
	// calculate left and right Cauchy-Green tensor
	mat3ds b = pt.LeftCauchyGreen();
	mat3ds c = pt.RightCauchyGreen();
	mat3ds c2 = c*c;
	mat3dd I(1.);
	
	for (i=0; i<3; i++) {	// Perform sum over all three texture directions
		// Copy the texture direction in the reference configuration to a0
		a0[i].x = pt.Q[0][i]; a0[i].y = pt.Q[1][i]; a0[i].z = pt.Q[2][i];
		K[i] = a0[i]*(c*a0[i]);
		L[i] = a0[i]*(c2*a0[i]);
		a[i] = F*a0[i]/sqrt(K[i]);	// Evaluate the texture direction in the current configuration
		A[i] = dyad(a[i]);			// Evaluate the texture tensor in the current configuration
	}
	
	// Evaluate exp(Q)
	eQ = 0;
	for (i=0; i<3; i++) {
		eQ += 2*mu[i]*(L[i]-2*K[i]+1);
		for (j=0; j<3; j++)
			eQ += lam[i][j]*(K[i]-1)*(K[j]-1);
	}
	eQ = exp(eQ/(4.*m_c));
	
	// Evaluate the Cauchy stress
	mat3ds s;
	s.zero();		// Initialize for summation
	//	tens4ds C;
	//	C.zero();
	tens4ds C(0.0);
	bmi = b - I;
	for (i=0; i<3; i++) {
		s += mu[i]*K[i]*(A[i]*bmi + bmi*A[i]);
		for (j=0; j<3; j++)
			s += lam[i][j]*((K[i]-1)*K[j]*A[j]+(K[j]-1)*K[i]*A[i])/2.;
		C += mu[i]*K[i]*dyad4s(A[i],b);
		for (j=0; j<3; j++)
			C += lam[i][j]*K[i]*K[j]*dyad1s(A[i],A[j])/2.;
	}
	s *= eQ/(2.0*J);
	
	// Elasticity tensor
	C = (eQ/J)*C + (2*J/m_c/eQ)*dyad1s(s);
	C += (dyad1s(I) - dyad4s(I)*(2*log(J)))*(m_k/J);
	
	return C;
}
