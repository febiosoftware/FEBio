#include "stdafx.h"
#include "FEOgdenUnconstrained.h"

BEGIN_PARAMETER_LIST(FEOgdenUnconstrained, FEElasticMaterial);
	ADD_PARAMETER(m_p   , FE_PARAM_DOUBLE, "cp");
	ADD_PARAMETER(m_c[0], FE_PARAM_DOUBLE, "c1");
	ADD_PARAMETER(m_c[1], FE_PARAM_DOUBLE, "c2");
	ADD_PARAMETER(m_c[2], FE_PARAM_DOUBLE, "c3");
	ADD_PARAMETER(m_c[3], FE_PARAM_DOUBLE, "c4");
	ADD_PARAMETER(m_c[4], FE_PARAM_DOUBLE, "c5");
	ADD_PARAMETER(m_c[5], FE_PARAM_DOUBLE, "c6");
	ADD_PARAMETER(m_m[0], FE_PARAM_DOUBLE, "m1");
	ADD_PARAMETER(m_m[1], FE_PARAM_DOUBLE, "m2");
	ADD_PARAMETER(m_m[2], FE_PARAM_DOUBLE, "m3");
	ADD_PARAMETER(m_m[3], FE_PARAM_DOUBLE, "m4");
	ADD_PARAMETER(m_m[4], FE_PARAM_DOUBLE, "m5");
	ADD_PARAMETER(m_m[5], FE_PARAM_DOUBLE, "m6");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! constructor
FEOgdenUnconstrained::FEOgdenUnconstrained()
{
	for (int i=0; i<MAX_TERMS; ++i)
	{
		m_c[i] = 0;
		m_m[i] = 1;
	}
	
	m_eps = 1e-12;
}

//-----------------------------------------------------------------------------
//! data initialization and checking
void FEOgdenUnconstrained::Init()
{
	FEElasticMaterial::Init();
	for (int i=0; i<MAX_TERMS; ++i) 
		if (m_m[i] == 0) throw MaterialError("Invalid value for m%d", i+1);
}

//-----------------------------------------------------------------------------

void FEOgdenUnconstrained::EigenValues(mat3ds& A, double l[3], vec3d r[3], const double eps)
{
	A.eigen(l, r);
	
	// correct for numerical inaccuracy
	double d01 = fabs(l[0] - l[1]);
	double d12 = fabs(l[1] - l[2]);
	double d02 = fabs(l[0] - l[2]);
	
	if (d01 < eps) l[1] = l[0]; //= 0.5*(l[0]+l[1]);
	if (d02 < eps) l[2] = l[0]; //= 0.5*(l[0]+l[2]);
	if (d12 < eps) l[2] = l[1]; //= 0.5*(l[1]+l[2]);
	
}

//-----------------------------------------------------------------------------
//! Calculates the Cauchy stress
mat3ds FEOgdenUnconstrained::Stress(FEMaterialPoint &mp)
{
	// extract elastic material data
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// jacobian
	double J = pt.m_J;
	
	// get the left Cauchy-Green tensor
	mat3ds b = pt.LeftCauchyGreen();
	
	// get the eigenvalues and eigenvectors of b
	double lam2[3];	// these are the squares of the eigenvalues of V
	vec3d ev[3];
	EigenValues(b, lam2, ev, m_eps);
	
	// get the eigenvalues of V
	double lam[3];
	lam[0] = sqrt(lam2[0]);
	lam[1] = sqrt(lam2[1]);
	lam[2] = sqrt(lam2[2]);
	
	// stress
	mat3ds s;
	s.zero();
	double T;
	for (int i=0; i<3; ++i) {
		T = m_p*(J-1);
		for (int j=0; j<MAX_TERMS; ++j)
			T += m_c[j]/m_m[j]*(pow(lam[i], m_m[j]) - 1)/J;
		s += dyad(ev[i])*T;
	}

	return s;
}

//-----------------------------------------------------------------------------
//! Calculates the spatial tangent
tens4ds FEOgdenUnconstrained::Tangent(FEMaterialPoint& mp)
{
	int i,j,k;
	
	// extract elastic material data
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// jacobian
	double J = pt.m_J;
	
	// get the left Cauchy-Green tensor
	mat3ds b = pt.LeftCauchyGreen();
	
	// get the eigenvalues and eigenvectors of b
	double lam2[3];	// these are the squares of the eigenvalues of V
	vec3d ev[3];
	EigenValues(b, lam2, ev, m_eps);
	
	// get the eigenvalues of V
	double lam[3];
	mat3ds N[3];
	for (i=0; i<3; ++i) {
		lam[i] = sqrt(lam2[i]);
		N[i] = dyad(ev[i]);
	}
	
	// calculate the powers of eigenvalues
	double lamp[3][MAX_TERMS];
	for (j=0; j<MAX_TERMS; ++j)
	{
		lamp[0][j] = pow(lam[0], m_m[j]);
		lamp[1][j] = pow(lam[1], m_m[j]);
		lamp[2][j] = pow(lam[2], m_m[j]);
	}
	
	// principal stresses
	mat3ds s;
	s.zero();
	double T[3];
	for (i=0; i<3; ++i) {
		T[i] = m_p*(J-1);
		for (j=0; j<MAX_TERMS; ++j)
			T[i] += m_c[j]/m_m[j]*(lamp[i][j] - 1)/J;
		s += N[i]*T[i];
	}
	
	// coefficients appearing in elasticity tensor
	double D[3][3],E[3][3];
	for (j=0; j<3; ++j) {
		D[j][j] = m_p;
		for (k=0; k<MAX_TERMS; ++k)
			D[j][j] += m_c[k]/m_m[k]*((m_m[k]-2)*lamp[j][k]+2)/J;
		for (i=j+1; i<3; ++i) {
			D[i][j] = m_p*(2*J-1);
			if (lam2[j] != lam2[i])
				E[i][j] = 2*(lam2[j]*T[i] - lam2[i]*T[j])/(lam2[i]-lam2[j]);
			else {
				E[i][j] = 2*m_p*(1-J);
				for (k=0; k<MAX_TERMS; ++k)
					E[i][j] += m_c[k]/m_m[k]*((m_m[k]-2)*lamp[j][k]+2)/J;
			}
		}
	}
	
	// spatial elasticity tensor
	tens4ds c(0.0);
	mat3dd I(1.0);
	for (j=0; j<3; ++j) {
		c += dyad1s(N[j])*D[j][j];
		for (i=j+1; i<3; ++i) {
			c += dyad1s(N[i],N[j])*D[i][j];
			c += dyad4s(N[i],N[j])*E[i][j];
		}
	}
	
	return c;
}
