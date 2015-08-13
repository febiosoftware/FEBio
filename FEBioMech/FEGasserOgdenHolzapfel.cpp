#include "stdafx.h"
#include "FEGasserOgdenHolzapfel.h"

// define the material parameters
BEGIN_PARAMETER_LIST(FEGasserOgdenHolzapfel, FEElasticMaterial)
	ADD_PARAMETER2(m_c, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "c");
	ADD_PARAMETER2(m_k1, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "k1");
	ADD_PARAMETER2(m_k2, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "k2");
	ADD_PARAMETER2(m_kappa, FE_PARAM_DOUBLE, FE_RANGE_CLOSED(0.0, 1.0/3.0), "kappa");
	ADD_PARAMETER(m_g, FE_PARAM_DOUBLE, "gamma");
    ADD_PARAMETER2(m_k, FE_PARAM_DOUBLE, FE_RANGE_GREATER(0.0), "k");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! Calculates the deviatoric stress
mat3ds FEGasserOgdenHolzapfel::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	vec3d n[3];			// local element basis directions
	vec3d a[2];			// structural direction in current configuration
	mat3ds h[2];		// structural tensor in current configuration
	double E[2];		// fiber strain

	// determinant of deformation gradient
	double J = pt.m_J;
	
	// Evaluate the deformation gradient
	mat3d F = pt.m_F;
	
	// calculate left Cauchy-Green tensor: b = F*Ft
	mat3ds b = pt.LeftCauchyGreen();
	
	// Copy the local element basis directions to n
	n[0].x = pt.m_Q[0][0]; n[0].y = pt.m_Q[1][0]; n[0].z = pt.m_Q[2][0];
	n[1].x = pt.m_Q[0][1]; n[1].y = pt.m_Q[1][1]; n[1].z = pt.m_Q[2][1];

	// Evaluate the structural direction in the current configuration
	double cg = cos(m_g); double sg = sin(m_g);
	a[0] = F*(n[0]*cg + n[1]*sg);
	a[1] = F*(n[0]*cg - n[1]*sg);
	
	// Evaluate the ground matrix stress
	mat3ds s = m_c/J*b + m_k*log(J)/J*mat3dd(1);
	
	// Evaluate the structural tensors in the current configuration
	// and the fiber strains and stress contributions
	h[0] = h[1] = m_kappa*b;
	if (a[0]*a[0] > 1) {
		h[0] += (1-3*m_kappa)*dyad(a[0]);
		E[0] = h[0].tr() - 1;
		s += 2./J*m_k1*E[0]*exp(m_k2*E[0]*E[0])*h[0];
	}
	if (a[1]*a[1] > 1) {
		h[1] += (1-3*m_kappa)*dyad(a[1]);
		E[1] = h[1].tr() - 1;
		s += 2./J*m_k1*E[1]*exp(m_k2*E[1]*E[1])*h[1];
	}
	
	return s;
}

//-----------------------------------------------------------------------------
//! Calculates the deviatoric tangent
tens4ds FEGasserOgdenHolzapfel::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	vec3d n[3];			// local element basis directions
	vec3d a[2];			// structural direction in current configuration
	mat3ds h[2];		// structural tensor in current configuration
	double E[2];		// fiber strain
	
	// determinant of deformation gradient
	double J = pt.m_J;
	
	// Evaluate the deformation gradient
	mat3d F = pt.m_F;
	
	// calculate left Cauchy-Green tensor: b = F*Ft
	mat3ds b = pt.LeftCauchyGreen();
	
	// Copy the local element basis directions to n
	n[0].x = pt.m_Q[0][0]; n[0].y = pt.m_Q[1][0]; n[0].z = pt.m_Q[2][0];
	n[1].x = pt.m_Q[0][1]; n[1].y = pt.m_Q[1][1]; n[1].z = pt.m_Q[2][1];
	
	// Evaluate the structural direction in the current configuration
	double cg = cos(m_g); double sg = sin(m_g);
	a[0] = F*(n[0]*cg + n[1]*sg);
	a[1] = F*(n[0]*cg - n[1]*sg);
	
	// Evaluate the ground matrix stress
	mat3ds s = m_c/J*b;
	
	// Evaluate the structural tensors in the current configuration
	// and the fiber strains and stress contributions
	h[0] = h[1] = m_kappa*b;
	if (a[0]*a[0] > 1) {
		h[0] += (1-3*m_kappa)*dyad(a[0]);
		E[0] = h[0].tr() - 1;
		s += 2./J*m_k1*E[0]*exp(m_k2*E[0]*E[0])*h[0];
	}
	if (a[1]*a[1] > 1) {
		h[1] += (1-3*m_kappa)*dyad(a[1]);
		E[1] = h[1].tr() - 1;
		s += 2./J*m_k1*E[1]*exp(m_k2*E[1]*E[1])*h[1];
	}
	
	// Evaluate the elasticity tensor
	mat3dd I(1);
	tens4ds IxI = dyad1s(I);
    tens4ds I4  = dyad4s(I);
    tens4ds C = (IxI - I4*(2*log(J)))*(m_k/J);
	if (a[0]*a[0] > 1)
		C += 4./J*m_k1*(1+2*m_k2*E[0]*E[0])*exp(m_k2*E[0]*E[0])*
		(dyad1s(h[0]) - h[0].tr()/3.*(dyad1s(h[0],I) - h[0].tr()/3.*IxI));
	if (a[1]*a[1] > 1)
		C += 4./J*m_k1*(1+2*m_k2*E[1]*E[1])*exp(m_k2*E[1]*E[1])*
		(dyad1s(h[1]) - h[1].tr()/3.*(dyad1s(h[1],I) - h[1].tr()/3.*IxI));
	
	return C;
}

//-----------------------------------------------------------------------------
//! Calculates the strain energy density
double FEGasserOgdenHolzapfel::StrainEnergyDensity(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	vec3d n[3];			// local element basis directions
	vec3d a[2];			// structural direction in current configuration
	mat3ds h[2];		// structural tensor in current configuration
	double E[2];		// fiber strain
    
	// Evaluate the deformation gradient
	mat3d F = pt.m_F;
	
    // determinant of deformation gradient
    double J = pt.m_J;
    double lnJ = log(J);
    
	// calculate left Cauchy-Green tensor: b = F*Ft
	mat3ds b = pt.LeftCauchyGreen();
    double I1 = b.tr();
	
	// Copy the local element basis directions to n
	n[0].x = pt.m_Q[0][0]; n[0].y = pt.m_Q[1][0]; n[0].z = pt.m_Q[2][0];
	n[1].x = pt.m_Q[0][1]; n[1].y = pt.m_Q[1][1]; n[1].z = pt.m_Q[2][1];
    
	// Evaluate the structural direction in the current configuration
	double cg = cos(m_g); double sg = sin(m_g);
	a[0] = F*(n[0]*cg + n[1]*sg);
	a[1] = F*(n[0]*cg - n[1]*sg);
	
	// Evaluate the ground matrix strain energy density
    double sed = 0.5*(m_c*(I1 - 3) + m_k*lnJ*lnJ);
	
	// Evaluate the structural tensors in the current configuration
	// and the fiber strains and stress contributions
	h[0] = h[1] = m_kappa*b;
	if (a[0]*a[0] > 1) {
		h[0] += (1-3*m_kappa)*dyad(a[0]);
		E[0] = h[0].tr() - 1;
        sed += 0.5*m_k1/m_k2*(exp(m_k2*E[0]*E[0])-1);
	}
	if (a[1]*a[1] > 1) {
		h[1] += (1-3*m_kappa)*dyad(a[1]);
		E[1] = h[1].tr() - 1;
		sed += 0.5*m_k1/m_k2*(exp(m_k2*E[1]*E[1])-1);
	}
	
	return sed;
}
