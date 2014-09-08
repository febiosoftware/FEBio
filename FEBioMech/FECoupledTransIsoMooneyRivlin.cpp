// FEMooneyRivlin.cpp: implementation of the FEMooneyRivlin class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FECoupledTransIsoMooneyRivlin.h"
#ifdef HAVE_GSL
#include "gsl/gsl_sf_expint.h"
#endif

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_PARAMETER_LIST(FECoupledTransIsoMooneyRivlin, FEElasticMaterial)
	ADD_PARAMETER(m_c1  , FE_PARAM_DOUBLE, "c1");
	ADD_PARAMETER(m_c2  , FE_PARAM_DOUBLE, "c2");
	ADD_PARAMETER(m_c3  , FE_PARAM_DOUBLE, "c3");
	ADD_PARAMETER(m_c4  , FE_PARAM_DOUBLE, "c4");
	ADD_PARAMETER(m_c5  , FE_PARAM_DOUBLE, "c5");
	ADD_PARAMETER(m_flam, FE_PARAM_DOUBLE, "lambda");
	ADD_PARAMETER(m_K   , FE_PARAM_DOUBLE, "k");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
void FECoupledTransIsoMooneyRivlin::Init()
{
	if (m_c1 <= 0) throw MaterialError("c1 must be positive");
	FEElasticMaterial::Init();
}

//-----------------------------------------------------------------------------
//! Calculate the Cauchy stress
mat3ds FECoupledTransIsoMooneyRivlin::Stress(FEMaterialPoint& mp)
{
	// get the material point data
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// calculate left Cauchy-Green tensor
	mat3ds B = pt.LeftCauchyGreen();

	// calculate square of B
	mat3ds B2 = B*B;

	// get the material fiber axis
	vec3d a0;
	a0.x = pt.m_Q[0][0];
	a0.y = pt.m_Q[1][0];
	a0.z = pt.m_Q[2][0];

	// get the spatial fiber axis
	vec3d a = pt.m_F*a0;
	double l = a.unit();

	// Invariants of B (= invariants of C)
	double I1 = B.tr();
	double J = pt.m_J;

	// some useful tensors
	mat3dd I(1.0);
	mat3ds A = dyad(a);

	// a. define the matrix stress
	//-----------------------------
	// W = c1*(I1 - 3) + c2*(I2 - 3)
	// Wi = dW/dIi
	double W1 = m_c1;
	double W2 = m_c2;

	mat3ds s = (B*(W1 + I1*W2) - B2*W2 - I*(W1 + 2.0*W2))*(2.0/J);

	// b. define fiber stress
	//-------------------------
	if (l > 1.0)
	{
		double Wl = 0.0;
		if (l < m_flam)
		{
			Wl = m_c3*(exp(m_c4*(l - 1.0)) - 1.0);
		}
		else
		{
			double c6 = m_c3*(exp(m_c4*(m_flam - 1.0)) - 1.0) - m_c5*m_flam;
			Wl = m_c5*l + c6;
		}
		s += A*(Wl/J);
	}

	// c. define dilational stress
	//------------------------------
	// U(J) = 1/2*k*(lnJ)^2
	double UJ = m_K*log(J)/J;
	s += I*UJ;

	return s;
}

//-----------------------------------------------------------------------------
//! Calculate the spatial elasticity tangent
tens4ds FECoupledTransIsoMooneyRivlin::Tangent(FEMaterialPoint& mp)
{
	// get material point data
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// calculate left Cauchy-Green tensor: B = F*Ft
	mat3ds B = pt.LeftCauchyGreen();

	// get the material fiber axis
	vec3d a0;
	a0.x = pt.m_Q[0][0];
	a0.y = pt.m_Q[1][0];
	a0.z = pt.m_Q[2][0];

	// get the spatial fiber axis
	vec3d a = pt.m_F*a0;
	double l = a.unit();

	// Invariants of B (= invariants of C)
	double J = pt.m_J;
	double I4 = l*l;

	// some useful tensors
	mat3dd I(1.0);
	mat3ds A = dyad(a);
	tens4ds IxI = dyad1s(I);
	tens4ds IoI = dyad4s(I);
	tens4ds BxB = dyad1s(B);
	tens4ds BoB = dyad4s(B);
	tens4ds AxA = dyad1s(A);

	// a. matrix tangent
	//----------------------------------
	double W1 = m_c1;
	double W2 = m_c2;

	tens4ds c = BxB*(4.0*W2/J) - BoB*(4.0*W2/J)  + IoI*(4.0*(W1+2.0*W2)/J);

	// b. fiber tangent
	// ---------------------------------
	if (l > 1.0)
	{
		double Fl = 0.0, Fll = 0.0;
		if (l < m_flam)
		{
			Fl = m_c3*(exp(m_c4*(l - 1.0)) - 1.0)/l;
			Fll = -m_c3*(exp(m_c4*(l-1.0)) - 1.0)/(l*l) + m_c3*m_c4*exp(m_c4*(l-1.0))/l;
		}
		else
		{
			double c6 = m_c3*(exp(m_c4*(m_flam - 1.0)) - 1.0) - m_c5*m_flam;
			Fl = m_c5 + c6 / l;
			Fll = -c6/(l*l);
		}

		double W44 = (Fll - Fl/l)/(4*l*l);

		c += AxA*(4.0*W44*I4*I4/J);
	}

	// c. dilational tangent
	// ---------------------------------
	double UJ = m_K*log(J)/J;
	double UJJ = m_K*(1 - log(J))/(J*J);
	c += IxI*(UJ + J*UJJ) - IoI*(2.0*UJ);

	return c;
}

#ifdef HAVE_GSL
//-----------------------------------------------------------------------------
//! calculate strain energy density at material point
double FECoupledTransIsoMooneyRivlin::StrainEnergyDensity(FEMaterialPoint& mp)
{
	// get the material point data
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
	// calculate left Cauchy-Green tensor
	mat3ds B = pt.LeftCauchyGreen();
    
	// calculate square of B
	mat3ds B2 = B*B;
    
	// get the material fiber axis
	vec3d a0;
	a0.x = pt.m_Q[0][0];
	a0.y = pt.m_Q[1][0];
	a0.z = pt.m_Q[2][0];
    
	// get the spatial fiber axis
	vec3d a = pt.m_F*a0;
	double l = a.unit();
    
	// Invariants of B (= invariants of C)
	double I1 = B.tr();
    double I2 = 0.5*(I1*I1 - B2.tr());
	double J = pt.m_J;
    double lnJ = log(J);
    
	// a. define the matrix sed
	//-----------------------------
	// W = c1*(I1 - 3) + c2*(I2 - 3)
    double sed = m_c1*(I1 - 3) + m_c2*(I2 - 3) - 2*(m_c1 + 2*m_c2)*lnJ;
    
	// b. define fiber strain energy density
	//-------------------------
	if (l > 1.0)
	{
		if (l < m_flam)
		{
            sed += m_c3*(exp(-m_c4)*(gsl_sf_expint_Ei(m_c4*l) - gsl_sf_expint_Ei(m_c4))-log(l));
		}
		else
		{
			double c6 = m_c3*(exp(m_c4*(m_flam - 1.0)) - 1.0) - m_c5*m_flam;
			sed += m_c5*(l-1) +c6*log(l);
		}
	}
    
	// c. define dilational stress
	//------------------------------
	// U(J) = 1/2*k*(lnJ)^2
	sed += m_K*lnJ*lnJ/2;
    
    return sed;
}
#endif
