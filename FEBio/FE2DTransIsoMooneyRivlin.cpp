#include "stdafx.h"
#include "FE2DTransIsoMooneyRivlin.h"

// register the material with the framework
REGISTER_MATERIAL(FE2DTransIsoMooneyRivlin, "2D trans iso Mooney-Rivlin");

// define the material parameters
BEGIN_PARAMETER_LIST(FE2DTransIsoMooneyRivlin, FETransverselyIsotropic)
	ADD_PARAMETER(m_c1, FE_PARAM_DOUBLE, "c1");
	ADD_PARAMETER(m_c2, FE_PARAM_DOUBLE, "c2");
	ADD_PARAMETERV(m_a, FE_PARAM_DOUBLEV, 2, "a");
	ADD_PARAMETER(m_ac, FE_PARAM_DOUBLE, "active_contraction");
	ADD_PARAMETERV(m_w, FE_PARAM_DOUBLEV, 2, "w");
END_PARAMETER_LIST();

double FE2DTransIsoMooneyRivlin::m_cth[FE2DTransIsoMooneyRivlin::NSTEPS];
double FE2DTransIsoMooneyRivlin::m_sth[FE2DTransIsoMooneyRivlin::NSTEPS];

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

//////////////////////////////////////////////////////////////////////
// FE2DTransIsoMooneyRivlin
//////////////////////////////////////////////////////////////////////

FE2DTransIsoMooneyRivlin::FE2DTransIsoMooneyRivlin()
{
	static bool bfirst = true;

	if (bfirst)
	{
		double ph;
		const double PI = 4.0*atan(1.0);
		for (int n=0; n<NSTEPS; ++n)
		{
			ph = 2.0*PI*n / (double) NSTEPS;
			m_cth[n] = cos(ph);
			m_sth[n] = sin(ph);
		}
		bfirst = false;
	}

	m_c1 = 0;
	m_c2 = 0;

	m_a[0] = m_a[1] = 1;
	m_ac = 0;

	m_w[0] = m_w[1] = 1;
}

//-----------------------------------------------------------------------------
//! Calculates the stress for this material.
//! \param pt material point at which to evaluate the stress
mat3ds FE2DTransIsoMooneyRivlin::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	const double third = 1.0/3.0;

	// deformation gradient
	mat3d &F = pt.F;
	double J = pt.J;
	double Jm13 = pow(J, -1.0/3.0);

	// calculate deviatoric left Cauchy-Green tensor
	mat3ds B = pt.DevLeftCauchyGreen();

	// calculate square of B
	mat3ds B2 = B*B;

	// Invariants of B (= invariants of C)
	// Note that these are the invariants of deviatoric B, not of B!
	double I1 = B.tr();
	double I2 = 0.5*(I1*I1 - B2.tr());

	// strain energy derivatives
	double W1 = m_c1;
	double W2 = m_c2;

	// --- M A T R I X   C O N T R I B U T I O N ---

	// calculate T = F*dW/dC*Ft
	mat3ds T = B*(W1 + W2*I1) - B2*W2;

	// --- F I B E R   C O N T R I B U T I O N ---

	// Next, we calculate the fiber contribution. For this material
	// the fibers lie randomly in a plane that is perpendicular to the transverse
	// axis. We therefor need to integrate over this plane.
	const double PI = 4.0*atan(1.0);
	double w, wtot = 0;
	vec3d a0, a, v;
	double lam, lamd, I4, W4;
	mat3ds Tf; Tf.zero();
	mat3ds N;
	for (int n=0; n<NSTEPS; ++n)
	{
		// calculate the local material fiber vector
		v.y = m_cth[n];
		v.z = m_sth[n];
		v.x = 0;

		// calculate the global material fiber vector
		a0 = pt.Q*v;

		// calculate the global spatial fiber vector
		a.x = F[0][0]*a0.x + F[0][1]*a0.y + F[0][2]*a0.z;
		a.y = F[1][0]*a0.x + F[1][1]*a0.y + F[1][2]*a0.z;
		a.z = F[2][0]*a0.x + F[2][1]*a0.y + F[2][2]*a0.z;

		// normalize material axis and store fiber stretch
		lam = a.unit();
		lamd = lam*Jm13; // i.e. lambda tilde
	
		// fourth invariant of C
		I4 = lamd*lamd;

		if (lamd > 1)
		{
			double lamdi = 1.0/lamd;
			double Wl;
			if (lamd < lam1)
			{
				Wl = lamdi*c3*(exp(c4*(lamd - 1)) - 1);
			}
			else
			{
				double c6 = c3*(exp(c4*(lam1-1))-1) - c5*lam1;
				Wl = lamdi*(c5*lamd + c6);
			}
			W4  = 0.5*lamdi*Wl;
		}
		else 
		{
			W4 = 0;
		}

		w = 1.0/sqrt((v.y/m_w[0])*(v.y/m_w[0]) + (v.z/m_w[1])*(v.z/m_w[1]));
		wtot += w;

		// calculate the stress
		N = dyad(a);
		Tf += N*(W4*I4*w);

		// add active contraction stuff
		if (m_ac > 0)
		{
			// The .5 is to compensate for the 2 multiplier later.
			double at = 0.5*w*m_ac /sqrt(SQR(a.x/m_a[0]) + SQR(a.y / m_a[1]) + SQR(a.z / m_a[2]));
			Tf += N*at;
		}
	}

	// add fiber to total
	T += Tf/wtot;

	// average element pressure
	double p = pt.avgp;

	// set the stress
	mat3dd I(1);
	mat3ds s = I*p + T.dev()*(2.0/J);

	return s;
}

//-----------------------------------------------------------------------------
//! Calculates the elasticity tensor for this material.
//! \param D elasticity tensor
//! \param pt material point at which to evaulate the elasticity tensor
tens4ds FE2DTransIsoMooneyRivlin::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// deformation gradient
	mat3d &F = pt.F;
	double J = pt.J;
	double Jm13 = pow(J, -1.0/3.0);
	double Jm23 = Jm13*Jm13;
	double Ji = 1.0/J;

	// deviatoric cauchy-stress, trs = trace[s]/3
	mat3ds& devs = pt.s.dev();

	// mean pressure
	double p = pt.avgp;

	// deviatoric right Cauchy-Green tensor: C = Ft*F
	mat3ds C = pt.DevRightCauchyGreen();

	// square of C
	mat3ds C2 = C*C;

	// Invariants of C
	double I1 = C.tr();
	double I2 = 0.5*(I1*I1 - C2.tr());

	// calculate left Cauchy-Green tensor: B = F*Ft
	mat3ds B = pt.DevLeftCauchyGreen();

	// calculate square of B
	// (we commented out the components we don't need)
	mat3ds B2 = B*B;

	// strain energy derivatives
	double W1 = m_c1;
	double W2 = m_c2;

	// --- M A T R I X   C O N T R I B U T I O N ---

	// calculate dWdC:C
	double WC = W1*I1 + 2*W2*I2;

	// calculate C:d2WdCdC:C
	double CWWC = 2*I2*W2;

	mat3dd I(1);	// Identity

	tens4ds IxI = dyad1s(I);
	tens4ds I4  = dyad4s(I);
	tens4ds BxB = dyad1s(B);
	tens4ds B4  = dyad4s(B);

	// d2W/dCdC:C
	mat3ds WCCxC = B*(W2*I1) - B2*W2;

	tens4ds cw = (BxB - B4)*(W2*4.0*Ji) - dyad1s(WCCxC, I)*(4.0/3.0*Ji) + IxI*(4.0/9.0*Ji*CWWC);

	tens4ds c = (IxI - I4*2)*p - dyad1s(devs, I)*(2.0/3.0) + (I4 - IxI/3.0)*(4.0/3.0*Ji*WC) + cw;

	// --- F I B E R   C O N T R I B U T I O N ---

	// Next, we add the fiber contribution. Since the fibers lie
	// randomly perpendicular to the transverse axis, we need
	// to integrate over that plane
	const double PI = 4.0*atan(1.0);
	double lam, lamd;
	double In, Wl, Wll;
	vec3d a0, a, v;
	double w, wtot = 0;
	tens4ds cf, cfw; cf.zero();
	mat3ds N2;
	tens4ds N4;
	tens4ds I4mIxId3 = I4 - IxI/3.0;
	for (int n=0; n<NSTEPS; ++n)
	{
		// calculate the local material fiber vector
		v.y = m_cth[n];
		v.z = m_sth[n];
		v.x = 0;

		// calculate the global material fiber vector
		a0 = pt.Q*v;

		// calculate the global spatial fiber vector
		a.x = F[0][0]*a0.x + F[0][1]*a0.y + F[0][2]*a0.z;
		a.y = F[1][0]*a0.x + F[1][1]*a0.y + F[1][2]*a0.z;
		a.z = F[2][0]*a0.x + F[2][1]*a0.y + F[2][2]*a0.z;

		// normalize material axis and store fiber stretch
		lam = a.unit();
		lamd = lam*Jm13; // i.e. lambda tilde
	
		// fourth invariant of C
		In = lamd*lamd;

		// Wi = dW/dIi
		if (lamd >= 1)
		{
			double lamdi = 1.0/lamd;
			double W4, W44;
			if (lamd < lam1)
			{
				W4  = lamdi*c3*(exp(c4*(lamd - 1)) - 1);
				W44 = c3*lamdi*(c4*exp(c4*(lamd - 1)) - lamdi*(exp(c4*(lamd-1))-1));
			}
			else
			{
				double c6 = c3*(exp(c4*(lam1-1))-1) - c5*lam1;
				W4  = lamdi*(c5*lamd + c6);
				W44 = -c6*lamdi*lamdi;
			}
			Wl  = 0.5*lamdi*W4;
			Wll = 0.25*lamdi*lamdi*(W44 - lamdi*W4);
		}
		else 
		{
			Wl = 0;
			Wll = 0;
		}

		w = 1.0/sqrt((v.y/m_w[0])*(v.y/m_w[0]) + (v.z/m_w[1])*(v.z/m_w[1]));
		wtot += w;

		// calculate dWdC:C
		double WC = Wl*In;

		// calculate C:d2WdCdC:C
		double CWWC = Wll*In*In;

		N2 = dyad(a);
		N4 = dyad1s(N2);

		WCCxC = N2*(Wll*In);

		cfw = N4*(4.0*Wll) - dyad1s(WCCxC, I)*(4.0/3.0) + IxI*(4.0/9.0*CWWC);

		cf += (I4mIxId3)*(4.0/3.0*Ji*WC*w) + cfw*(Ji*w);
	}

	// add fiber to total
	c += cf/wtot;

	return c;
}
