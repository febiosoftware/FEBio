#include "stdafx.h"
#include "FEEFDVerondaWestmann.h"

// The following file contains the integration points and weights
// for the integration over a unit sphere in spherical coordinates
#include "geodesic.h"

// we store the cos and sin of the angles here
int FEEFDVerondaWestmann::m_nres = 0;
double FEEFDVerondaWestmann::m_cth[NSTH];
double FEEFDVerondaWestmann::m_sth[NSTH];
double FEEFDVerondaWestmann::m_cph[NSTH];
double FEEFDVerondaWestmann::m_sph[NSTH];
double FEEFDVerondaWestmann::m_w[NSTH];

// register the material with the framework
REGISTER_MATERIAL(FEEFDVerondaWestmann, "EFD Veronda-Westmann");

// define the material parameters
BEGIN_PARAMETER_LIST(FEEFDVerondaWestmann, FEUncoupledMaterial)
	ADD_PARAMETER(m_c1, FE_PARAM_DOUBLE, "c1");
	ADD_PARAMETER(m_c2, FE_PARAM_DOUBLE, "c2");
	ADD_PARAMETERV(m_beta, FE_PARAM_DOUBLEV, 3, "beta");
	ADD_PARAMETERV(m_ksi , FE_PARAM_DOUBLEV, 3, "ksi" );
	ADD_PARAMETERV(m_a, FE_PARAM_DOUBLEV, 3, "a");
	ADD_PARAMETER(m_ac, FE_PARAM_DOUBLE, "active_contraction");
END_PARAMETER_LIST();

#ifndef SQR
	#define SQR(x) ((x)*(x))
#endif

//////////////////////////////////////////////////////////////////////
// FEEFDVerondaWestmann
//////////////////////////////////////////////////////////////////////

FEEFDVerondaWestmann::FEEFDVerondaWestmann()
{
	static bool bfirst = true;

	if (bfirst)
	{
		// select the integration rule
		const int nint = (m_nres == 0? NSTL  : NSTH  );
		const double* phi = (m_nres == 0? PHIL  : PHIH  );
		const double* the = (m_nres == 0? THETAL: THETAH);
		const double* w   = (m_nres == 0? AREAL : AREAH );

		for (int n=0; n<nint; ++n)
		{
			m_cth[n] = cos(the[n]);
			m_sth[n] = sin(the[n]);
			m_cph[n] = cos(phi[n]);
			m_sph[n] = sin(phi[n]);
			m_w[n] = w[n];
		}

		bfirst = false;
	}

	m_a[0] = m_a[1] = m_a[2] = 0;
	m_ac = 0;
}

//-----------------------------------------------------------------------------
void FEEFDVerondaWestmann::Init()
{
	FEUncoupledMaterial::Init();
	
	if (m_ksi[0] <= 0) throw MaterialError("ksi1 must be positive.");
	if (m_ksi[1] <= 0) throw MaterialError("ksi2 must be positive.");
	if (m_ksi[2] <= 0) throw MaterialError("ksi3 must be positive.");
	if (m_beta[0] <= 2) throw MaterialError("beta1 must be bigger than 2.");
	if (m_beta[1] <= 2) throw MaterialError("beta1 must be bigger than 2.");
	if (m_beta[2] <= 2) throw MaterialError("beta1 must be bigger than 2.");
}

//-----------------------------------------------------------------------------
mat3ds FEEFDVerondaWestmann::DevStress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// deformation gradient and its determinant
	mat3d &F = pt.F;
	double J = pt.J;

	// calculate deviatoric left Cauchy-Green tensor
	mat3ds B = pt.DevLeftCauchyGreen();

	// calculate square of B
	mat3ds B2 = B*B;

	// Invariants of B (= invariants of C)
	// Note that these are the invariants of Btilde, not of B!
	double I1 = B.tr();
	double I2 = 0.5*(I1*I1 - B2.tr());

	// --- TODO: put strain energy derivatives here ---
	//
	// W = C1*(exp(C2*(I1-3)-1)-0.5*C1*C2*(I2 - 3)
	//
	// Wi = dW/dIi
	double W1 = m_c1*m_c2*exp(m_c2*(I1-3));
	double W2 = -0.5*m_c1*m_c2;
	// ---

	// calculate T = F*dW/dC*Ft
	mat3ds T = B*(W1 + W2*I1) - B2*W2;

	// Add fiber contribution
	T += FiberStress(mp);

	return T.dev()*(2.0/J);
}

//-----------------------------------------------------------------------------
//! calculate the fiber contribution to the deviatoric stress
mat3ds FEEFDVerondaWestmann::FiberStress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// deformation gradient
	mat3d &F = pt.F;
	double J = pt.J;
	double Jm13 = pow(J, -1.0/3.0);

	// calculate deviatoric right Cauchy-Green tensor
	mat3ds C = pt.DevRightCauchyGreen();

	// get the element's local coordinate system
	mat3d& Q = pt.Q;

	// loop over all integration points
	double ksi, beta;
	vec3d nr, n0, nt;
	double In, Wl;
	const double eps = 0;
	const int nint = (m_nres == 0? NSTL  : NSTH  );
	mat3ds Tf; Tf.zero();
	for (int n=0; n<nint; ++n)
	{
		// set the local fiber direction
		nr.x = m_cth[n]*m_sph[n];
		nr.y = m_sth[n]*m_sph[n];
		nr.z = m_cph[n];

		// get the global material fiber direction
		n0 = Q*nr;

		// get the global spatial fiber direction
		nt = (F*n0)*Jm13;

		// Calculate In = nr*C*nr
		In = nt*nt;

		// calculate the outer product of nt
		mat3ds N = dyad(nt);

		// only take fibers in tension into consideration
		if (In > 1. + eps)
		{
			// calculate material coefficients
			ksi  = 1.0 / sqrt(SQR(nr.x / m_ksi [0]) + SQR(nr.y / m_ksi [1]) + SQR(nr.z / m_ksi [2]));
			beta = 1.0 / sqrt(SQR(nr.x / m_beta[0]) + SQR(nr.y / m_beta[1]) + SQR(nr.z / m_beta[2]));

			// calculate strain energy derivative
			Wl = beta*ksi*pow(In - 1.0, beta-1.0);

			// calculate the stress
			Tf += N*(Wl*m_w[n]);
		}

		// add active contraction stuff
		if (m_ac > 0)
		{
			// The .5 is to compensate for the 2 multiplier later.
			double at = 0.5*m_w[n]*m_ac /sqrt(SQR(nr.x/m_a[0]) + SQR(nr.y / m_a[1]) + SQR(nr.z / m_a[2]));
			Tf += N*at;
		}
	}

	return Tf;
}

//-----------------------------------------------------------------------------
tens4ds FEEFDVerondaWestmann::DevTangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// deformation gradient
	mat3d &F = pt.F;
	double J = pt.J;
	double Ji = 1.0/J;

	// calculate deviatoric left Cauchy-Green tensor: B = F*Ft
	mat3ds B = pt.DevLeftCauchyGreen();

	// calculate square of B
	mat3ds B2 = B*B;

	// Invariants of B (= invariants of C)
	double I1 = B.tr();
	double I2 = 0.5*(I1*I1 - B2.tr());

	// --- TODO: put strain energy derivatives here ---
	// Wi = dW/dIi
	double W1, W2, W11;
	W1 = m_c1*m_c2*exp(m_c2*(I1-3));
	W2 = -0.5*m_c1*m_c2;
	W11 = m_c2*W1;
	// ---

	// calculate dWdC:C
	double WC = W1*I1 + 2*W2*I2;

	// calculate C:d2WdCdC:C
	double CWWC = W11*I1*I1+2*I2*W2;

	// deviatoric cauchy-stress, trs = trace[s]/3
	mat3ds devs = pt.s.dev();

	mat3ds I(1,1,1,0,0,0);	// Identity

	tens4ds IxI = dyad1s(I);
	tens4ds I4  = dyad4s(I);
	tens4ds BxB = dyad1s(B);
	tens4ds B4  = dyad4s(B);

	// d2W/dCdC:C
	mat3ds WCCxC = B*(I1*(W11 + W2)) - B2*W2;

	tens4ds cw = BxB*((W11 + W2)*4.0*Ji) - B4*(W2*4.0*Ji) - dyad1s(WCCxC, I)*(4.0/3.0*Ji) + IxI*(4.0/9.0*Ji*CWWC);

	tens4ds c = dyad1s(devs, I)*(-2.0/3.0) + (I4 - IxI/3.0)*(4.0/3.0*Ji*WC) + cw;

	// add fiber contribution
	c += FiberTangent(mp);
	
	return c;
}


//-----------------------------------------------------------------------------
//! Calculate fiber tangent
tens4ds FEEFDVerondaWestmann::FiberTangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// deformation gradient
	mat3d &F = pt.F;
	double J = pt.J;
	double Jm13 = pow(J, -1.0/3.0);
	double Ji = 1.0/J;

	// get the element's local coordinate system
	mat3d& Q = pt.Q;

	mat3dd I(1);	// Identity

	tens4ds IxI = dyad1s(I);
	tens4ds I4  = dyad4s(I);

	// loop over all integration points
	double ksi, beta;
	vec3d nr, n0, nt;
	double In, Wl, Wll;
	const double eps = 0;
	tens4ds cf, cfw; cf.zero();
	mat3ds N2;
	tens4ds N4;
	tens4ds I4mIxId3 = I4 - IxI/3.0;
	const int nint = (m_nres == 0? NSTL  : NSTH  );
	for (int n=0; n<nint; ++n)
	{
		// set the local fiber direction
		nr.x = m_cth[n]*m_sph[n];
		nr.y = m_sth[n]*m_sph[n];
		nr.z = m_cph[n];

		// get the global material fiber direction
		n0 = Q*nr;

		// get the global spatial fiber direction
		nt = (F*n0)*Jm13;

		// Calculate In = nr*C*nr
		In = nt*nt;

		// only take fibers in tension into consideration
		if (In > 1. + eps)
		{
			// calculate material coefficients
			ksi  = 1.0 / sqrt(SQR(nr.x / m_ksi [0]) + SQR(nr.y / m_ksi [1]) + SQR(nr.z / m_ksi [2]));
			beta = 1.0 / sqrt(SQR(nr.x / m_beta[0]) + SQR(nr.y / m_beta[1]) + SQR(nr.z / m_beta[2]));

			// calculate strain energy derivative
			Wl  = beta*ksi*pow(In - 1.0, beta - 1.0);
			Wll = beta*(beta-1.0)*ksi*pow(In - 1.0, beta-2.0);

			// calculate dWdC:C
			// TODO: is this correct? the nt vector is not normalized!!!
			double WC = Wl*In;

			// calculate C:d2WdCdC:C
			double CWWC = Wll*In*In;

			N2 = dyad(nt);
			N4 = dyad1s(N2);

			mat3ds WCCxC = N2*(Wll*In);

			cfw = N4*(4.0*Wll) - dyad1s(WCCxC, I)*(4.0/3.0) + IxI*(4.0/9.0*CWWC);

			cf += (I4mIxId3)*(4.0/3.0*Ji*WC*m_w[n]) + cfw*(Ji*m_w[n]);
		}
	}

	return cf;
}
