#include "stdafx.h"
#include "FERandomFiberMaterial.h"

// TODO: there is a possible bug in the fiber stress and or fiber stiffness

// The following file contains the integration points and weights
// for the integration over a unit sphere in spherical coordinates
#include "geodesic.h"

// we store the cos and sin of the angles here
int FERandomFiberMooneyRivlin::m_nres = 0;
double FERandomFiberMooneyRivlin::m_cth[NSTH];
double FERandomFiberMooneyRivlin::m_sth[NSTH];
double FERandomFiberMooneyRivlin::m_cph[NSTH];
double FERandomFiberMooneyRivlin::m_sph[NSTH];
double FERandomFiberMooneyRivlin::m_w[NSTH];

// register the material with the framework
REGISTER_MATERIAL(FERandomFiberMooneyRivlin, "random fiber Mooney-Rivlin");

// define the material parameters
BEGIN_PARAMETER_LIST(FERandomFiberMooneyRivlin, FEIncompressibleMaterial)
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
// FERandomFiberMooneyRivlin
//////////////////////////////////////////////////////////////////////

FERandomFiberMooneyRivlin::FERandomFiberMooneyRivlin()
{
	static bool bfirst = true;

	if (bfirst)
	{
		// select the integration rule
		m_nint = (m_nres == 0? NSTL  : NSTH  );
		const double* phi = (m_nres == 0? PHIL  : PHIH  );
		const double* the = (m_nres == 0? THETAL: THETAH);
		const double* w   = (m_nres == 0? AREAL : AREAH );

		for (int n=0; n<m_nint; ++n)
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

void FERandomFiberMooneyRivlin::Init()
{
	FEElasticMaterial::Init();

	if (m_ksi[0] <= 0) throw MaterialError("ksi1 must be positive.");
	if (m_ksi[1] <= 0) throw MaterialError("ksi2 must be positive.");
	if (m_ksi[2] <= 0) throw MaterialError("ksi3 must be positive.");
	if (m_beta[0] <= 2) throw MaterialError("beta1 must be bigger than 2.");
	if (m_beta[1] <= 2) throw MaterialError("beta1 must be bigger than 2.");
	if (m_beta[2] <= 2) throw MaterialError("beta1 must be bigger than 2.");
}

mat3ds FERandomFiberMooneyRivlin::Stress(FEMaterialPoint& mp)
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

	// calculate deviatoric right Cauchy-Green tensor
	mat3ds C = pt.DevRightCauchyGreen();

	// get the element's local coordinate system
	mat3d& Q = pt.Q;

	// loop over all integration points
	double ksi, beta;
	vec3d nr, n0, nt;
	double In, Wl;
	const double eps = 0;
	for (int n=0; n<m_nint; ++n)
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
			T += N*(Wl*m_w[n]);
		}

		// add active contraction stuff
		if (m_ac > 0)
		{
			// The .5 is to compensate for the 2 multiplier later.
			double at = 0.5*m_w[n]*m_ac *sqrt(SQR(m_a[0]*nr.x) + SQR(m_a[1]*nr.y) + SQR(m_a[2]*nr.z));
			T += N*at;
		}
	}

	// --- END FIBER CONTRIBUTION --

	// average element pressure
	double p = pt.avgp;

	// set the stress
	mat3dd I(1);
	mat3ds s = I*p + T.dev()*(2.0/J);

	return s;
}


void FERandomFiberMooneyRivlin::Tangent(double D[6][6], FEMaterialPoint& mp)
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

	// D[0][0] = c(0,0,0,0)
	D[0][0] = -p - (4.0/3.0)*devs(0,0) + (8.0/9.0)*Ji*WC;
	D[0][0] -= (8.0/3.0)*Ji*(W2*I1*B(0,0) - W2*B2(0,0));
	D[0][0] += (4.0/9.0)*Ji*CWWC;

	// D[1][1] = c(1,1,1,1)
	D[1][1] = -p - (4.0/3.0)*devs(1,1) + (8.0/9.0)*Ji*WC;
	D[1][1] -= (8.0/3.0)*Ji*(W2*I1*B(1,1) - W2*B2(1,1));
	D[1][1] += (4.0/9.0)*Ji*CWWC;

	// D[2][2] = c(2,2,2,2)
	D[2][2] = -p - (4.0/3.0)*devs(2,2) + (8.0/9.0)*Ji*WC;
	D[2][2] -= (8.0/3.0)*Ji*(W2*I1*B(2,2) - W2*B2(2,2));
	D[2][2] += (4.0/9.0)*Ji*CWWC;



	// D[0][1] = D[1][0] = c(0,0,1,1)
	D[0][1] = p - (2.0/3.0)*(devs(0,0) + devs(1,1)) - (4.0/9.0)*Ji*WC;
	D[0][1] += 4.0*Ji*W2*(B(0,0)*B(1,1) - B(0,1)*B(0,1));
	D[0][1] += (4.0/9.0)*Ji*CWWC;
	D[0][1] -= (4.0/3.0)*Ji*(W2*(I1*B(0,0) - B2(0,0)));
	D[0][1] -= (4.0/3.0)*Ji*(W2*(I1*B(1,1) - B2(1,1)));

	// D[1][2] = D[2][1] = c(1,1,2,2)
	D[1][2] = p - (2.0/3.0)*(devs(1,1) + devs(2,2)) - (4.0/9.0)*Ji*WC;
	D[1][2] += 4.0*Ji*W2*(B(1,1)*B(2,2) - B(1,2)*B(1,2));
	D[1][2] += (4.0/9.0)*Ji*CWWC;
	D[1][2] -= (4.0/3.0)*Ji*(W2*(I1*B(1,1) - B2(1,1)));
	D[1][2] -= (4.0/3.0)*Ji*(W2*(I1*B(2,2) - B2(2,2)));

	// D[0][2] = D[2][0] = c(0,0,2,2)
	D[0][2] = p - (2.0/3.0)*(devs(0,0) + devs(2,2)) - (4.0/9.0)*Ji*WC;
	D[0][2] += 4.0*Ji*W2*(B(0,0)*B(2,2) - B(0,2)*B(0,2));
	D[0][2] += (4.0/9.0)*Ji*CWWC;
	D[0][2] -= (4.0/3.0)*Ji*(W2*(I1*B(0,0) - B2(0,0)));
	D[0][2] -= (4.0/3.0)*Ji*(W2*(I1*B(2,2) - B2(2,2)));



	// D[3][3] = 0.5*(c(0,1,0,1) + c(0,1,1,0))
	D[3][3] = -p + (2.0/3.0)*Ji*WC;
	D[3][3] += 2.0*Ji*W2*(B(0,1)*B(0,1) - B(0,0)*B(1,1));

	// D[4][4] = 0.5*(c(1,2,1,2) + c(1,2,2,1))
	D[4][4] = -p + (2.0/3.0)*Ji*WC;
	D[4][4] += 2.0*Ji*W2*(B(1,2)*B(1,2) - B(1,1)*B(2,2));

	// D[5][5] = 0.5*(c(0,2,0,2) + c(0,2,2,0))
	D[5][5] = -p + (2.0/3.0)*Ji*WC;
	D[5][5] += 2.0*Ji*W2*(B(0,2)*B(0,2) - B(0,0)*B(2,2));



	// D[0][3] = 0.5*(c(0,0,0,1) + c(0,0,1,0))
	D[0][3] =  -(2.0/3.0)*devs(0,1);
	D[0][3] -= (4.0/3.0)*Ji*(W2*(I1*B(0,1) - B2(0,1)));

	// D[0][4] = 0.5*(c(0,0,1,2) + c(0,0,2,1))
	D[0][4] =  -(2.0/3.0)*devs(1,2);
	D[0][4] += 4.0*Ji*W2*(B(0,0)*B(1,2) - B(0,1)*B(0,2));
	D[0][4] -= (4.0/3.0)*Ji*(W2*(I1*B(1,2) - B2(1,2)));

	// D[0][5] = 0.5*(c(0,0,0,2) + c(0,0,2,0))
	D[0][5] =  -(2.0/3.0)*devs(0,2);
	D[0][5] -= (4.0/3.0)*Ji*(W2*(I1*B(0,2) - B2(0,2)));

	// D[1][3] = 0.5*(c(1,1,0,1) + c(1,1,1,0))
	D[1][3] =  -(2.0/3.0)*devs(0,1);
	D[1][3] -= (4.0/3.0)*Ji*(W2*(I1*B(0,1) - B2(0,1)));

	// D[1][4] = 0.5*(c(1,1,1,2) + c(1,1,2,1))
	D[1][4] =  -(2.0/3.0)*devs(1,2);
	D[1][4] -= (4.0/3.0)*Ji*(W2*(I1*B(1,2) - B2(1,2)));

	// D[1][5] = 0.5*(c(1,1,0,2) + c(1,1,2,0))
	D[1][5] =  -(2.0/3.0)*devs(0,2);
	D[1][5] += 4.0*Ji*W2*(B(1,1)*B(0,2) - B(0,1)*B(1,2));
	D[1][5] -= (4.0/3.0)*Ji*(W2*(I1*B(0,2) - B2(0,2)));

	// D[2][3] = 0.5*(c(2,2,0,1) + c(2,2,1,0))
	D[2][3] =  -(2.0/3.0)*devs(0,1);
	D[2][3] += 4.0*Ji*W2*(B(2,2)*B(0,1) - B(0,2)*B(1,2));
	D[2][3] -= (4.0/3.0)*Ji*(W2*(I1*B(0,1) - B2(0,1)));

	// D[2][4] = 0.5*(c(2,2,1,2) + c(2,2,2,1))
	D[2][4] =  -(2.0/3.0)*devs(1,2);
	D[2][4] -= (4.0/3.0)*Ji*(W2*(I1*B(1,2) - B2(1,2)));

	// D[2][5] = 0.5*(c(2,2,0,2) + c(2,2,2,0))
	D[2][5] =  -(2.0/3.0)*devs(0,2);
	D[2][5] -= (4.0/3.0)*Ji*(W2*(I1*B(0,2) - B2(0,2)));



	// D[3][4] = 0.5*(c(0,1,1,2) + c(0,1,2,1))
	D[3][4] = 2.0*Ji*W2*(B(0,1)*B(1,2) - B(0,2)*B(1,1));

	// D[3][5] = 0.5*(c(0,1,0,2) + c(0,1,2,0))
	D[3][5] = 2.0*Ji*W2*(B(0,1)*B(0,2) - B(0,0)*B(1,2));

	// D[4][5] = 0.5*(c(1,2,0,2) + c(1,2,2,0))
	D[4][5] = 2.0*Ji*W2*(B(1,2)*B(0,2) - B(0,1)*B(2,2));

	// --- F I B E R   C O N T R I B U T I O N ---

	// get the element's local coordinate system
	mat3d& Q = pt.Q;

	// loop over all integration points
	double ksi, beta;
	vec3d nr, n0, nt;
	double In, Wl, Wll;
	const double eps = 0;
	for (int n=0; n<m_nint; ++n)
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
			double WC = Wl*In;

			// calculate C:d2WdCdC:C
			double CWWC = Wll*In*In;

			// D[0][0] = c(0,0,0,0)
			D[0][0] += m_w[n]*(8.0/9.0)*Ji*WC;
			D[0][0] += m_w[n]*4.0*Ji*Wll*nt.x*nt.x*nt.x*nt.x;
			D[0][0] += m_w[n]*(4.0/9.0)*Ji*CWWC;
			D[0][0] -= m_w[n]*(8.0/3.0)*Ji*(Wll*In*nt.x*nt.x);

			// D[1][1] = c(1,1,1,1)
			D[1][1] += m_w[n]*8.0/9.0*Ji*WC;
			D[1][1] += m_w[n]*4.0*Ji*Wll*nt.y*nt.y*nt.y*nt.y;
			D[1][1] += m_w[n]*(4.0/9.0)*Ji*CWWC;
			D[1][1] -= m_w[n]*(8.0/3.0)*Ji*(Wll*In*nt.y*nt.y);

			// D[2][2] = c(2,2,2,2)
			D[2][2] += m_w[n]*(8.0/9.0)*Ji*WC;
			D[2][2] += m_w[n]*4.0*Ji*Wll*nt.z*nt.z*nt.z*nt.z;
			D[2][2] += m_w[n]*(4.0/9.0)*Ji*CWWC;
			D[2][2] -= m_w[n]*(8.0/3.0)*Ji*(Wll*In*nt.z*nt.z);

			// D[0][1] = D[1][0] = c(0,0,1,1)
			D[0][1] -= m_w[n]*(4.0/9.0)*Ji*WC;
			D[0][1] += m_w[n]*4.0*Ji*Wll*nt.x*nt.x*nt.y*nt.y;
			D[0][1] += m_w[n]*(4.0/9.0)*Ji*CWWC;
			D[0][1] -= m_w[n]*(4.0/3.0)*Ji*(Wll*In*nt.x*nt.x);
			D[0][1] -= m_w[n]*(4.0/3.0)*Ji*(Wll*In*nt.y*nt.y);

			// D[1][2] = D[2][1] = c(1,1,2,2)
			D[1][2] -= m_w[n]*(4.0/9.0)*Ji*WC;
			D[1][2] += m_w[n]*4.0*Ji*Wll*nt.y*nt.y*nt.z*nt.z;
			D[1][2] += m_w[n]*(4.0/9.0)*Ji*CWWC;
			D[1][2] -= m_w[n]*(4.0/3.0)*Ji*(Wll*In*nt.y*nt.y);
			D[1][2] -= m_w[n]*(4.0/3.0)*Ji*(Wll*In*nt.z*nt.z);

			// D[0][2] = D[2][0] = c(0,0,2,2)
			D[0][2] -= m_w[n]*(4.0/9.0)*Ji*WC;
			D[0][2] += m_w[n]*4.0*Ji*Wll*nt.x*nt.x*nt.z*nt.z;
			D[0][2] += m_w[n]*(4.0/9.0)*Ji*CWWC;
			D[0][2] -= m_w[n]*(4.0/3.0)*Ji*(Wll*In*nt.x*nt.x);
			D[0][2] -= m_w[n]*(4.0/3.0)*Ji*(Wll*In*nt.z*nt.z);



			// D[3][3] = 0.5*(c(0,1,0,1) + c(0,1,1,0))
			D[3][3] += m_w[n]*(2.0/3.0)*Ji*WC;
			D[3][3] += m_w[n]*4.0*Ji*Wll*nt.x*nt.y*nt.x*nt.y;

			// D[4][4] = 0.5*(c(1,2,1,2) + c(1,2,2,1))
			D[4][4] += m_w[n]*(2.0/3.0)*Ji*WC;
			D[4][4] += m_w[n]*4.0*Ji*Wll*nt.y*nt.z*nt.y*nt.z;

			// D[5][5] = 0.5*(c(0,2,0,2) + c(0,2,2,0))
			D[5][5] += m_w[n]*(2.0/3.0)*Ji*WC;
			D[5][5] += m_w[n]*4.0*Ji*Wll*nt.x*nt.z*nt.x*nt.z;



			// D[0][3] = 0.5*(c(0,0,0,1) + c(0,0,1,0))
			D[0][3] += m_w[n]*4.0*Ji*Wll*nt.x*nt.x*nt.x*nt.y;
			D[0][3] -= m_w[n]*(4.0/3.0)*Ji*(Wll*In*nt.x*nt.y);

			// D[0][4] = 0.5*(c(0,0,1,2) + c(0,0,2,1))
			D[0][4] += m_w[n]*4.0*Ji*Wll*nt.x*nt.x*nt.y*nt.z;
			D[0][4] -= m_w[n]*(4.0/3.0)*Ji*(Wll*In*nt.y*nt.z);

			// D[0][5] = 0.5*(c(0,0,0,2) + c(0,0,2,0))
			D[0][5] += m_w[n]*4.0*Ji*Wll*nt.x*nt.x*nt.x*nt.z;
			D[0][5] -= m_w[n]*(4.0/3.0)*Ji*(Wll*In*nt.x*nt.z);

			// D[1][3] = 0.5*(c(1,1,0,1) + c(1,1,1,0))
			D[1][3] += m_w[n]*4.0*Ji*Wll*nt.y*nt.y*nt.x*nt.y;
			D[1][3] -= m_w[n]*(4.0/3.0)*Ji*(Wll*In*nt.x*nt.y);

			// D[1][4] = 0.5*(c(1,1,1,2) + c(1,1,2,1))
			D[1][4] += m_w[n]*4.0*Ji*Wll*nt.y*nt.y*nt.y*nt.z;
			D[1][4] -= m_w[n]*(4.0/3.0)*Ji*(Wll*In*nt.y*nt.z);

			// D[1][5] = 0.5*(c(1,1,0,2) + c(1,1,2,0))
			D[1][5] += m_w[n]*4.0*Ji*Wll*nt.y*nt.y*nt.x*nt.z;
			D[1][5] -= m_w[n]*(4.0/3.0)*Ji*(Wll*In*nt.x*nt.z);

			// D[2][3] = 0.5*(c(2,2,0,1) + c(2,2,1,0))
			D[2][3] += m_w[n]*4.0*Ji*Wll*nt.z*nt.z*nt.x*nt.y;
			D[2][3] -= m_w[n]*(4.0/3.0)*Ji*(Wll*In*nt.x*nt.y);

			// D[2][4] = 0.5*(c(2,2,1,2) + c(2,2,2,1))
			D[2][4] += m_w[n]*4.0*Ji*Wll*nt.z*nt.z*nt.y*nt.z;
			D[2][4] -= m_w[n]*(4.0/3.0)*Ji*(Wll*In*nt.y*nt.z);

			// D[2][5] = 0.5*(c(2,2,0,2) + c(2,2,2,0))
			D[2][5] += m_w[n]*4.0*Ji*Wll*nt.z*nt.z*nt.x*nt.z;
			D[2][5] -= m_w[n]*(4.0/3.0)*Ji*(Wll*In*nt.x*nt.z);



			// D[3][4] = 0.5*(c(0,1,1,2) + c(0,1,2,1))
			D[3][4] += m_w[n]*4.0*Ji*Wll*nt.x*nt.y*nt.y*nt.z;

			// D[3][5] = 0.5*(c(0,1,0,2) + c(0,1,2,0))
			D[3][5] += m_w[n]*4.0*Ji*Wll*nt.x*nt.y*nt.x*nt.z;

			// D[4][5] = 0.5*(c(1,2,0,2) + c(1,2,2,0))
			D[4][5] += m_w[n]*4.0*Ji*Wll*nt.y*nt.z*nt.x*nt.z;
		}
	}

	// set symmetric components
	D[1][0] = D[0][1]; D[2][0] = D[0][2]; D[3][0] = D[0][3]; D[4][0] = D[0][4]; D[5][0] = D[0][5];
	D[2][1] = D[1][2]; D[3][1] = D[1][3]; D[4][1] = D[1][4]; D[5][1] = D[1][5];
	D[3][2] = D[2][3]; D[4][2] = D[2][4]; D[5][2] = D[2][5];
	D[4][3] = D[3][4]; D[5][3] = D[3][5];
	D[5][4] = D[4][5];
}
