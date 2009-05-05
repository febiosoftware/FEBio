#include "stdafx.h"
#include "FERandomFiberNeoHookean.h"

#include "geodesic.h"

// we store the cos and sin of the angles here
int FERandomFiberNeoHookean::m_nres = 0;
double FERandomFiberNeoHookean::m_cth[NSTH];
double FERandomFiberNeoHookean::m_sth[NSTH];
double FERandomFiberNeoHookean::m_cph[NSTH];
double FERandomFiberNeoHookean::m_sph[NSTH];

// register the material with the framework
REGISTER_MATERIAL(FERandomFiberNeoHookean, "random fiber neo-Hookean");

// define the material parameters
BEGIN_PARAMETER_LIST(FERandomFiberNeoHookean, FEElasticMaterial)
	ADD_PARAMETER(m_E, FE_PARAM_DOUBLE, "E");
	ADD_PARAMETER(m_v, FE_PARAM_DOUBLE, "v");
	ADD_PARAMETERV(m_beta, FE_PARAM_DOUBLEV, 3, "beta");
	ADD_PARAMETERV(m_ksi , FE_PARAM_DOUBLEV, 3, "ksi" );
END_PARAMETER_LIST();

#ifndef SQR
	#define SQR(x) ((x)*(x))
#endif


//////////////////////////////////////////////////////////////////////
// FERandomFiberNeoHookean
//////////////////////////////////////////////////////////////////////

void FERandomFiberNeoHookean::Init()
{
	FEElasticMaterial::Init();

	if (m_E <= 0) throw MaterialError("Invalid value for E");
	if (!INRANGE(m_v, -1.0, 0.5)) throw MaterialError("Invalid value for v");

	if (m_ksi[0] <= 0) throw MaterialError("ksi1 must be positive.");
	if (m_ksi[1] <= 0) throw MaterialError("ksi2 must be positive.");
	if (m_ksi[2] <= 0) throw MaterialError("ksi3 must be positive.");
	if (m_beta[0] <= 2) throw MaterialError("beta1 must be bigger than 2.");
	if (m_beta[1] <= 2) throw MaterialError("beta1 must be bigger than 2.");
	if (m_beta[2] <= 2) throw MaterialError("beta1 must be bigger than 2.");

	static bool bfirst = true;

	if (bfirst)
	{
		// select the integration rule
		const int nint    = (m_nres == 0? NSTL  : NSTH  );
		const double* phi = (m_nres == 0? PHIL  : PHIH  );
		const double* the = (m_nres == 0? THETAL: THETAH);
		const double* w   = (m_nres == 0? AREAL : AREAH );

		for (int n=0; n<nint; ++n)
		{
			m_cth[n] = cos(the[n]); 
			m_sth[n] = sin(the[n]);
			m_cph[n] = cos(phi[n]);
			m_sph[n] = sin(phi[n]);
		}

		bfirst = false;
	}
}

mat3ds FERandomFiberNeoHookean::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	mat3d &F = pt.F;
	double detF = pt.J;
	double detFi = 1.0/detF;
	double lndetF = log(detF);

	// calculate left Cauchy-Green tensor
	// (we commented out the matrix components we do not need)
	double b[3][3], C[3][3];

	b[0][0] = F[0][0]*F[0][0]+F[0][1]*F[0][1]+F[0][2]*F[0][2];
	b[0][1] = F[0][0]*F[1][0]+F[0][1]*F[1][1]+F[0][2]*F[1][2];
	b[0][2] = F[0][0]*F[2][0]+F[0][1]*F[2][1]+F[0][2]*F[2][2];

//	b[1][0] = F[1][0]*F[0][0]+F[1][1]*F[0][1]+F[1][2]*F[0][2];
	b[1][1] = F[1][0]*F[1][0]+F[1][1]*F[1][1]+F[1][2]*F[1][2];
	b[1][2] = F[1][0]*F[2][0]+F[1][1]*F[2][1]+F[1][2]*F[2][2];

//	b[2][0] = F[2][0]*F[0][0]+F[2][1]*F[0][1]+F[2][2]*F[0][2];
//	b[2][1] = F[2][0]*F[1][0]+F[2][1]*F[1][1]+F[2][2]*F[1][2];
	b[2][2] = F[2][0]*F[2][0]+F[2][1]*F[2][1]+F[2][2]*F[2][2];

	C[0][0] = F[0][0]*F[0][0]+F[1][0]*F[1][0]+F[2][0]*F[2][2];
	C[0][1] = F[0][0]*F[0][1]+F[1][0]*F[1][1]+F[2][0]*F[2][1];
	C[0][2] = F[0][0]*F[0][2]+F[1][0]*F[1][2]+F[2][0]*F[2][0];

	C[1][0] = F[0][1]*F[0][0]+F[1][1]*F[1][0]+F[2][1]*F[2][0];
	C[1][1] = F[0][1]*F[0][1]+F[1][1]*F[1][1]+F[2][1]*F[2][1];
	C[1][2] = F[0][1]*F[0][2]+F[1][1]*F[1][2]+F[2][1]*F[2][2];

	C[2][0] = F[0][2]*F[0][0]+F[1][2]*F[1][0]+F[2][2]*F[2][0];
	C[2][1] = F[0][2]*F[0][1]+F[1][2]*F[1][1]+F[2][2]*F[2][1];
	C[2][2] = F[0][2]*F[0][2]+F[1][2]*F[1][2]+F[2][2]*F[2][2];


	// --- F I B E R   C O N T R I B U T I O N ---

	// select the integration rule
	const int nint    = (m_nres == 0? NSTL  : NSTH  );
	const double* phi = (m_nres == 0? PHIL  : PHIH  );
	const double* the = (m_nres == 0? THETAL: THETAH);
	const double* w   = (m_nres == 0? AREAL : AREAH );

	// get the element's local coordinate system
	mat3d& Q = pt.Q;

	double T[3][3] = {0};

	// loop over all integration points
	double ksi, beta;
	double nr[3], n0[3], nt[3];
	double In, Wl;
	const double eps = 0;
	for (int n=0; n<nint; ++n)
	{
		// set the local fiber direction
		nr[0] = m_cth[n]*m_sph[n];
		nr[1] = m_sth[n]*m_sph[n];
		nr[2] = m_cph[n];

		// get the global material fiber direction
		n0[0] = Q[0][0]*nr[0] + Q[0][1]*nr[1] + Q[0][2]*nr[2];
		n0[1] = Q[1][0]*nr[0] + Q[1][1]*nr[1] + Q[1][2]*nr[2];
		n0[2] = Q[2][0]*nr[0] + Q[2][1]*nr[1] + Q[2][2]*nr[2];

		// Calculate In = nr*C*nr
		In = n0[0]*( C[0][0]*n0[0] + C[0][1]*n0[1] + C[0][2]*n0[2]) +
			 n0[1]*( C[1][0]*n0[0] + C[1][1]*n0[1] + C[1][2]*n0[2]) +
			 n0[2]*( C[2][0]*n0[0] + C[2][1]*n0[1] + C[2][2]*n0[2]);

		// only take fibers in tension into consideration
		if (In > 1. + eps)
		{
			// calculate material coefficients
			ksi  = 1.0 / sqrt(SQR(nr[0] / m_ksi [0]) + SQR(nr[1] / m_ksi [1]) + SQR(nr[2] / m_ksi [2]));
			beta = 1.0 / sqrt(SQR(nr[0] / m_beta[0]) + SQR(nr[1] / m_beta[1]) + SQR(nr[2] / m_beta[2]));

			// calculate strain energy derivative
			Wl = beta*ksi*pow(In - 1.0, beta-1.0);

			// get the global spatial fiber direction
			nt[0] = F[0][0]*n0[0] + F[0][1]*n0[1] + F[0][2]*n0[2];
			nt[1] = F[1][0]*n0[0] + F[1][1]*n0[1] + F[1][2]*n0[2];
			nt[2] = F[2][0]*n0[0] + F[2][1]*n0[1] + F[2][2]*n0[2];

			// calculate the stress
			T[0][0] += Wl*nt[0]*nt[0]*w[n];
			T[0][1] += Wl*nt[0]*nt[1]*w[n];
			T[0][2] += Wl*nt[0]*nt[2]*w[n];

//			T[1][0] += Wl*nt[1]*nt[0]*w[n];
			T[1][1] += Wl*nt[1]*nt[1]*w[n];
			T[1][2] += Wl*nt[1]*nt[2]*w[n];

//			T[2][0] += Wl*nt[2]*nt[0]*w[n];
//			T[2][1] += Wl*nt[2]*nt[1]*w[n];
			T[2][2] += Wl*nt[2]*nt[2]*w[n];
		}
	}


	// lame parameters
	double lam = m_v*m_E/((1+m_v)*(1-2*m_v));
	double mu  = 0.5*m_E/(1+m_v);

	// calculate stress
	mat3ds s;

	s.xx() = (mu*(b[0][0] - 1) + lam*lndetF)*detFi + 2*detFi*T[0][0];
	s.yy() = (mu*(b[1][1] - 1) + lam*lndetF)*detFi + 2*detFi*T[1][1];
	s.zz() = (mu*(b[2][2] - 1) + lam*lndetF)*detFi + 2*detFi*T[2][2];
	s.xy() = mu*b[0][1]*detFi + 2*detFi*T[0][1];
	s.yz() = mu*b[1][2]*detFi + 2*detFi*T[1][2];
	s.xz() = mu*b[0][2]*detFi + 2*detFi*T[0][2];

	return s;
}

tens4ds FERandomFiberNeoHookean::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	mat3d &F = pt.F;
	double detF = pt.J;
	double detFi = 1.0/detF;
	double lndetF = log(detF);

	// calculate right Cauchy-Green tensor
	// (we commented out the matrix components we do not need)
	double C[3][3];
	C[0][0] = F[0][0]*F[0][0]+F[1][0]*F[1][0]+F[2][0]*F[2][2];
	C[0][1] = F[0][0]*F[0][1]+F[1][0]*F[1][1]+F[2][0]*F[2][1];
	C[0][2] = F[0][0]*F[0][2]+F[1][0]*F[1][2]+F[2][0]*F[2][0];

	C[1][0] = F[0][1]*F[0][0]+F[1][1]*F[1][0]+F[2][1]*F[2][0];
	C[1][1] = F[0][1]*F[0][1]+F[1][1]*F[1][1]+F[2][1]*F[2][1];
	C[1][2] = F[0][1]*F[0][2]+F[1][1]*F[1][2]+F[2][1]*F[2][2];

	C[2][0] = F[0][2]*F[0][0]+F[1][2]*F[1][0]+F[2][2]*F[2][0];
	C[2][1] = F[0][2]*F[0][1]+F[1][2]*F[1][1]+F[2][2]*F[2][1];
	C[2][2] = F[0][2]*F[0][2]+F[1][2]*F[1][2]+F[2][2]*F[2][2];

	// lame parameters
	double lam = m_v*m_E/((1+m_v)*(1-2*m_v));
	double mu  = 0.5*m_E/(1+m_v);

	double lam1 = lam / detF;
	double mu1  = (mu - lam*log(detF)) / detF;

	double D[6][6] = {0};

	D[0][0] = lam1+2.*mu1; D[0][1] = lam1       ; D[0][2] = lam1       ;
	D[1][0] = lam1       ; D[1][1] = lam1+2.*mu1; D[1][2] = lam1       ;
	D[2][0] = lam1       ; D[2][1] = lam1       ; D[2][2] = lam1+2.*mu1;
	D[3][3] = mu1;
	D[4][4] = mu1;
	D[5][5] = mu1;

	// --- F I B E R   C O N T R I B U T I O N ---

	// select the integration rule
	const int nint    = (m_nres == 0? NSTL  : NSTH  );
	const double* phi = (m_nres == 0? PHIL  : PHIH  );
	const double* the = (m_nres == 0? THETAL: THETAH);
	const double* w   = (m_nres == 0? AREAL : AREAH );

	// get the element's local coordinate system
	mat3d& Q = pt.Q;

	// loop over all integration points
	double ksi, beta;
	double nr[3], n0[3], nt[3];
	double In, Wl, Wll;
	const double eps = 0;
	for (int n=0; n<nint; ++n)
	{
		// set the local fiber direction
		nr[0] = m_cth[n]*m_sph[n];
		nr[1] = m_sth[n]*m_sph[n];
		nr[2] = m_cph[n];

		// get the global material fiber direction
		n0[0] = Q[0][0]*nr[0] + Q[0][1]*nr[1] + Q[0][2]*nr[2];
		n0[1] = Q[1][0]*nr[0] + Q[1][1]*nr[1] + Q[1][2]*nr[2];
		n0[2] = Q[2][0]*nr[0] + Q[2][1]*nr[1] + Q[2][2]*nr[2];

		// Calculate In = nr*C*nr
		In = n0[0]*( C[0][0]*n0[0] + C[0][1]*n0[1] + C[0][2]*n0[2]) +
			 n0[1]*( C[1][0]*n0[0] + C[1][1]*n0[1] + C[1][2]*n0[2]) +
			 n0[2]*( C[2][0]*n0[0] + C[2][1]*n0[1] + C[2][2]*n0[2]);

		// only take fibers in tension into consideration
		if (In > 1. + eps)
		{
			// calculate material coefficients
			ksi  = 1.0 / sqrt(SQR(nr[0] / m_ksi [0]) + SQR(nr[1] / m_ksi [1]) + SQR(nr[2] / m_ksi [2]));
			beta = 1.0 / sqrt(SQR(nr[0] / m_beta[0]) + SQR(nr[1] / m_beta[1]) + SQR(nr[2] / m_beta[2]));

			// calculate strain energy derivative
			Wl  = beta*ksi*pow(In - 1.0, beta - 1.0);
			Wll = beta*(beta-1.0)*ksi*pow(In - 1.0, beta-2.0);

			// get the global spatial fiber direction
			nt[0] = F[0][0]*n0[0] + F[0][1]*n0[1] + F[0][2]*n0[2];
			nt[1] = F[1][0]*n0[0] + F[1][1]*n0[1] + F[1][2]*n0[2];
			nt[2] = F[2][0]*n0[0] + F[2][1]*n0[1] + F[2][2]*n0[2];

			D[0][0] += w[n]*4.0*detFi*Wll*nt[0]*nt[0]*nt[0]*nt[0];
			D[1][1] += w[n]*4.0*detFi*Wll*nt[1]*nt[1]*nt[1]*nt[1];
			D[2][2] += w[n]*4.0*detFi*Wll*nt[2]*nt[2]*nt[2]*nt[2];
			D[0][1] += w[n]*4.0*detFi*Wll*nt[0]*nt[0]*nt[1]*nt[1];
			D[1][2] += w[n]*4.0*detFi*Wll*nt[1]*nt[1]*nt[2]*nt[2];
			D[0][2] += w[n]*4.0*detFi*Wll*nt[0]*nt[0]*nt[2]*nt[2];
			D[3][3] += w[n]*4.0*detFi*Wll*nt[0]*nt[1]*nt[0]*nt[1];
			D[4][4] += w[n]*4.0*detFi*Wll*nt[1]*nt[2]*nt[1]*nt[2];
			D[5][5] += w[n]*4.0*detFi*Wll*nt[0]*nt[2]*nt[0]*nt[2];
			D[0][3] += w[n]*4.0*detFi*Wll*nt[0]*nt[0]*nt[0]*nt[1];
			D[0][4] += w[n]*4.0*detFi*Wll*nt[0]*nt[0]*nt[1]*nt[2];
			D[0][5] += w[n]*4.0*detFi*Wll*nt[0]*nt[0]*nt[0]*nt[2];
			D[1][3] += w[n]*4.0*detFi*Wll*nt[1]*nt[1]*nt[0]*nt[1];
			D[1][4] += w[n]*4.0*detFi*Wll*nt[1]*nt[1]*nt[1]*nt[2];
			D[1][5] += w[n]*4.0*detFi*Wll*nt[1]*nt[1]*nt[0]*nt[2];
			D[2][3] += w[n]*4.0*detFi*Wll*nt[2]*nt[2]*nt[0]*nt[1];
			D[2][4] += w[n]*4.0*detFi*Wll*nt[2]*nt[2]*nt[1]*nt[2];
			D[2][5] += w[n]*4.0*detFi*Wll*nt[2]*nt[2]*nt[0]*nt[2];
			D[3][4] += w[n]*4.0*detFi*Wll*nt[0]*nt[1]*nt[1]*nt[2];
			D[3][5] += w[n]*4.0*detFi*Wll*nt[0]*nt[1]*nt[0]*nt[2];
			D[4][5] += w[n]*4.0*detFi*Wll*nt[1]*nt[2]*nt[0]*nt[2];
		}
	}

	// set symmetric components
	D[1][0] = D[0][1]; D[2][0] = D[0][2]; D[3][0] = D[0][3]; D[4][0] = D[0][4]; D[5][0] = D[0][5];
	D[2][1] = D[1][2]; D[3][1] = D[1][3]; D[4][1] = D[1][4]; D[5][1] = D[1][5];
	D[3][2] = D[2][3]; D[4][2] = D[2][4]; D[5][2] = D[2][5];
	D[4][3] = D[3][4]; D[5][3] = D[3][5];
	D[5][4] = D[4][5];

	return tens4ds(D);
}
