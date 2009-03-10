#include "StdAfx.h"
#include "FEFiberNeoHookean.h"

// The following file contains the integration points and weights
// for the integration over a unit sphere in spherical coordinates
#include "geodesic.h"

// we store the cos and sin of the angles here
int FEFiberNeoHookean::m_nres = 0;
double FEFiberNeoHookean::m_cth[NSTH];
double FEFiberNeoHookean::m_sth[NSTH];
double FEFiberNeoHookean::m_cph[NSTH];
double FEFiberNeoHookean::m_sph[NSTH];

// register the material with the framework
REGISTER_MATERIAL(FEFiberNeoHookean, "fiber neo-Hookean");

// define the material parameters
BEGIN_PARAMETER_LIST(FEFiberNeoHookean, FEElasticMaterial)
	ADD_PARAMETER(m_E, FE_PARAM_DOUBLE, "E");
	ADD_PARAMETER(m_v, FE_PARAM_DOUBLE, "v");
	ADD_PARAMETERV(m_a, FE_PARAM_DOUBLEV, 3, "a");
	ADD_PARAMETER(m_ac, FE_PARAM_DOUBLE, "active_contraction");
END_PARAMETER_LIST();

#ifndef SQR
	#define SQR(x) ((x)*(x))
#endif

//////////////////////////////////////////////////////////////////////
// FEFiberNeoHookean
//////////////////////////////////////////////////////////////////////


FEFiberNeoHookean::FEFiberNeoHookean()
{
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

	m_a[0] = m_a[1] = m_a[2] = 0;
	m_ac = 0;
}

void FEFiberNeoHookean::Init()
{
	FEElasticMaterial::Init();

	if (m_E <= 0) throw MaterialError("Invalid value for E");
	if (!INRANGE(m_v, -1.0, 0.5)) throw MaterialError("Invalid value for v");
}

mat3ds FEFiberNeoHookean::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	mat3d &F = pt.F;
	double detF = pt.J;
	double detFi = 1.0/detF;
	double lndetF = log(detF);

	// calculate left Cauchy-Green tensor
	// (we commented out the matrix components we do not need)
	double b[3][3];

	b[0][0] = F[0][0]*F[0][0]+F[0][1]*F[0][1]+F[0][2]*F[0][2];
	b[0][1] = F[0][0]*F[1][0]+F[0][1]*F[1][1]+F[0][2]*F[1][2];
	b[0][2] = F[0][0]*F[2][0]+F[0][1]*F[2][1]+F[0][2]*F[2][2];

//	b[1][0] = F[1][0]*F[0][0]+F[1][1]*F[0][1]+F[1][2]*F[0][2];
	b[1][1] = F[1][0]*F[1][0]+F[1][1]*F[1][1]+F[1][2]*F[1][2];
	b[1][2] = F[1][0]*F[2][0]+F[1][1]*F[2][1]+F[1][2]*F[2][2];

//	b[2][0] = F[2][0]*F[0][0]+F[2][1]*F[0][1]+F[2][2]*F[0][2];
//	b[2][1] = F[2][0]*F[1][0]+F[2][1]*F[1][1]+F[2][2]*F[1][2];
	b[2][2] = F[2][0]*F[2][0]+F[2][1]*F[2][1]+F[2][2]*F[2][2];

	// lame parameters
	double lam = m_v*m_E/((1+m_v)*(1-2*m_v));
	double mu  = 0.5*m_E/(1+m_v);

	// calculate stress
	mat3ds s;

	s.xx() = (mu*(b[0][0] - 1) + lam*lndetF)*detFi;
	s.yy() = (mu*(b[1][1] - 1) + lam*lndetF)*detFi;
	s.zz() = (mu*(b[2][2] - 1) + lam*lndetF)*detFi;
	s.xy() = mu*b[0][1]*detFi;
	s.yz() = mu*b[1][2]*detFi;
	s.xz() = mu*b[0][2]*detFi;

	// --- F I B E R   C O N T R I B U T I O N ---
	if (m_ac)
	{
		// select the integration rule
		const int nint    = (m_nres == 0? NSTL  : NSTH  );
		const double* phi = (m_nres == 0? PHIL  : PHIH  );
		const double* the = (m_nres == 0? THETAL: THETAH);
		const double* w   = (m_nres == 0? AREAL : AREAH );

		// get the element's local coordinate system
		mat3d& Q = pt.Q;

		// loop over all integration points
		double nr[3], n0[3], nt[3];
		double at;
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

			// get the global spatial fiber direction
			nt[0] = F[0][0]*n0[0] + F[0][1]*n0[1] + F[0][2]*n0[2];
			nt[1] = F[1][0]*n0[0] + F[1][1]*n0[1] + F[1][2]*n0[2];
			nt[2] = F[2][0]*n0[0] + F[2][1]*n0[1] + F[2][2]*n0[2];

			// add active contraction stuff
			at = m_ac *sqrt(SQR(m_a[0]*nr[0]) + SQR(m_a[1]*nr[1]) + SQR(m_a[2]*nr[2]));

			s.xx() += at*nt[0]*nt[0];
			s.yy() += at*nt[1]*nt[1];
			s.zz() += at*nt[2]*nt[2];
			s.xy() += at*nt[0]*nt[1];
			s.yz() += at*nt[1]*nt[2];
			s.xz() += at*nt[0]*nt[2];
		}
	}

	return s;
}

void FEFiberNeoHookean::Tangent(double D[6][6], FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// deformation gradient
	mat3d &F = pt.F;
	double detF = pt.J;

	// lame parameters
	double lam = m_v*m_E/((1+m_v)*(1-2*m_v));
	double mu  = 0.5*m_E/(1+m_v);

	double lam1 = lam / detF;
	double mu1  = (mu - lam*log(detF)) / detF;
	
	D[0][0] = lam1+2.*mu1; D[0][1] = lam1       ; D[0][2] = lam1       ;
	D[1][0] = lam1       ; D[1][1] = lam1+2.*mu1; D[1][2] = lam1       ;
	D[2][0] = lam1       ; D[2][1] = lam1       ; D[2][2] = lam1+2.*mu1;
	D[3][3] = mu1;
	D[4][4] = mu1;
	D[5][5] = mu1;
}
