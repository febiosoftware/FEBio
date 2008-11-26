#include "StdAfx.h"
#include "FERandomFiberMaterial.h"

// The following file contains the integration points and weights
// for the integration over a unit sphere in spherical coordinates
#include "geodesic.h"

// we store the cos and sin of the angles here
int FERandomFiberMaterial::m_nres = 0;
double FERandomFiberMaterial::m_cth[NSTH];
double FERandomFiberMaterial::m_sth[NSTH];
double FERandomFiberMaterial::m_cph[NSTH];
double FERandomFiberMaterial::m_sph[NSTH];

// register the material with the framework
REGISTER_MATERIAL(FERandomFiberMaterial, "random fiber Mooney-Rivlin");

// define the material parameters
BEGIN_PARAMETER_LIST(FERandomFiberMaterial, FEIncompressibleMaterial)
	ADD_PARAMETER(m_c1, FE_PARAM_DOUBLE, "c1");
	ADD_PARAMETER(m_c2, FE_PARAM_DOUBLE, "c2");
	ADD_PARAMETERV(m_beta, FE_PARAM_DOUBLEV, 3, "beta");
	ADD_PARAMETERV(m_ksi , FE_PARAM_DOUBLEV, 3, "ksi" );
END_PARAMETER_LIST();

#ifndef SQR
	#define SQR(x) ((x)*(x))
#endif

//////////////////////////////////////////////////////////////////////
// FERandomFiberMaterial
//////////////////////////////////////////////////////////////////////

FERandomFiberMaterial::FERandomFiberMaterial() : FEIncompressibleMaterial(FE_RAND_FIBER_MATERIAL)
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
}

mat3ds FERandomFiberMaterial::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// deformation gradient
	mat3d &F = pt.F;
	double J = pt.J;
	double Ji = 1.0 / J;
	double Jm13 = pow(J, -1.0/3.0);
	double Jm23 = Jm13*Jm13;
	double twoJi = 2.0*Ji;

	// average element pressure
	double p = pt.avgp;

	// calculate deviatoric left Cauchy-Green tensor
	double B[3][3];
	B[0][0] = Jm23*(F[0][0]*F[0][0]+F[0][1]*F[0][1]+F[0][2]*F[0][2]);
	B[0][1] = Jm23*(F[0][0]*F[1][0]+F[0][1]*F[1][1]+F[0][2]*F[1][2]);
	B[0][2] = Jm23*(F[0][0]*F[2][0]+F[0][1]*F[2][1]+F[0][2]*F[2][2]);

	B[1][0] = Jm23*(F[1][0]*F[0][0]+F[1][1]*F[0][1]+F[1][2]*F[0][2]);
	B[1][1] = Jm23*(F[1][0]*F[1][0]+F[1][1]*F[1][1]+F[1][2]*F[1][2]);
	B[1][2] = Jm23*(F[1][0]*F[2][0]+F[1][1]*F[2][1]+F[1][2]*F[2][2]);

	B[2][0] = Jm23*(F[2][0]*F[0][0]+F[2][1]*F[0][1]+F[2][2]*F[0][2]);
	B[2][1] = Jm23*(F[2][0]*F[1][0]+F[2][1]*F[1][1]+F[2][2]*F[1][2]);
	B[2][2] = Jm23*(F[2][0]*F[2][0]+F[2][1]*F[2][1]+F[2][2]*F[2][2]);

	// calculate deviatoric right Cauchy-Green tensor (= transpose of B)
	double C[3][3];
	C[0][0] = B[0][0]; C[0][1] = B[1][0]; C[0][2] = B[2][0];
	C[1][0] = B[0][1]; C[1][1] = B[1][1]; C[1][2] = B[2][1];
	C[2][0] = B[0][2]; C[2][1] = B[1][2]; C[2][2] = B[2][2];

	// calculate square of B
	// (we commented out the matrix components we do not need)
	double B2[3][3];
	B2[0][0] = B[0][0]*B[0][0]+B[0][1]*B[1][0]+B[0][2]*B[2][0];
	B2[0][1] = B[0][0]*B[0][1]+B[0][1]*B[1][1]+B[0][2]*B[2][1];
	B2[0][2] = B[0][0]*B[0][2]+B[0][1]*B[1][2]+B[0][2]*B[2][2];

//	B2[1][0] = B[1][0]*B[0][0]+B[1][1]*B[1][0]+B[1][2]*B[2][0];
	B2[1][1] = B[1][0]*B[0][1]+B[1][1]*B[1][1]+B[1][2]*B[2][1];
	B2[1][2] = B[1][0]*B[0][2]+B[1][1]*B[1][2]+B[1][2]*B[2][2];

//	B2[2][0] = B[2][0]*B[0][0]+B[2][1]*B[1][0]+B[2][2]*B[2][0];
//	B2[2][1] = B[2][0]*B[0][1]+B[2][1]*B[1][1]+B[2][2]*B[2][1];
	B2[2][2] = B[2][0]*B[0][2]+B[2][1]*B[1][2]+B[2][2]*B[2][2];

	// Invariants of B (= invariants of C)
	// Note that these are the invariants of Btilde, not of B!
	double I1, I2;
	I1 = B[0][0]+B[1][1]+B[2][2];
	I2 = 0.5*(I1*I1 - ( B2[0][0] + B2[1][1] + B2[2][2]) );

	// strain energy derivatives
	double W1 = m_c1;
	double W2 = m_c2;

	// --- M A T R I X   C O N T R I B U T I O N ---

	// calculate T = F*dW/dC*Ft
	// (we commented out the matrix components we do not need)
	double T[3][3];
	T[0][0] = (W1 + W2*I1)*B[0][0] - W2*B2[0][0];
	T[0][1] = (W1 + W2*I1)*B[0][1] - W2*B2[0][1];
	T[0][2] = (W1 + W2*I1)*B[0][2] - W2*B2[0][2];

//	T[1][0] = (W1 + W2*I1)*B[1][0] - W2*B2[1][0];
	T[1][1] = (W1 + W2*I1)*B[1][1] - W2*B2[1][1];
	T[1][2] = (W1 + W2*I1)*B[1][2] - W2*B2[1][2];

//	T[2][0] = (W1 + W2*I1)*B[2][0] - W2*B2[2][0];
//	T[2][1] = (W1 + W2*I1)*B[2][1] - W2*B2[2][1];
	T[2][2] = (W1 + W2*I1)*B[2][2] - W2*B2[2][2];

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
	double In, Wl;
	const double eps = 1.0e-9;
	for (int n=0; n<nint; ++n)
	{
		// set the local fiber direction
		nr[0] = m_cth[n]*m_sph[n];
		nr[1] = m_sth[n]*m_sph[n];
		nr[2] = m_cph[n];

		// Calculate In = nr*C*nr
		In = nr[0]*( C[0][0]*nr[0] + C[0][1]*nr[1] + C[0][2]*nr[2]) +
			 nr[1]*( C[1][0]*nr[0] + C[1][1]*nr[1] + C[1][2]*nr[2]) +
			 nr[2]*( C[2][0]*nr[0] + C[2][1]*nr[1] + C[2][2]*nr[2]);

		// only take fibers in tension into consideration
		if (In > 1. + eps)
		{
			// calculate material coefficients
			ksi  = 1.0 / sqrt(SQR(nr[0] / m_ksi [0]) + SQR(nr[1] / m_ksi [1]) + SQR(nr[2] / m_ksi [2]));
			beta = 1.0 / sqrt(SQR(nr[0] / m_beta[0]) + SQR(nr[1] / m_beta[1]) + SQR(nr[2] / m_beta[2]));

			// calculate strain energy derivative
			Wl = beta*ksi*pow(In - 1.0, beta-1.0);

			// get the global material fiber direction
			n0[0] = Q[0][0]*nr[0] + Q[0][1]*nr[1] + Q[0][2]*nr[2];
			n0[1] = Q[1][0]*nr[0] + Q[1][1]*nr[1] + Q[1][2]*nr[2];
			n0[2] = Q[2][0]*nr[0] + Q[2][1]*nr[1] + Q[2][2]*nr[2];

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

	// trace of T/3
	double TrT = (T[0][0] + T[1][1] + T[2][2])/3.0;

	// set the stress
	mat3ds s;
	s.xx() = p + twoJi*(T[0][0] - TrT);
	s.yy() = p + twoJi*(T[1][1] - TrT);
	s.zz() = p + twoJi*(T[2][2] - TrT);
	s.xy() = twoJi*T[0][1];
	s.yz() = twoJi*T[1][2];
	s.xz() = twoJi*T[0][1];

	return s;
}

void FERandomFiberMaterial::Tangent(double D[6][6], FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// deformation gradient
	mat3d &F = pt.F;
	double J = pt.J;
	double Jm13 = pow(J, -1.0/3.0);
	double Jm23 = Jm13*Jm13;
	double Ji = 1.0/J;

	// deviatoric cauchy-stress, trs = trace[s]/3
	double s[3][3], trs;
	mat3ds& es = pt.s;
	s[0][0] = es.xx();
	s[1][1] = es.yy();
	s[2][2] = es.zz();
	s[0][1] = s[1][0] = es.xy();
	s[1][2] = s[2][1] = es.yz();
	s[0][2] = s[2][0] = es.xz();

	trs = (s[0][0] + s[1][1] + s[2][2])/3.0;
	s[0][0] -= trs;	
	s[1][1] -= trs;	
	s[2][2] -= trs;

	// mean pressure
	double p = pt.avgp;

	// deviatoric right Cauchy-Green tensor: C = Ft*F
	double C[3][3];
	C[0][0] = Jm23*(F[0][0]*F[0][0]+F[1][0]*F[1][0]+F[2][0]*F[2][0]);
	C[0][1] = Jm23*(F[0][0]*F[0][1]+F[1][0]*F[1][1]+F[2][0]*F[2][1]);
	C[0][2] = Jm23*(F[0][0]*F[0][2]+F[1][0]*F[1][2]+F[2][0]*F[2][2]);

	C[1][0] = Jm23*(F[0][1]*F[0][0]+F[1][1]*F[1][0]+F[2][1]*F[2][0]);
	C[1][1] = Jm23*(F[0][1]*F[0][1]+F[1][1]*F[1][1]+F[2][1]*F[2][1]);
	C[1][2] = Jm23*(F[0][1]*F[0][2]+F[1][1]*F[1][2]+F[2][1]*F[2][2]);

	C[2][0] = Jm23*(F[0][2]*F[0][0]+F[1][2]*F[1][0]+F[2][2]*F[2][0]);
	C[2][1] = Jm23*(F[0][2]*F[0][1]+F[1][2]*F[1][1]+F[2][2]*F[2][1]);
	C[2][2] = Jm23*(F[0][2]*F[0][2]+F[1][2]*F[1][2]+F[2][2]*F[2][2]);

	// square of C
	// (we commented out the components we don't need)
	double C2[3][3];
	C2[0][0] = C[0][0]*C[0][0]+C[0][1]*C[1][0]+C[0][2]*C[2][0];
//	C2[0][1] = C[0][0]*C[0][1]+C[0][1]*C[1][1]+C[0][2]*C[2][1];
//	C2[0][2] = C[0][0]*C[0][2]+C[0][1]*C[1][2]+C[0][2]*C[2][2];

//	C2[1][0] = C[1][0]*C[0][0]+C[1][1]*C[1][0]+C[1][2]*C[2][0];
	C2[1][1] = C[1][0]*C[0][1]+C[1][1]*C[1][1]+C[1][2]*C[2][1];
//	C2[1][2] = C[1][0]*C[0][2]+C[1][1]*C[1][2]+C[1][2]*C[2][2];

//	C2[2][0] = C[2][0]*C[0][0]+C[2][1]*C[1][0]+C[2][2]*C[2][0];
//	C2[2][1] = C[2][0]*C[0][1]+C[2][1]*C[1][1]+C[2][2]*C[2][1];
	C2[2][2] = C[2][0]*C[0][2]+C[2][1]*C[1][2]+C[2][2]*C[2][2];

	// Invariants of C
	double I1 = C[0][0] + C[1][1] + C[2][2];
	double I2 = 0.5*(I1*I1 - (C2[0][0] + C2[1][1] + C2[2][2]));

	// calculate left Cauchy-Green tensor: B = F*Ft
	double B[3][3];
	B[0][0] = Jm23*(F[0][0]*F[0][0]+F[0][1]*F[0][1]+F[0][2]*F[0][2]);
	B[0][1] = Jm23*(F[0][0]*F[1][0]+F[0][1]*F[1][1]+F[0][2]*F[1][2]);
	B[0][2] = Jm23*(F[0][0]*F[2][0]+F[0][1]*F[2][1]+F[0][2]*F[2][2]);

	B[1][0] = Jm23*(F[1][0]*F[0][0]+F[1][1]*F[0][1]+F[1][2]*F[0][2]);
	B[1][1] = Jm23*(F[1][0]*F[1][0]+F[1][1]*F[1][1]+F[1][2]*F[1][2]);
	B[1][2] = Jm23*(F[1][0]*F[2][0]+F[1][1]*F[2][1]+F[1][2]*F[2][2]);

	B[2][0] = Jm23*(F[2][0]*F[0][0]+F[2][1]*F[0][1]+F[2][2]*F[0][2]);
	B[2][1] = Jm23*(F[2][0]*F[1][0]+F[2][1]*F[1][1]+F[2][2]*F[1][2]);
	B[2][2] = Jm23*(F[2][0]*F[2][0]+F[2][1]*F[2][1]+F[2][2]*F[2][2]);

	// calculate square of B
	// (we commented out the components we don't need)
	double B2[3][3];
	B2[0][0] = B[0][0]*B[0][0]+B[0][1]*B[1][0]+B[0][2]*B[2][0];
	B2[0][1] = B[0][0]*B[0][1]+B[0][1]*B[1][1]+B[0][2]*B[2][1];
	B2[0][2] = B[0][0]*B[0][2]+B[0][1]*B[1][2]+B[0][2]*B[2][2];

//	B2[1][0] = B[1][0]*B[0][0]+B[1][1]*B[1][0]+B[1][2]*B[2][0];
	B2[1][1] = B[1][0]*B[0][1]+B[1][1]*B[1][1]+B[1][2]*B[2][1];
	B2[1][2] = B[1][0]*B[0][2]+B[1][1]*B[1][2]+B[1][2]*B[2][2];

//	B2[2][0] = B[2][0]*B[0][0]+B[2][1]*B[1][0]+B[2][2]*B[2][0];
//	B2[2][1] = B[2][0]*B[0][1]+B[2][1]*B[1][1]+B[2][2]*B[2][1];
	B2[2][2] = B[2][0]*B[0][2]+B[2][1]*B[1][2]+B[2][2]*B[2][2];

	// strain energy derivatives
	double W1 = m_c1;
	double W2 = m_c2;

	// --- M A T R I X   C O N T R I B U T I O N ---

	// calculate dWdC:C
	double WC = W1*I1 + 2*W2*I2;

	// calculate C:d2WdCdC:C
	double CWWC = 2*I2*W2;

	// D[0][0] = c(0,0,0,0)
	D[0][0] = -p - (4.0/3.0)*s[0][0] + (8.0/9.0)*Ji*WC;
	D[0][0] -= (8.0/3.0)*Ji*(W2*I1*B[0][0] - W2*B2[0][0]);
	D[0][0] += (4.0/9.0)*Ji*CWWC;

	// D[1][1] = c(1,1,1,1)
	D[1][1] = -p - (4.0/3.0)*s[1][1] + (8.0/9.0)*Ji*WC;
	D[1][1] -= (8.0/3.0)*Ji*(W2*I1*B[1][1] - W2*B2[1][1]);
	D[1][1] += (4.0/9.0)*Ji*CWWC;

	// D[2][2] = c(2,2,2,2)
	D[2][2] = -p - (4.0/3.0)*s[2][2] + (8.0/9.0)*Ji*WC;
	D[2][2] -= (8.0/3.0)*Ji*(W2*I1*B[2][2] - W2*B2[2][2]);
	D[2][2] += (4.0/9.0)*Ji*CWWC;



	// D[0][1] = D[1][0] = c(0,0,1,1)
	D[0][1] = p - (2.0/3.0)*(s[0][0] + s[1][1]) - (4.0/9.0)*Ji*WC;
	D[0][1] += 4.0*Ji*W2*(B[0][0]*B[1][1] - B[0][1]*B[0][1]);
	D[0][1] += (4.0/9.0)*Ji*CWWC;
	D[0][1] -= (4.0/3.0)*Ji*(W2*(I1*B[0][0] - B2[0][0]));
	D[0][1] -= (4.0/3.0)*Ji*(W2*(I1*B[1][1] - B2[1][1]));

	// D[1][2] = D[2][1] = c(1,1,2,2)
	D[1][2] = p - (2.0/3.0)*(s[1][1] + s[2][2]) - (4.0/9.0)*Ji*WC;
	D[1][2] += 4.0*Ji*W2*(B[1][1]*B[2][2] - B[1][2]*B[1][2]);
	D[1][2] += (4.0/9.0)*Ji*CWWC;
	D[1][2] -= (4.0/3.0)*Ji*(W2*(I1*B[1][1] - B2[1][1]));
	D[1][2] -= (4.0/3.0)*Ji*(W2*(I1*B[2][2] - B2[2][2]));

	// D[0][2] = D[2][0] = c(0,0,2,2)
	D[0][2] = p - (2.0/3.0)*(s[0][0] + s[2][2]) - (4.0/9.0)*Ji*WC;
	D[0][2] += 4.0*Ji*W2*(B[0][0]*B[2][2] - B[0][2]*B[0][2]);
	D[0][2] += (4.0/9.0)*Ji*CWWC;
	D[0][2] -= (4.0/3.0)*Ji*(W2*(I1*B[0][0] - B2[0][0]));
	D[0][2] -= (4.0/3.0)*Ji*(W2*(I1*B[2][2] - B2[2][2]));



	// D[3][3] = 0.5*(c(0,1,0,1) + c(0,1,1,0))
	D[3][3] = -p + (2.0/3.0)*Ji*WC;
	D[3][3] += 2.0*Ji*W2*(B[0][1]*B[0][1] - B[0][0]*B[1][1]);

	// D[4][4] = 0.5*(c(1,2,1,2) + c(1,2,2,1))
	D[4][4] = -p + (2.0/3.0)*Ji*WC;
	D[4][4] += 2.0*Ji*W2*(B[1][2]*B[1][2] - B[1][1]*B[2][2]);

	// D[5][5] = 0.5*(c(0,2,0,2) + c(0,2,2,0))
	D[5][5] = -p + (2.0/3.0)*Ji*WC;
	D[5][5] += 2.0*Ji*W2*(B[0][2]*B[0][2] - B[0][0]*B[2][2]);



	// D[0][3] = 0.5*(c(0,0,0,1) + c(0,0,1,0))
	D[0][3] =  -(2.0/3.0)*s[0][1];
	D[0][3] -= (4.0/3.0)*Ji*(W2*(I1*B[0][1] - B2[0][1]));

	// D[0][4] = 0.5*(c(0,0,1,2) + c(0,0,2,1))
	D[0][4] =  -(2.0/3.0)*s[1][2];
	D[0][4] += 4.0*Ji*W2*(B[0][0]*B[1][2] - B[0][1]*B[0][2]);
	D[0][4] -= (4.0/3.0)*Ji*(W2*(I1*B[1][2] - B2[1][2]));

	// D[0][5] = 0.5*(c(0,0,0,2) + c(0,0,2,0))
	D[0][5] =  -(2.0/3.0)*s[0][2];
	D[0][5] -= (4.0/3.0)*Ji*(W2*(I1*B[0][2] - B2[0][2]));

	// D[1][3] = 0.5*(c(1,1,0,1) + c(1,1,1,0))
	D[1][3] =  -(2.0/3.0)*s[0][1];
	D[1][3] -= (4.0/3.0)*Ji*(W2*(I1*B[0][1] - B2[0][1]));

	// D[1][4] = 0.5*(c(1,1,1,2) + c(1,1,2,1))
	D[1][4] =  -(2.0/3.0)*s[1][2];
	D[1][4] -= (4.0/3.0)*Ji*(W2*(I1*B[1][2] - B2[1][2]));

	// D[1][5] = 0.5*(c(1,1,0,2) + c(1,1,2,0))
	D[1][5] =  -(2.0/3.0)*s[0][2];
	D[1][5] += 4.0*Ji*W2*(B[1][1]*B[0][2] - B[0][1]*B[1][2]);
	D[1][5] -= (4.0/3.0)*Ji*(W2*(I1*B[0][2] - B2[0][2]));

	// D[2][3] = 0.5*(c(2,2,0,1) + c(2,2,1,0))
	D[2][3] =  -(2.0/3.0)*s[0][1];
	D[2][3] += 4.0*Ji*W2*(B[2][2]*B[0][1] - B[0][2]*B[1][2]);
	D[2][3] -= (4.0/3.0)*Ji*(W2*(I1*B[0][1] - B2[0][1]));

	// D[2][4] = 0.5*(c(2,2,1,2) + c(2,2,2,1))
	D[2][4] =  -(2.0/3.0)*s[1][2];
	D[2][4] -= (4.0/3.0)*Ji*(W2*(I1*B[1][2] - B2[1][2]));

	// D[2][5] = 0.5*(c(2,2,0,2) + c(2,2,2,0))
	D[2][5] =  -(2.0/3.0)*s[0][2];
	D[2][5] -= (4.0/3.0)*Ji*(W2*(I1*B[0][2] - B2[0][2]));



	// D[3][4] = 0.5*(c(0,1,1,2) + c(0,1,2,1))
	D[3][4] = 2.0*Ji*W2*(B[0][1]*B[1][2] - B[0][2]*B[1][1]);

	// D[3][5] = 0.5*(c(0,1,0,2) + c(0,1,2,0))
	D[3][5] = 2.0*Ji*W2*(B[0][1]*B[0][2] - B[0][0]*B[1][2]);

	// D[4][5] = 0.5*(c(1,2,0,2) + c(1,2,2,0))
	D[4][5] = 2.0*Ji*W2*(B[1][2]*B[0][2] - B[0][1]*B[2][2]);

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
	const double eps = 1.0e-9;
	for (int n=0; n<nint; ++n)
	{
		// set the local fiber direction
		nr[0] = m_cth[n]*m_sph[n];
		nr[1] = m_sth[n]*m_sph[n];
		nr[2] = m_cph[n];

		// Calculate In = nr*C*nr
		In = nr[0]*( C[0][0]*nr[0] + C[0][1]*nr[1] + C[0][2]*nr[2]) +
			 nr[1]*( C[1][0]*nr[0] + C[1][1]*nr[1] + C[1][2]*nr[2]) +
			 nr[2]*( C[2][0]*nr[0] + C[2][1]*nr[1] + C[2][2]*nr[2]);

		// only take fibers in tension into consideration
		if (In > 1. + eps)
		{
			// calculate material coefficients
			ksi  = 1.0 / sqrt(SQR(nr[0] / m_ksi [0]) + SQR(nr[1] / m_ksi [1]) + SQR(nr[2] / m_ksi [2]));
			beta = 1.0 / sqrt(SQR(nr[0] / m_beta[0]) + SQR(nr[1] / m_beta[1]) + SQR(nr[2] / m_beta[2]));

			// calculate strain energy derivative
			Wl  = beta*ksi*pow(In - 1.0, beta - 1.0);
			Wll = beta*(beta-1.0)*ksi*pow(In - 1.0, beta-2.0);

			// get the global material fiber direction
			n0[0] = Q[0][0]*nr[0] + Q[0][1]*nr[1] + Q[0][2]*nr[2];
			n0[1] = Q[1][0]*nr[0] + Q[1][1]*nr[1] + Q[1][2]*nr[2];
			n0[2] = Q[2][0]*nr[0] + Q[2][1]*nr[1] + Q[2][2]*nr[2];

			// get the global spatial fiber direction
			nt[0] = F[0][0]*n0[0] + F[0][1]*n0[1] + F[0][2]*n0[2];
			nt[1] = F[1][0]*n0[0] + F[1][1]*n0[1] + F[1][2]*n0[2];
			nt[2] = F[2][0]*n0[0] + F[2][1]*n0[1] + F[2][2]*n0[2];

			// calculate dWdC:C
			double WC = Wl*In;

			// calculate C:d2WdCdC:C
			double CWWC = Wll*In*In;

			// D[0][0] = c(0,0,0,0)
			D[0][0] += w[n]*(8.0/9.0)*Ji*WC;
			D[0][0] += w[n]*4.0*Ji*Wll*nt[0]*nt[0]*nt[0]*nt[0];
			D[0][0] += w[n]*(4.0/9.0)*Ji*CWWC;
			D[0][0] -= w[n]*(8.0/3.0)*Ji*(Wll*In*nt[0]*nt[0]);

			// D[1][1] = c(1,1,1,1)
			D[1][1] += w[n]*8.0/9.0*Ji*WC;
			D[1][1] += w[n]*4.0*Ji*Wll*nt[1]*nt[1]*nt[1]*nt[1];
			D[1][1] += w[n]*(4.0/9.0)*Ji*CWWC;
			D[1][1] -= w[n]*(8.0/3.0)*Ji*(Wll*In*nt[1]*nt[1]);

			// D[2][2] = c(2,2,2,2)
			D[2][2] += w[n]*(8.0/9.0)*Ji*WC;
			D[2][2] += w[n]*4.0*Ji*Wll*nt[2]*nt[2]*nt[2]*nt[2];
			D[2][2] += w[n]*(4.0/9.0)*Ji*CWWC;
			D[2][2] -= w[n]*(8.0/3.0)*Ji*(Wll*In*nt[2]*nt[2]);

			// D[0][1] = D[1][0] = c(0,0,1,1)
			D[0][1] -= w[n]*(4.0/9.0)*Ji*WC;
			D[0][1] += w[n]*4.0*Ji*Wll*nt[0]*nt[0]*nt[1]*nt[1];
			D[0][1] += w[n]*(4.0/9.0)*Ji*CWWC;
			D[0][1] -= w[n]*(4.0/3.0)*Ji*(Wll*In*nt[0]*nt[0]);
			D[0][1] -= w[n]*(4.0/3.0)*Ji*(Wll*In*nt[1]*nt[1]);

			// D[1][2] = D[2][1] = c(1,1,2,2)
			D[1][2] -= w[n]*(4.0/9.0)*Ji*WC;
			D[1][2] += w[n]*4.0*Ji*Wll*nt[1]*nt[1]*nt[2]*nt[2];
			D[1][2] += w[n]*(4.0/9.0)*Ji*CWWC;
			D[1][2] -= w[n]*(4.0/3.0)*Ji*(Wll*In*nt[1]*nt[1]);
			D[1][2] -= w[n]*(4.0/3.0)*Ji*(Wll*In*nt[2]*nt[2]);

			// D[0][2] = D[2][0] = c(0,0,2,2)
			D[0][2] -= w[n]*(4.0/9.0)*Ji*WC;
			D[0][2] += w[n]*4.0*Ji*Wll*nt[0]*nt[0]*nt[2]*nt[2];
			D[0][2] += w[n]*(4.0/9.0)*Ji*CWWC;
			D[0][2] -= w[n]*(4.0/3.0)*Ji*(Wll*In*nt[0]*nt[0]);
			D[0][2] -= w[n]*(4.0/3.0)*Ji*(Wll*In*nt[2]*nt[2]);



			// D[3][3] = 0.5*(c(0,1,0,1) + c(0,1,1,0))
			D[3][3] += w[n]*(2.0/3.0)*Ji*WC;
			D[3][3] += w[n]*4.0*Ji*Wll*nt[0]*nt[1]*nt[0]*nt[1];

			// D[4][4] = 0.5*(c(1,2,1,2) + c(1,2,2,1))
			D[4][4] += w[n]*(2.0/3.0)*Ji*WC;
			D[4][4] += w[n]*4.0*Ji*Wll*nt[1]*nt[2]*nt[1]*nt[2];

			// D[5][5] = 0.5*(c(0,2,0,2) + c(0,2,2,0))
			D[5][5] += w[n]*(2.0/3.0)*Ji*WC;
			D[5][5] += w[n]*4.0*Ji*Wll*nt[0]*nt[2]*nt[0]*nt[2];



			// D[0][3] = 0.5*(c(0,0,0,1) + c(0,0,1,0))
			D[0][3] += w[n]*4.0*Ji*Wll*nt[0]*nt[0]*nt[0]*nt[1];
			D[0][3] -= w[n]*(4.0/3.0)*Ji*(Wll*In*nt[0]*nt[1]);

			// D[0][4] = 0.5*(c(0,0,1,2) + c(0,0,2,1))
			D[0][4] += w[n]*4.0*Ji*Wll*nt[0]*nt[0]*nt[1]*nt[2];
			D[0][4] -= w[n]*(4.0/3.0)*Ji*(Wll*In*nt[1]*nt[2]);

			// D[0][5] = 0.5*(c(0,0,0,2) + c(0,0,2,0))
			D[0][5] += w[n]*4.0*Ji*Wll*nt[0]*nt[0]*nt[0]*nt[2];
			D[0][5] -= w[n]*(4.0/3.0)*Ji*(Wll*In*nt[0]*nt[2]);

			// D[1][3] = 0.5*(c(1,1,0,1) + c(1,1,1,0))
			D[1][3] += w[n]*4.0*Ji*Wll*nt[1]*nt[1]*nt[0]*nt[1];
			D[1][3] -= w[n]*(4.0/3.0)*Ji*(Wll*In*nt[0]*nt[1]);

			// D[1][4] = 0.5*(c(1,1,1,2) + c(1,1,2,1))
			D[1][4] += w[n]*4.0*Ji*Wll*nt[1]*nt[1]*nt[1]*nt[2];
			D[1][4] -= w[n]*(4.0/3.0)*Ji*(Wll*In*nt[1]*nt[2]);

			// D[1][5] = 0.5*(c(1,1,0,2) + c(1,1,2,0))
			D[1][5] += w[n]*4.0*Ji*Wll*nt[1]*nt[1]*nt[0]*nt[2];
			D[1][5] -= w[n]*(4.0/3.0)*Ji*(Wll*In*nt[0]*nt[2]);

			// D[2][3] = 0.5*(c(2,2,0,1) + c(2,2,1,0))
			D[2][3] += w[n]*4.0*Ji*Wll*nt[2]*nt[2]*nt[0]*nt[1];
			D[2][3] -= w[n]*(4.0/3.0)*Ji*(Wll*In*nt[0]*nt[1]);

			// D[2][4] = 0.5*(c(2,2,1,2) + c(2,2,2,1))
			D[2][4] += w[n]*4.0*Ji*Wll*nt[2]*nt[2]*nt[1]*nt[2];
			D[2][4] -= w[n]*(4.0/3.0)*Ji*(Wll*In*nt[1]*nt[2]);

			// D[2][5] = 0.5*(c(2,2,0,2) + c(2,2,2,0))
			D[2][5] += w[n]*4.0*Ji*Wll*nt[2]*nt[2]*nt[0]*nt[2];
			D[2][5] -= w[n]*(4.0/3.0)*Ji*(Wll*In*nt[0]*nt[2]);



			// D[3][4] = 0.5*(c(0,1,1,2) + c(0,1,2,1))
			D[3][4] = w[n]*4.0*Ji*Wll*nt[0]*nt[1]*nt[1]*nt[2];

			// D[3][5] = 0.5*(c(0,1,0,2) + c(0,1,2,0))
			D[3][5] = w[n]*4.0*Ji*Wll*nt[0]*nt[1]*nt[0]*nt[2];

			// D[4][5] = 0.5*(c(1,2,0,2) + c(1,2,2,0))
			D[4][5] = w[n]*4.0*Ji*Wll*nt[1]*nt[2]*nt[0]*nt[2];
		}
	}

	// set symmetric components
	D[1][0] = D[0][1]; D[2][0] = D[0][2]; D[3][0] = D[0][3]; D[4][0] = D[0][4]; D[5][0] = D[0][5];
	D[2][1] = D[1][2]; D[3][1] = D[1][3]; D[4][1] = D[1][4]; D[5][1] = D[1][5];
	D[3][2] = D[2][3]; D[4][2] = D[2][4]; D[5][2] = D[2][5];
	D[4][3] = D[3][4]; D[5][3] = D[3][5];
	D[5][4] = D[4][5];
}
