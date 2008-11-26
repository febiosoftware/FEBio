#include "stdafx.h"
#include "FE2DTransIsoMooneyRivlin.h"

// register the material with the framework
REGISTER_MATERIAL(FE2DTransIsoMooneyRivlin, "2D trans iso Mooney-Rivlin");

// define the material parameters
BEGIN_PARAMETER_LIST(FE2DTransIsoMooneyRivlin, FETransverselyIsotropic)
	ADD_PARAMETER(m_c1, FE_PARAM_DOUBLE, "c1");
	ADD_PARAMETER(m_c2, FE_PARAM_DOUBLE, "c2");
END_PARAMETER_LIST();

double FE2DTransIsoMooneyRivlin::m_cth[FE2DTransIsoMooneyRivlin::NSTEPS];
double FE2DTransIsoMooneyRivlin::m_sth[FE2DTransIsoMooneyRivlin::NSTEPS];

//////////////////////////////////////////////////////////////////////
// FE2DTransIsoMooneyRivlin
//////////////////////////////////////////////////////////////////////

FE2DTransIsoMooneyRivlin::FE2DTransIsoMooneyRivlin() : FETransverselyIsotropic(FE_TISO2D_MOONEY_RIVLIN)
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
}

//-----------------------------------------------------------------------------
//! Calculates the stress for this material.
//! \param pt material point at which to evaluate the stress
mat3ds FE2DTransIsoMooneyRivlin::Stress(FEMaterialPoint& mp)
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

	// left Cauchy-Green tensor and its square
	double B[3][3], B2[3][3];

	// calculate deviatoric left Cauchy-Green tensor
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
	// (we commented out the matrix components we do not need)
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

	// --- M A T R I X   C O N T R I B U T I O N ---

	// strain energy derivatives
	double W1, W2;
	W1 = m_c1;
	W2 = m_c2;

	// T = F*dW/dC*Ft
	double T[3][3];
	double w1pw2i1 = W1 + W2*I1;
	T[0][0] = w1pw2i1*B[0][0] - W2*B2[0][0];
	T[0][1] = w1pw2i1*B[0][1] - W2*B2[0][1];
	T[0][2] = w1pw2i1*B[0][2] - W2*B2[0][2];

//	T[1][0] = w1pw2i1*B[1][0] - W2*B2[1][0];
	T[1][1] = w1pw2i1*B[1][1] - W2*B2[1][1];
	T[1][2] = w1pw2i1*B[1][2] - W2*B2[1][2];

//	T[2][0] = w1pw2i1*B[2][0] - W2*B2[2][0];
//	T[2][1] = w1pw2i1*B[2][1] - W2*B2[2][1];
	T[2][2] = w1pw2i1*B[2][2] - W2*B2[2][2];

	// --- F I B E R   C O N T R I B U T I O N ---

	// Next, we calculate the fiber contribution. For this material
	// the fibers lie randomly in a plane that is perpendicular to the transverse
	// axis. We therefor need to integrate over this plane.
	const double PI = 4.0*atan(1.0);
	double wa = 1.0 / (double) NSTEPS;
	vec3d a0, a, v;
	double lam, lamd, I4, W4;
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

		// Add fiber contribution to T
		T[0][0] += wa*W4*I4*a.x*a.x;
		T[0][1] += wa*W4*I4*a.x*a.y;
		T[0][2] += wa*W4*I4*a.x*a.z;

//		T[1][0] += wa*W4*I4*a.y*a.x;
		T[1][1] += wa*W4*I4*a.y*a.y;
		T[1][2] += wa*W4*I4*a.y*a.z;

//		T[2][0] += wa*W4*I4*a.z*a.x;
//		T[2][1] += wa*W4*I4*a.z*a.y;
		T[2][2] += wa*W4*I4*a.z*a.z;
	}

	// (trace of T)/3
	double trT = (T[0][0] + T[1][1] + T[2][2])/3.0;

	// calculate stress
	mat3ds s(
		p + twoJi*(T[0][0] - trT),
		p + twoJi*(T[1][1] - trT),
		p + twoJi*(T[2][2] - trT),
		twoJi*T[0][1],
		twoJi*T[1][2],
		twoJi*T[0][2]);

	return s;
}

//-----------------------------------------------------------------------------
//! Calculates the elasticity tensor for this material.
//! \param D elasticity tensor
//! \param pt material point at which to evaulate the elasticity tensor
void FE2DTransIsoMooneyRivlin::Tangent(double D[6][6], FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// deformation gradient
	mat3d &F = pt.F;
	double J = pt.J;
	double Jm13 = pow(J, -1.0/3.0);
	double Jm23 = Jm13*Jm13;
	double Ji = 1.0/J;

	// deviatoric cauchy-stress, trs = trace[s]/3
	mat3ds& es = pt.s;
	double s[3][3], trs;
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

	// --- M A T R I X   C O N T R I B U T I O N ---

	// strain energy derivatives
	double W1 = m_c1;
	double W2 = m_c2;

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

	// Next, we add the fiber contribution. Since the fibers lie
	// randomly perpendicular to the transverse axis, we need
	// to integrate over that plane
	const double PI = 4.0*atan(1.0);
	double lam, lamd;
	double I4, W4, W44;
	vec3d a0, a, v;
	double wa = 1.0 / (double) NSTEPS;
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

		// Wi = dW/dIi
		if (lamd >= 1)
		{
			double lamdi = 1.0/lamd;
			double Wl, Wll;
			if (lamd < lam1)
			{
				Wl  = lamdi*c3*(exp(c4*(lamd - 1)) - 1);
				Wll = c3*lamdi*(c4*exp(c4*(lamd - 1)) - lamdi*(exp(c4*(lamd-1))-1));
			}
			else
			{
				double c6 = c3*(exp(c4*(lam1-1))-1) - c5*lam1;
				Wl  = lamdi*(c5*lamd + c6);
				Wll = -c6*lamdi*lamdi;
			}
			W4  = 0.5*lamdi*Wl;
			W44 = 0.25*lamdi*lamdi*(Wll - lamdi*Wl);
		}
		else 
		{
			W4 = 0;
			W44 = 0;
		}

		// calculate dWdC:C
		double WC = W4*I4;

		// calculate C:d2WdCdC:C
		double CWWC = W44*I4*I4;

		// D[0][0] = c(0,0,0,0)
		D[0][0] += wa*(8.0/9.0)*Ji*WC;
		D[0][0] += wa*4.0*Ji*W44*I4*I4*a.x*a.x*a.x*a.x;
		D[0][0] += wa*(4.0/9.0)*Ji*CWWC;
		D[0][0] -= wa*(8.0/3.0)*Ji*(W44*I4*I4*a.x*a.x);

		// D[1][1] = c(1,1,1,1)
		D[1][1] += wa*(8.0/9.0)*Ji*WC;
		D[1][1] += wa*4.0*Ji*W44*I4*I4*a.y*a.y*a.y*a.y;
		D[1][1] += wa*(4.0/9.0)*Ji*CWWC;
		D[1][1] -= wa*(8.0/3.0)*Ji*(W44*I4*I4*a.y*a.y);

		// D[2][2] = c(2,2,2,2)
		D[2][2] += wa*(8.0/9.0)*Ji*WC;
		D[2][2] += wa*4.0*Ji*W44*I4*I4*a.z*a.z*a.z*a.z;
		D[2][2] += wa*(4.0/9.0)*Ji*CWWC;
		D[2][2] -= wa*(8.0/3.0)*Ji*(W44*I4*I4*a.z*a.z);



		// D[0][1] = D[1][0] = c(0,0,1,1)
		D[0][1] -= wa*(4.0/9.0)*Ji*WC;
		D[0][1] += wa*4.0*Ji*W44*I4*I4*a.x*a.x*a.y*a.y;
		D[0][1] += wa*(4.0/9.0)*Ji*CWWC;
		D[0][1] -= wa*(4.0/3.0)*Ji*(W44*I4*I4*a.x*a.x);
		D[0][1] -= wa*(4.0/3.0)*Ji*(W44*I4*I4*a.y*a.y);

		// D[1][2] = D[2][1] = c(1,1,2,2)
		D[1][2] -= wa*(4.0/9.0)*Ji*WC;
		D[1][2] += wa*4.0*Ji*W44*I4*I4*a.y*a.y*a.z*a.z;
		D[1][2] += wa*(4.0/9.0)*Ji*CWWC;
		D[1][2] -= wa*(4.0/3.0)*Ji*(W44*I4*I4*a.y*a.y);
		D[1][2] -= wa*(4.0/3.0)*Ji*(W44*I4*I4*a.z*a.z);

		// D[0][2] = D[2][0] = c(0,0,2,2)
		D[0][2] -= wa*(4.0/9.0)*Ji*WC;
		D[0][2] += wa*4.0*Ji*W44*I4*I4*a.x*a.x*a.z*a.z;
		D[0][2] += wa*(4.0/9.0)*Ji*CWWC;
		D[0][2] -= wa*(4.0/3.0)*Ji*(W44*I4*I4*a.x*a.x);
		D[0][2] -= wa*(4.0/3.0)*Ji*(W44*I4*I4*a.z*a.z);


		// D[3][3] = 0.5*(c(0,1,0,1) + c(0,1,1,0))
		D[3][3] += wa*(2.0/3.0)*Ji*WC;
		D[3][3] += wa*4.0*Ji*W44*I4*I4*a.x*a.y*a.x*a.y;

		// D[4][4] = 0.5*(c(1,2,1,2) + c(1,2,2,1))
		D[4][4] += wa*(2.0/3.0)*Ji*WC;
		D[4][4] += wa*4.0*Ji*W44*I4*I4*a.y*a.z*a.y*a.z;

		// D[5][5] = 0.5*(c(0,2,0,2) + c(0,2,2,0))
		D[5][5] += wa*(2.0/3.0)*Ji*WC;
		D[5][5] += wa*4.0*Ji*W44*I4*I4*a.x*a.z*a.x*a.z;



		// D[0][3] = 0.5*(c(0,0,0,1) + c(0,0,1,0))
		D[0][3] += wa*4.0*Ji*W44*I4*I4*a.x*a.x*a.x*a.y;
		D[0][3] -= wa*(4.0/3.0)*Ji*(W44*I4*I4*a.x*a.y);

		// D[0][4] = 0.5*(c(0,0,1,2) + c(0,0,2,1))
		D[0][4] += wa*4.0*Ji*W44*I4*I4*a.x*a.x*a.y*a.z;
		D[0][4] -= wa*(4.0/3.0)*Ji*(W44*I4*I4*a.y*a.z);

		// D[0][5] = 0.5*(c(0,0,0,2) + c(0,0,2,0))
		D[0][5] += wa*4.0*Ji*W44*I4*I4*a.x*a.x*a.x*a.z;
		D[0][5] -= wa*(4.0/3.0)*Ji*(W44*I4*I4*a.x*a.z);

		// D[1][3] = 0.5*(c(1,1,0,1) + c(1,1,1,0))
		D[1][3] += wa*4.0*Ji*W44*I4*I4*a.y*a.y*a.x*a.y;
		D[1][3] -= wa*(4.0/3.0)*Ji*(W44*I4*I4*a.x*a.y);

		// D[1][4] = 0.5*(c(1,1,1,2) + c(1,1,2,1))
		D[1][4] += wa*4.0*Ji*W44*I4*I4*a.y*a.y*a.y*a.z;
		D[1][4] -= wa*(4.0/3.0)*Ji*(W44*I4*I4*a.y*a.z);

		// D[1][5] = 0.5*(c(1,1,0,2) + c(1,1,2,0))
		D[1][5] += wa*4.0*Ji*W44*I4*I4*a.y*a.y*a.x*a.z;
		D[1][5] -= wa*(4.0/3.0)*Ji*(W44*I4*I4*a.x*a.z);

		// D[2][3] = 0.5*(c(2,2,0,1) + c(2,2,1,0))
		D[2][3] += wa*4.0*Ji*W44*I4*I4*a.z*a.z*a.x*a.y;
		D[2][3] -= wa*(4.0/3.0)*Ji*(W44*I4*I4*a.x*a.y);

		// D[2][4] = 0.5*(c(2,2,1,2) + c(2,2,2,1))
		D[2][4] += wa*4.0*Ji*W44*I4*I4*a.z*a.z*a.y*a.z;
		D[2][4] -= wa*(4.0/3.0)*Ji*(W44*I4*I4*a.y*a.z);

		// D[2][5] = 0.5*(c(2,2,0,2) + c(2,2,2,0))
		D[2][5] += wa*4.0*Ji*W44*I4*I4*a.z*a.z*a.x*a.z;
		D[2][5] -= wa*(4.0/3.0)*Ji*(W44*I4*I4*a.x*a.z);



		// D[3][4] = 0.5*(c(0,1,1,2) + c(0,1,2,1))
		D[3][4] += wa*4.0*Ji*W44*I4*I4*a.x*a.y*a.y*a.z;

		// D[3][5] = 0.5*(c(0,1,0,2) + c(0,1,2,0))
		D[3][5] += wa*4.0*Ji*W44*I4*I4*a.x*a.y*a.x*a.z;

		// D[4][5] = 0.5*(c(1,2,0,2) + c(1,2,2,0))
		D[4][5] += wa*4.0*Ji*W44*I4*I4*a.y*a.z*a.x*a.z;
	}

	// set symmetric components
	D[1][0] = D[0][1]; D[2][0] = D[0][2]; D[3][0] = D[0][3]; D[4][0] = D[0][4]; D[5][0] = D[0][5];
	D[2][1] = D[1][2]; D[3][1] = D[1][3]; D[4][1] = D[1][4]; D[5][1] = D[1][5];
	D[3][2] = D[2][3]; D[4][2] = D[2][4]; D[5][2] = D[2][5];
	D[4][3] = D[3][4]; D[5][3] = D[3][5];
	D[5][4] = D[4][5];
}
