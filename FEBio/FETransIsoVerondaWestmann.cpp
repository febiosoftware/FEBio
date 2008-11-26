// FETransIsoVerondaWestmann.cpp: implementation of the FETransIsoVerondaWestmann class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FETransIsoVerondaWestmann.h"

// register the material with the framework
REGISTER_MATERIAL(FETransIsoVerondaWestmann, "trans iso Veronda-Westmann");

// define the material parameters
BEGIN_PARAMETER_LIST(FETransIsoVerondaWestmann, FETransverselyIsotropic)
	ADD_PARAMETER(c1, FE_PARAM_DOUBLE, "c1");
	ADD_PARAMETER(c2, FE_PARAM_DOUBLE, "c2");
END_PARAMETER_LIST();

//////////////////////////////////////////////////////////////////////
// FETransIsoVerondaWestmann
//////////////////////////////////////////////////////////////////////

mat3ds FETransIsoVerondaWestmann::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// deformation gradient
	mat3d& F = pt.F;
	double J = pt.J;
	double Ji = 1.0 / J;
	double Jm13 = pow(J, -1.0/3.0);
	double Jm23 = Jm13*Jm13;
	double twoJi = 2.0*Ji;

	// left Cauchy-Green tensor and its square
	double B[3][3], B2[3][3];

	// get average pressure
	double p = pt.avgp;

	// invariants of B
	double I1, I2, I4;

	// strain energy derivatives
	double W1, W2, W4;

	// T = F*dW/dC*Ft, trT = trace[T]
	double T[3][3], trT;

	// current local material axis
	vec3d a, a0;
	double lam, lamd;

	double w1pw2i1; // = W1 + W2*I1

	const double third = 1.0/3.0;

	// get initial fiber axis
	a0.x = pt.Q[0][0];
	a0.y = pt.Q[1][0];
	a0.z = pt.Q[2][0];

	// calculate the current material axis lam*a = F*a0;
	a.x = F[0][0]*a0.x + F[0][1]*a0.y + F[0][2]*a0.z;
	a.y = F[1][0]*a0.x + F[1][1]*a0.y + F[1][2]*a0.z;
	a.z = F[2][0]*a0.x + F[2][1]*a0.y + F[2][2]*a0.z;

	lam = a.unit();
	lamd = lam*Jm13; // i.e. lambda tilde

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
	I1 = B[0][0]+B[1][1]+B[2][2];
	I2 = 0.5*(I1*I1 - ( B2[0][0] + B2[1][1] + B2[2][2]) );
	I4 = lamd*lamd;

	// --- TODO: put strain energy derivatives here ---
	// Wi = dW/dIi
	W1 = c1*c2*exp(c2*(I1-3));
	W2 = -0.5*c1*c2;
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
	// ---

	// calculate T = F*dW/dC*Ft
	// (we commented out the matrix components we do not need)
	w1pw2i1 = W1 + W2*I1;
	T[0][0] = w1pw2i1*B[0][0] - W2*B2[0][0] + W4*I4*a.x*a.x;
	T[0][1] = w1pw2i1*B[0][1] - W2*B2[0][1] + W4*I4*a.x*a.y;
	T[0][2] = w1pw2i1*B[0][2] - W2*B2[0][2] + W4*I4*a.x*a.z;

//	T[1][0] = w1pw2i1*B[1][0] - W2*B2[1][0] + W4*I4*a.y*a.x;
	T[1][1] = w1pw2i1*B[1][1] - W2*B2[1][1] + W4*I4*a.y*a.y;
	T[1][2] = w1pw2i1*B[1][2] - W2*B2[1][2] + W4*I4*a.y*a.z;

//	T[2][0] = w1pw2i1*B[2][0] - W2*B2[2][0] + W4*I4*a.z*a.x;
//	T[2][1] = w1pw2i1*B[2][1] - W2*B2[2][1] + W4*I4*a.z*a.y;
	T[2][2] = w1pw2i1*B[2][2] - W2*B2[2][2] + W4*I4*a.z*a.z;

	// trace of T/3
	trT = (T[0][0] + T[1][1] + T[2][2])*third;

	// calculate stress: 
	mat3ds s;

	s.xx() = p + twoJi*(T[0][0] - trT);
	s.yy() = p + twoJi*(T[1][1] - trT);
	s.zz() = p + twoJi*(T[2][2] - trT);
	s.xy() = twoJi*T[0][1];
	s.yz() = twoJi*T[1][2];
	s.xz() = twoJi*T[0][2];

	// --- active contraction contribution ---
	if (lcna >= 0)
	{
		double ctenslm = m_plc->Value();

		// current sarcomere length
		double strl = refl*lamd;

		// sarcomere length change
		double dl = strl - l0;

		if (dl >= 0)
		{
			double eca50i = (exp(beta*dl) - 1);

			// active fiber stress
			double saf = ctenslm*eca50i / ( eca50i + 1 );

			s.xx() += saf*a.x*a.x;
			s.yy() += saf*a.y*a.y;
			s.zz() += saf*a.z*a.z;
			s.xy() += saf*a.x*a.y;
			s.yz() += saf*a.y*a.z;
			s.xz() += saf*a.x*a.z;
		}
	}
	//-----------------------------------------

	return s;
}

void FETransIsoVerondaWestmann::Tangent(double D[6][6], FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	const double third = 1.0 / 3.0;
	// deformation gradient
	mat3d& F = pt.F;
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

	trs = (s[0][0] + s[1][1] + s[2][2])*third;
	s[0][0] -= trs;	
	s[1][1] -= trs;	
	s[2][2] -= trs;

	// mean pressure
	double p = pt.avgp;

	// current material axis lam*a = F*a0;
	double lam;
	vec3d a, a0;

	a0.x = pt.Q[0][0];
	a0.y = pt.Q[1][0];
	a0.z = pt.Q[2][0];

	// calculate current material axes
	a.x = F[0][0]*a0.x + F[0][1]*a0.y + F[0][2]*a0.z;
	a.y = F[1][0]*a0.x + F[1][1]*a0.y + F[1][2]*a0.z;
	a.z = F[2][0]*a0.x + F[2][1]*a0.y + F[2][2]*a0.z;

	lam = a.unit();

	// deviatoric stretch
	double lamd = lam*Jm13;

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
	double I4 = lamd*lamd;

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

	// --- TODO: put strain energy derivatives here ---
	// Wi = dW/dIi
	double W1, W2, W11, W4, W44;
	W1 = c1*c2*exp(c2*(I1-3));
	W2 = -0.5*c1*c2;
	W11 = c2*W1;
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
	// ---

	// calculate dWdC:C
	double WC = W1*I1 + 2*W2*I2 + W4*I4;

	// calculate C:d2WdCdC:C
	double CWWC = W11*I1*I1 + 2*I2*W2 + W44*I4*I4;

	// D[0][0] = c(0,0,0,0)
	D[0][0] = -p - (4.0/3.0)*s[0][0] + (8.0/9.0)*Ji*WC;
	D[0][0] += 4.0*Ji*W11*B[0][0]*B[0][0];
	D[0][0] += 4.0*Ji*W44*I4*I4*a.x*a.x*a.x*a.x;
	D[0][0] += (4.0/9.0)*Ji*CWWC;
	D[0][0] -= (8.0/3.0)*Ji*((W11+W2)*I1*B[0][0] - W2*B2[0][0] + W44*I4*I4*a.x*a.x);

	// D[1][1] = c(1,1,1,1)
	D[1][1] = -p - (4.0/3.0)*s[1][1] + (8.0/9.0)*Ji*WC;
	D[1][1] += 4.0*Ji*W11*B[1][1]*B[1][1];
	D[1][1] += 4.0*Ji*W44*I4*I4*a.y*a.y*a.y*a.y;
	D[1][1] += (4.0/9.0)*Ji*CWWC;
	D[1][1] -= (8.0/3.0)*Ji*((W11+W2)*I1*B[1][1] - W2*B2[1][1] + W44*I4*I4*a.y*a.y);

	// D[2][2] = c(2,2,2,2)
	D[2][2] = -p - (4.0/3.0)*s[2][2] + (8.0/9.0)*Ji*WC;
	D[2][2] += 4.0*Ji*W11*B[2][2]*B[2][2];
	D[2][2] += 4.0*Ji*W44*I4*I4*a.z*a.z*a.z*a.z;
	D[2][2] += (4.0/9.0)*Ji*CWWC;
	D[2][2] -= (8.0/3.0)*Ji*((W11+W2)*I1*B[2][2] - W2*B2[2][2] + W44*I4*I4*a.z*a.z);



	// D[0][1] = D[1][0] = c(0,0,1,1)
	D[0][1] = p - (2.0/3.0)*(s[0][0] + s[1][1]) - (4.0/9.0)*Ji*WC;
	D[0][1] += 4.0*Ji*W11*B[0][0]*B[1][1];
	D[0][1] += 4.0*Ji*W44*I4*I4*a.x*a.x*a.y*a.y;
	D[0][1] += 4.0*Ji*W2*(B[0][0]*B[1][1] - B[0][1]*B[0][1]);
	D[0][1] += (4.0/9.0)*Ji*CWWC;
	D[0][1] -= (4.0/3.0)*Ji*((W11+W2)*I1*B[0][0] - W2*B2[0][0] + W44*I4*I4*a.x*a.x);
	D[0][1] -= (4.0/3.0)*Ji*((W11+W2)*I1*B[1][1] - W2*B2[1][1] + W44*I4*I4*a.y*a.y);

	// D[1][2] = D[2][1] = c(1,1,2,2)
	D[1][2] = p - (2.0/3.0)*(s[1][1] + s[2][2]) - (4.0/9.0)*Ji*WC;
	D[1][2] += 4.0*Ji*W11*B[1][1]*B[2][2];
	D[1][2] += 4.0*Ji*W44*I4*I4*a.y*a.y*a.z*a.z;
	D[1][2] += 4.0*Ji*W2*(B[1][1]*B[2][2] - B[1][2]*B[1][2]);
	D[1][2] += (4.0/9.0)*Ji*CWWC;
	D[1][2] -= (4.0/3.0)*Ji*((W11+W2)*I1*B[1][1] - W2*B2[1][1] + W44*I4*I4*a.y*a.y);
	D[1][2] -= (4.0/3.0)*Ji*((W11+W2)*I1*B[2][2] - W2*B2[2][2] + W44*I4*I4*a.z*a.z);

	// D[0][2] = D[2][0] = c(0,0,2,2)
	D[0][2] = p - (2.0/3.0)*(s[0][0] + s[2][2]) - (4.0/9.0)*Ji*WC;
	D[0][2] += 4.0*Ji*W11*B[0][0]*B[2][2];
	D[0][2] += 4.0*Ji*W44*I4*I4*a.x*a.x*a.z*a.z;
	D[0][2] += 4.0*Ji*W2*(B[0][0]*B[2][2] - B[0][2]*B[0][2]);
	D[0][2] += (4.0/9.0)*Ji*CWWC;
	D[0][2] -= (4.0/3.0)*Ji*((W11+W2)*I1*B[0][0] - W2*B2[0][0] + W44*I4*I4*a.x*a.x);
	D[0][2] -= (4.0/3.0)*Ji*((W11+W2)*I1*B[2][2] - W2*B2[2][2] + W44*I4*I4*a.z*a.z);



	// D[3][3] = 0.5*(c(0,1,0,1) + c(0,1,1,0))
	D[3][3] = -p + (2.0/3.0)*Ji*WC;
	D[3][3] += 4.0*Ji*W11*B[0][1]*B[0][1];
	D[3][3] += 4.0*Ji*W44*I4*I4*a.x*a.y*a.x*a.y;
	D[3][3] += 2.0*Ji*W2*(B[0][1]*B[0][1] - B[0][0]*B[1][1]);

	// D[4][4] = 0.5*(c(1,2,1,2) + c(1,2,2,1))
	D[4][4] = -p + (2.0/3.0)*Ji*WC;
	D[4][4] += 4.0*Ji*W11*B[1][2]*B[1][2];
	D[4][4] += 4.0*Ji*W44*I4*I4*a.y*a.z*a.y*a.z;
	D[4][4] += 2.0*Ji*W2*(B[1][2]*B[1][2] - B[1][1]*B[2][2]);

	// D[5][5] = 0.5*(c(0,2,0,2) + c(0,2,2,0))
	D[5][5] = -p + (2.0/3.0)*Ji*WC;
	D[5][5] += 4.0*Ji*W11*B[0][2]*B[0][2];
	D[5][5] += 4.0*Ji*W44*I4*I4*a.x*a.z*a.x*a.z;
	D[5][5] += 2.0*Ji*W2*(B[0][2]*B[0][2] - B[0][0]*B[2][2]);



	// D[0][3] = 0.5*(c(0,0,0,1) + c(0,0,1,0))
	D[0][3] =  -(2.0/3.0)*s[0][1];
	D[0][3] += 4.0*Ji*W11*B[0][0]*B[0][1];
	D[0][3] += 4.0*Ji*W44*I4*I4*a.x*a.x*a.x*a.y;
	D[0][3] -= (4.0/3.0)*Ji*((W11+W2)*I1*B[0][1] - W2*B2[0][1] + W44*I4*I4*a.x*a.y);

	// D[0][4] = 0.5*(c(0,0,1,2) + c(0,0,2,1))
	D[0][4] =  -(2.0/3.0)*s[1][2];
	D[0][4] += 4.0*Ji*W11*B[0][0]*B[1][2];
	D[0][4] += 4.0*Ji*W44*I4*I4*a.x*a.x*a.y*a.z;
	D[0][4] += 4.0*Ji*W2*(B[0][0]*B[1][2] - B[0][1]*B[0][2]);
	D[0][4] -= (4.0/3.0)*Ji*((W11+W2)*I1*B[1][2] - W2*B2[1][2] + W44*I4*I4*a.y*a.z);

	// D[0][5] = 0.5*(c(0,0,0,2) + c(0,0,2,0))
	D[0][5] =  -(2.0/3.0)*s[0][2];
	D[0][5] += 4.0*Ji*W11*B[0][0]*B[0][2];
	D[0][5] += 4.0*Ji*W44*I4*I4*a.x*a.x*a.x*a.z;
	D[0][5] -= (4.0/3.0)*Ji*((W11+W2)*I1*B[0][2] - W2*B2[0][2] + W44*I4*I4*a.x*a.z);

	// D[1][3] = 0.5*(c(1,1,0,1) + c(1,1,1,0))
	D[1][3] =  -(2.0/3.0)*s[0][1];
	D[1][3] += 4.0*Ji*W11*B[1][1]*B[0][1];
	D[1][3] += 4.0*Ji*W44*I4*I4*a.y*a.y*a.x*a.y;
	D[1][3] -= (4.0/3.0)*Ji*((W11+W2)*I1*B[0][1] - W2*B2[0][1] + W44*I4*I4*a.x*a.y);

	// D[1][4] = 0.5*(c(1,1,1,2) + c(1,1,2,1))
	D[1][4] =  -(2.0/3.0)*s[1][2];
	D[1][4] += 4.0*Ji*W11*B[1][1]*B[1][2];
	D[1][4] += 4.0*Ji*W44*I4*I4*a.y*a.y*a.y*a.z;
	D[1][4] -= (4.0/3.0)*Ji*((W11+W2)*I1*B[1][2] - W2*B2[1][2] + W44*I4*I4*a.y*a.z);

	// D[1][5] = 0.5*(c(1,1,0,2) + c(1,1,2,0))
	D[1][5] =  -(2.0/3.0)*s[0][2];
	D[1][5] += 4.0*Ji*W11*B[1][1]*B[0][2];
	D[1][5] += 4.0*Ji*W44*I4*I4*a.y*a.y*a.x*a.z;
	D[1][5] += 4.0*Ji*W2*(B[1][1]*B[0][2] - B[0][1]*B[1][2]);
	D[1][5] -= (4.0/3.0)*Ji*((W11+W2)*I1*B[0][2] - W2*B2[0][2] + W44*I4*I4*a.x*a.z);

	// D[2][3] = 0.5*(c(2,2,0,1) + c(2,2,1,0))
	D[2][3] =  -(2.0/3.0)*s[0][1];
	D[2][3] += 4.0*Ji*W11*B[2][2]*B[0][1];
	D[2][3] += 4.0*Ji*W44*I4*I4*a.z*a.z*a.x*a.y;
	D[2][3] += 4.0*Ji*W2*(B[2][2]*B[0][1] - B[0][2]*B[1][2]);
	D[2][3] -= (4.0/3.0)*Ji*((W11+W2)*I1*B[0][1] - W2*B2[0][1] + W44*I4*I4*a.x*a.y);

	// D[2][4] = 0.5*(c(2,2,1,2) + c(2,2,2,1))
	D[2][4] =  -(2.0/3.0)*s[1][2];
	D[2][4] += 4.0*Ji*W11*B[2][2]*B[1][2];
	D[2][4] += 4.0*Ji*W44*I4*I4*a.z*a.z*a.y*a.z;
	D[2][4] -= (4.0/3.0)*Ji*((W11+W2)*I1*B[1][2] - W2*B2[1][2] + W44*I4*I4*a.y*a.z);

	// D[2][5] = 0.5*(c(2,2,0,2) + c(2,2,2,0))
	D[2][5] =  -(2.0/3.0)*s[0][2];
	D[2][5] += 4.0*Ji*W11*B[2][2]*B[0][2];
	D[2][5] += 4.0*Ji*W44*I4*I4*a.z*a.z*a.x*a.z;
	D[2][5] -= (4.0/3.0)*Ji*((W11+W2)*I1*B[0][2] - W2*B2[0][2] + W44*I4*I4*a.x*a.z);



	// D[3][4] = 0.5*(c(0,1,1,2) + c(0,1,2,1))
	D[3][4] = 4.0*Ji*W44*I4*I4*a.x*a.y*a.y*a.z;
	D[3][4] += 4.0*Ji*W11*B[0][1]*B[1][2];
	D[3][4] += 2.0*Ji*W2*(B[0][1]*B[1][2] - B[0][2]*B[1][1]);

	// D[3][5] = 0.5*(c(0,1,0,2) + c(0,1,2,0))
	D[3][5] = 4.0*Ji*W44*I4*I4*a.x*a.y*a.x*a.z;
	D[3][5] += 4.0*Ji*W11*B[0][1]*B[0][2];
	D[3][5] += 2.0*Ji*W2*(B[0][1]*B[0][2] - B[0][0]*B[1][2]);

	// D[4][5] = 0.5*(c(1,2,0,2) + c(1,2,2,0))
	D[4][5] = 4.0*Ji*W44*I4*I4*a.y*a.z*a.x*a.z;
	D[4][5] += 4.0*Ji*W11*B[1][2]*B[0][2];
	D[4][5] += 2.0*Ji*W2*(B[1][2]*B[0][2] - B[0][1]*B[2][2]);



	// set symmetric components
	D[1][0] = D[0][1]; D[2][0] = D[0][2]; D[3][0] = D[0][3]; D[4][0] = D[0][4]; D[5][0] = D[0][5];
	D[2][1] = D[1][2]; D[3][1] = D[1][3]; D[4][1] = D[1][4]; D[5][1] = D[1][5];
	D[3][2] = D[2][3]; D[4][2] = D[2][4]; D[5][2] = D[2][5];
	D[4][3] = D[3][4]; D[5][3] = D[3][5];
	D[5][4] = D[4][5];
}
