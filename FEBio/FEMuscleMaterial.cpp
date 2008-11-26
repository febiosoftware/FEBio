// FEMuscleMaterial.cpp: implementation of the FEMuscleMaterial class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FEMuscleMaterial.h"

#ifndef SQR
	#define SQR(x) ((x)*(x))
#endif

// register the material with the framework
REGISTER_MATERIAL(FEMuscleMaterial, "muscle material");

// define the material parameters
BEGIN_PARAMETER_LIST(FEMuscleMaterial, FETransverselyIsotropic)
	ADD_PARAMETER(m_G1, FE_PARAM_DOUBLE, "g1");
	ADD_PARAMETER(m_G2, FE_PARAM_DOUBLE, "g2");
	ADD_PARAMETER(m_G3, FE_PARAM_DOUBLE, "g3");
	ADD_PARAMETER(m_P1, FE_PARAM_DOUBLE, "p1");
	ADD_PARAMETER(m_P2, FE_PARAM_DOUBLE, "p2");
	ADD_PARAMETER(m_Lofl, FE_PARAM_DOUBLE, "Lofl");
	ADD_PARAMETER(m_smax, FE_PARAM_DOUBLE, "smax");
END_PARAMETER_LIST();

/////////////////////////////////////////////////////////////////////////
// FEMuscleMaterial
/////////////////////////////////////////////////////////////////////////

mat3ds FEMuscleMaterial::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// deformation gradient
	mat3d &F = pt.F;
	double J = pt.J;
	double Ji = 1.0/J;
	double Jm13 = pow(J, -1.0/3.0);
	double Jm23 = Jm13*Jm13;
	double twoJi = 2.0*Ji;

	// left Cauchy-Green tensor and its square
	double B[3][3], B2[3][3];

	// average element pressure
	double p = pt.avgp;

	// current local material axis
	vec3d a0, a;
	double la, lat;

	// get the initial fiber direction
	a0.x = pt.Q[0][0];
	a0.y = pt.Q[1][0];
	a0.z = pt.Q[2][0];

	// calculate the current material axis lam*a = F*a0;
	a.x = F[0][0]*a0.x + F[0][1]*a0.y + F[0][2]*a0.z;
	a.y = F[1][0]*a0.x + F[1][1]*a0.y + F[1][2]*a0.z;
	a.z = F[2][0]*a0.x + F[2][1]*a0.y + F[2][2]*a0.z;

	// normalize material axis and store fiber stretch
	la = a.unit();
	lat = la*Jm13; // i.e. lambda tilde = sqrt(I4)

	// calculate deviatoric left Cauchy-Green tensor: B = Jm23*F*Ft
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
	B2[0][0] = B[0][0]*B[0][0]+B[0][1]*B[1][0]+B[0][2]*B[2][0];
	B2[0][1] = B[0][0]*B[0][1]+B[0][1]*B[1][1]+B[0][2]*B[2][1];
	B2[0][2] = B[0][0]*B[0][2]+B[0][1]*B[1][2]+B[0][2]*B[2][2];

	B2[1][0] = B[1][0]*B[0][0]+B[1][1]*B[1][0]+B[1][2]*B[2][0];
	B2[1][1] = B[1][0]*B[0][1]+B[1][1]*B[1][1]+B[1][2]*B[2][1];
	B2[1][2] = B[1][0]*B[0][2]+B[1][1]*B[1][2]+B[1][2]*B[2][2];

	B2[2][0] = B[2][0]*B[0][0]+B[2][1]*B[1][0]+B[2][2]*B[2][0];
	B2[2][1] = B[2][0]*B[0][1]+B[2][1]*B[1][1]+B[2][2]*B[2][1];
	B2[2][2] = B[2][0]*B[0][2]+B[2][1]*B[1][2]+B[2][2]*B[2][2];

	// calculate Ba = B*a
	vec3d Ba;
	Ba.x = B[0][0]*a.x + B[0][1]*a.y + B[0][2]*a.z;
	Ba.y = B[1][0]*a.x + B[1][1]*a.y + B[1][2]*a.z;
	Ba.z = B[2][0]*a.x + B[2][1]*a.y + B[2][2]*a.z;

	// ----- invariants of B -----
	double I1, I2, I4, I5;

	I1 = B[0][0]+B[1][1]+B[2][2];
	I2 = 0.5*(I1*I1 - ( B2[0][0] + B2[1][1] + B2[2][2]) );
	I4 = lat*lat;
	I5 = I4*(a*Ba);

	// ----- calculate new invariants b1 and b2 of B ------

	// note that we need to be carefull about roundoff errors
	// since these functions may create numerical problems

	double g = I5/(I4*I4) - 1;
	double b1 = (g > 0 ? sqrt(g) : 0);
	
	double b2 = acosh(0.5*(I1*I4 - I5)/lat);

	// calculate omage (w)
	double I4r = sqrt(I4);
	double w = 0.5*(I1*I4 - I5)/sqrt(I4);

	// set beta and ksi to their limit values
	double beta = 1.0;
	double ksi = -1.0/3.0;

	// if w not equals unity, we can calculate beta and ksi
	if (w > 1.0001)
	{
		beta = b2/sqrt(w*w-1);
		ksi = (1.0/(w*w-1))*(1 - w*b2/sqrt(w*w-1));
	}

	// ----- calculate strain derivatives -----

	// We assume that W(I1, I4, I5, alpha) = F1(B1(I4, I5)) + F2(B2(I1,I4,I5)) + F3(lam(I4), alpha)
	double W1, W2, W4, W5;

	// calculate derivatives for F1
	double F1D4 = -2*m_G1*(I5/(I4*I4*I4));
	double F1D5 = m_G1/(I4*I4);

	// calculate derivatives for F2
	double F2D1 =  m_G2*beta*lat;
	double F2D4 =  m_G2*beta*(I1*I4 + I5)*0.5*pow(I4, -1.5);
	double F2D5 = -m_G2*beta/lat;

	// calculate derivatives for F3
	// these terms are proposed to fix the zero-stress problem
	double F3D4 = 9.0*m_G3*0.125*log(I4)/I4;

	// calculate passive fiber force
	double Fp;
/*	if (lat <= m_Lofl)
	{
		Fp = 0;
	}
	else if (lat < lam1)
*/	if (lat < lam1)
	{
		Fp = m_P1*(exp(m_P2*(lat/m_Lofl - 1)) - 1);
	}
	else
	{
		double P3, P4;

		P3 = m_P1*m_P2*exp(m_P2*(lam1/m_Lofl - 1));
		P4 = m_P1*(exp(m_P2*(lam1/m_Lofl - 1)) - 1) - P3*lam1/m_Lofl;

		Fp = P3*lat/m_Lofl + P4;
	}

	// calculate active fiber force
	double Fa;

	if ((lat <= 0.4*m_Lofl) || (lat >= 1.6*m_Lofl))
	{
		// we have added this part to make sure that 
		// Fa is zero outside the range [0.4, 1.6] *m_Lofl
		Fa = 0;
	}
	else
	{
		if (lat <= 0.6*m_Lofl)
		{
			Fa = 9*SQR(lat/m_Lofl - 0.4);
		}
		else if (lat >= 1.4*m_Lofl)
		{
			Fa = 9*SQR(lat/m_Lofl - 1.6);
		}
		else if ((lat >= 0.6*m_Lofl) && (lat <= 1.4*m_Lofl))
		{
			Fa = 1 - 4*SQR(1 - lat/m_Lofl);
		}
	}

	// activation level
	double alpha = m_ascl*(m_plc ? m_plc->Value():1);

	// calculate total fiber force
	double FfDl = m_smax*(Fp + alpha*Fa)/m_Lofl;
	double FfD4  = 0.5*FfDl/lat;

	// add all derivatives
	W1 = F2D1;
	W2 = 0;
	W4 = F1D4 + F2D4 + F3D4 + FfD4;
	W5 = F1D5 + F2D5;

	// ----- calculate stress -----

	// calculate T 
	double T[3][3];

	T[0][0] = (W1 + W2*I1)*B[0][0] - W2*B2[0][0] + W4*I4*a.x*a.x + W5*I4*(a.x*Ba.x + Ba.x*a.x);
	T[0][1] = (W1 + W2*I1)*B[0][1] - W2*B2[0][1] + W4*I4*a.x*a.y + W5*I4*(a.x*Ba.y + Ba.x*a.y);
	T[0][2] = (W1 + W2*I1)*B[0][2] - W2*B2[0][2] + W4*I4*a.x*a.z + W5*I4*(a.x*Ba.z + Ba.x*a.z);

//	T[1][0] = (W1 + W2*I1)*B[1][0] - W2*B2[1][0] + W4*I4*a.y*a.x + W5*I4*(a.y*Ba.x + Ba.y*a.x);
	T[1][1] = (W1 + W2*I1)*B[1][1] - W2*B2[1][1] + W4*I4*a.y*a.y + W5*I4*(a.y*Ba.y + Ba.y*a.y);
	T[1][2] = (W1 + W2*I1)*B[1][2] - W2*B2[1][2] + W4*I4*a.y*a.z + W5*I4*(a.y*Ba.z + Ba.y*a.z);
 
//	T[2][0] = (W1 + W2*I1)*B[2][0] - W2*B2[2][0] + W4*I4*a.z*a.x + W5*I4*(a.z*Ba.x + Ba.z*a.x);
//	T[2][1] = (W1 + W2*I1)*B[2][1] - W2*B2[2][1] + W4*I4*a.z*a.y + W5*I4*(a.z*Ba.y + Ba.z*a.y);
	T[2][2] = (W1 + W2*I1)*B[2][2] - W2*B2[2][2] + W4*I4*a.z*a.z + W5*I4*(a.z*Ba.z + Ba.z*a.z);

	// calculate 1/3 of trace of T
	double trT = (T[0][0] + T[1][1] + T[2][2]) / 3.0;

	// calculate stress
	mat3ds s;
	s.xx() = p + twoJi*(T[0][0] - trT);
	s.yy() = p + twoJi*(T[1][1] - trT);
	s.zz() = p + twoJi*(T[2][2] - trT);
	s.xy() = twoJi*T[0][1];
	s.yz() = twoJi*T[1][2];
	s.xz() = twoJi*T[0][2];

	return s;
}

void FEMuscleMaterial::Tangent(double D[6][6], FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	const double third = 1.0 / 3.0;

	// deformation gradient
	mat3d &F = pt.F;
	double J = pt.J;
	double Jm13 = pow(J, -1.0/3.0);
	double Jm23 = Jm13*Jm13;
	double Ji = 1.0/J;
	double twoJi = 2.0*Ji;

	// left Cauchy-Green tensor and its square
	double B[3][3], B2[3][3];

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

	// current local material axis
	vec3d a0, a;

	// get the initial fiber direction
	a0.x = pt.Q[0][0];
	a0.y = pt.Q[1][0];
	a0.z = pt.Q[2][0];

	// calculate the current material axis lam*a = F*a0;
	a.x = F[0][0]*a0.x + F[0][1]*a0.y + F[0][2]*a0.z;
	a.y = F[1][0]*a0.x + F[1][1]*a0.y + F[1][2]*a0.z;
	a.z = F[2][0]*a0.x + F[2][1]*a0.y + F[2][2]*a0.z;

	// normalize material axis and store fiber stretch
	double la  = a.unit();
	double lat = la*Jm13; // i.e. lambda tilde

	// calculate deviatoric left Cauchy-Green tensor (use symmetry to reduce FLOP)
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

	B2[1][0] = B[1][0]*B[0][0]+B[1][1]*B[1][0]+B[1][2]*B[2][0];
	B2[1][1] = B[1][0]*B[0][1]+B[1][1]*B[1][1]+B[1][2]*B[2][1];
	B2[1][2] = B[1][0]*B[0][2]+B[1][1]*B[1][2]+B[1][2]*B[2][2];

	B2[2][0] = B[2][0]*B[0][0]+B[2][1]*B[1][0]+B[2][2]*B[2][0];
	B2[2][1] = B[2][0]*B[0][1]+B[2][1]*B[1][1]+B[2][2]*B[2][1];
	B2[2][2] = B[2][0]*B[0][2]+B[2][1]*B[1][2]+B[2][2]*B[2][2];

	// calculate Ba
	vec3d Ba;
	Ba.x = B[0][0]*a.x + B[0][1]*a.y + B[0][2]*a.z;
	Ba.y = B[1][0]*a.x + B[1][1]*a.y + B[1][2]*a.z;
	Ba.z = B[2][0]*a.x + B[2][1]*a.y + B[2][2]*a.z;

	// invariants of B
	double I1, I2, I4, I5;

	I1 = B[0][0]+B[1][1]+B[2][2];
	I2 = 0.5*(I1*I1 - ( B2[0][0] + B2[1][1] + B2[2][2]) );
	I4 = lat*lat;
	I5 = I4*(a*Ba);

	// calculate new invariants
	double g = I5/(I4*I4) - 1;
	double b1 = (g > 0 ? sqrt(g) : 0);
	
	double b2 = acosh(0.5*(I1*I4 - I5)/lat);

	// calculate omage (w)
	double w = 0.5*(I1*I4 - I5)/lat;

	// set beta and ksi to their limit values
	double beta = 1.0;
	double ksi = -1.0/3.0;

	// if w not equals unity, we can calculate beta and ksi
	if (w > 1.0001)
	{
		beta = b2/sqrt(w*w-1);
		ksi = (1.0/(w*w-1))*(1 - w*b2/sqrt(w*w-1));
	}

	// --- strain energy derivatives ---
	// We assume that W(I1, I4, I5, alpha) = F1(B1(I4, I5)) + F2(B2(I1,I4,I5)) + F3(lam(I4), alpha)
	double W1, W2, W4, W5;

	// -- A. matrix contribution --
	// calculate derivatives for F1
	double F1D4 = -2*m_G1*(I5/(I4*I4*I4));
	double F1D5 = m_G1/(I4*I4);

	double F1D44 = 6*m_G1*(I5/(I4*I4*I4*I4));
	double F1D45 = -2*m_G1/(I4*I4*I4);

	// calculate derivatives for F2
	double F2D1 =  m_G2*beta*lat;
	double F2D4 =  m_G2*beta*(I1*I4 + I5)*0.5*pow(I4, -1.5);
	double F2D5 = -m_G2*beta/lat;

	double F2D11 = ksi*m_G2*I4*0.5;
	double F2D44 = 2.0*m_G2*ksi*pow(0.25*(I1*I4+I5)/pow(I4, 1.5), 2) - m_G2*beta*(0.25*(I1*I4 + 3*I5) / pow(I4, 2.5));
	double F2D55 = 0.5*m_G2*ksi/I4;
	double F2D14 = m_G2*beta*0.5/lat + m_G2*ksi*(I1*I4+I5)*0.25/I4;
	double F2D15 = -0.5*m_G2*ksi;
	double F2D45 = m_G2*beta*0.5*pow(I4, -1.5) - m_G2*ksi*0.25*(I1*I4+I5)/(I4*I4);

	// calculate derivatives for F3
	// these terms are proposed to fix the zero-stress problem
	double F3D4  = 9.0*m_G3*0.125*log(I4)/I4;
	double F3D44 = 9.0*m_G3*0.125*(1 - log(I4))/(I4*I4);

	// -- B. fiber contribution --

	// calculate passive fiber force
	double Fp, FpDl;
/*	if (lat <= m_Lofl)
	{
		Fp = 0;
		FpDl = 0;
	}
	else if (lat < lam1)
*/	if (lat < lam1)
	{
		Fp = m_P1*(exp(m_P2*(lat/m_Lofl - 1)) - 1);
		FpDl = m_P1*m_P2*exp(m_P2*(lat/m_Lofl-1))/m_Lofl;
	}
	else
	{
		double P3, P4;

		P3 = m_P1*m_P2*exp(m_P2*(lam1/m_Lofl - 1));
		P4 = m_P1*(exp(m_P2*(lam1/m_Lofl - 1)) - 1) - P3*lam1/m_Lofl;

		Fp = P3*lat/m_Lofl + P4;
		FpDl = P3/m_Lofl;
	}

	// calculate active fiber force
	double Fa, FaDl;

	if ((lat <= 0.4*m_Lofl) || (lat >= 1.6*m_Lofl))
	{
		// we have added this part to make sure that 
		// Fa is zero outside the range [0.4, 1.6] *m_Lofl
		Fa = 0;
		FaDl = 0;
	}
	else
	{
		if (lat <= 0.6*m_Lofl)
		{
			Fa = 9*SQR(lat/m_Lofl - 0.4);
			FaDl = 18*(lat/m_Lofl - 0.4)/m_Lofl;
		}
		else if (lat >= 1.4*m_Lofl)
		{
			Fa = 9*SQR(lat/m_Lofl - 1.6);
			FaDl = 18*(lat/m_Lofl - 1.6)/m_Lofl;
		}
		else if ((lat >= 0.6*m_Lofl) && (lat <= 1.4*m_Lofl))
		{
			Fa = 1 - 4*SQR(1 - lat/m_Lofl);
			FaDl = 8*(1 - lat/m_Lofl)/m_Lofl;
		}
	}


	// activation level
	double alpha = m_ascl*(m_plc ? m_plc->Value():1);

	// calculate total fiber force
	double FfDl = m_smax*(Fp + alpha*Fa)/m_Lofl;
	double FfD4  = 0.5*FfDl/lat;

	double FfDll = m_smax*(FpDl + alpha*FaDl)/m_Lofl;
	double FfD44 = 0.25*(FfDll - FfDl / lat)/I4;

	// add all derivatives
	W1 = F2D1;
	W2 = 0;
	W4 = F1D4 + F2D4 + F3D4 + FfD4;
	W5 = F1D5 + F2D5;

	// calculate second derivatives
	double W11, W12, W22, W14, W24, W15, W25, W44, W45, W55;

	W11 = F2D11;
	W12 = 0;
	W22 = 0;
	W14 = F2D14;
	W24 = 0;
	W15 = F2D15;
	W25 = 0;
	W44 = F1D44 + F2D44 + F3D44 + FfD44;
	W45 = F1D45 + F2D45;
	W55 = F2D55;

	// calculate dWdC:C
	double WCC = W1*I1 + 2*W2*I2 + W4*I4 + 2*W5*I5;

	// calculate C:d2WdCdC:C
	double CW2CCC = (W11*I1 + W12*I1*I1 + W2*I1 + 2*W12*I2 + 2*W22*I1*I2 + W14*I4 + W24*I1*I4 + 2*W15*I5 + 2*W25*I1*I5)*I1
		           -(W12*I1 + 2*W22*I2 + W2 + W24*I4 + 2*W25*I5)*(I1*I1 - 2*I2)
				   +(W14*I1 + 2*W24*I2 + W44*I4 + 2*W45*I5)*I4 + (W15*I1 + 2*W25*I2 + W45*I4 + 2*W55*I5)*2*I5
				   + 2*W5*I5;

	// calculate spatial tangent
	double A[3] = {a.x, a.y, a.z};
	double BA[3] = {Ba.x, Ba.y, Ba.z};

	// fourth-order identity IxI
	double IxI[6][6] = {
		{ 1, 1, 1, 0, 0, 0 },
		{ 1, 1, 1, 0, 0, 0 },
		{ 1, 1, 1, 0, 0, 0 },
		{ 0, 0, 0, 0, 0, 0 },
		{ 0, 0, 0, 0, 0, 0 },
		{ 0, 0, 0, 0, 0, 0 }};

	// second order identity tensor
	double ID[3][3] = {
		{ 1, 0, 0 },
		{ 0, 1, 0 },
		{ 0, 0, 1 }};

	// fourth-order identity I
	double  I[6][6] = {
		{ 1, 0, 0, 0, 0, 0 },
		{ 0, 1, 0, 0, 0, 0 },
		{ 0, 0, 1, 0, 0, 0 },
		{ 0, 0, 0, 0.5, 0, 0 },
		{ 0, 0, 0, 0, 0.5, 0 },
		{ 0, 0, 0, 0, 0, 0.5 }};

	// dev(S)x1
	double SxI[6][6] = {
		{ s[0][0], s[0][0], s[0][0], 0, 0, 0 },
		{ s[1][1], s[1][1], s[1][1], 0, 0, 0 },
		{ s[2][2], s[2][2], s[2][2], 0, 0, 0 },
		{ s[0][1], s[0][1], s[0][1], 0, 0, 0 },
		{ s[1][2], s[1][2], s[1][2], 0, 0, 0 },
		{ s[0][2], s[0][2], s[0][2], 0, 0, 0 }};

	//  1xdev(S)
	double IxS[6][6] = {
		{ s[0][0], s[1][1], s[2][2], s[0][1], s[1][2], s[0][2] },
		{ s[0][0], s[1][1], s[2][2], s[0][1], s[1][2], s[0][2] },
		{ s[0][0], s[1][1], s[2][2], s[0][1], s[1][2], s[0][2] },
		{ 0, 0, 0, 0, 0, 0 },
		{ 0, 0, 0, 0, 0, 0 },
		{ 0, 0, 0, 0, 0, 0 }};


	// --- calculate elasticity tensor ---
	int i, j, k, l, m, n;
	for (i=0; i<6; ++i)
		for (j=i; j<6; ++j)
		{
			D[i][j] = p*(IxI[i][j] - 2.0*I[i][j]) - (2.0/3.0)*(SxI[i][j] + IxS[i][j]) + (4.0/3.0*Ji)*WCC*(I[i][j] - IxI[i][j]/3.0); 
		}

	// now we add the cw part
	// first we need to define some additional quantities

	// calculate push-forward of dI5/dC
	double I5C[3][3];
	I5C[0][0] = I4*(A[0]*BA[0] + BA[0]*A[0]);
	I5C[0][1] = I4*(A[0]*BA[1] + BA[0]*A[1]);
	I5C[0][2] = I4*(A[0]*BA[2] + BA[0]*A[2]);

	I5C[1][0] = I4*(A[1]*BA[0] + BA[1]*A[0]);
	I5C[1][1] = I4*(A[1]*BA[1] + BA[1]*A[1]);
	I5C[1][2] = I4*(A[1]*BA[2] + BA[1]*A[2]);

	I5C[2][0] = I4*(A[2]*BA[0] + BA[2]*A[0]);
	I5C[2][1] = I4*(A[2]*BA[1] + BA[2]*A[1]);
	I5C[2][2] = I4*(A[2]*BA[2] + BA[2]*A[2]);

	// calculate push-forward of d2I5/dCdC
	double I5CC[6][6];

	I5CC[0][0] = 0.5*I4*(4*A[0]*A[0]*B[0][0]);	// (0,0) = [0,0,0,0]
	I5CC[1][1] = 0.5*I4*(4*A[1]*A[1]*B[1][1]);	// (1,1) = [1,1,1,1]
	I5CC[2][2] = 0.5*I4*(4*A[2]*A[2]*B[2][2]);	// (2,2) = [2,2,2,2]

	I5CC[0][1] = 0.5*I4*(4*A[0]*A[1]*B[0][1]);	// (0,1) = [0,0,1,1]
	I5CC[0][2] = 0.5*I4*(4*A[0]*A[2]*B[0][2]);	// (0,2) = [0,0,2,2]
	I5CC[1][2] = 0.5*I4*(4*A[1]*A[2]*B[1][2]);	// (1,2) = [1,1,2,2]

	I5CC[3][3] = 0.5*I4*(2*A[0]*A[1]*B[1][0] + A[0]*A[0]*B[1][1] + A[1]*A[1]*B[0][0]);	// (3,3) = [0,1,0,1]
	I5CC[4][4] = 0.5*I4*(2*A[1]*A[2]*B[2][1] + A[1]*A[1]*B[2][2] + A[2]*A[2]*B[1][1]);	// (4,4) = [1,2,1,2]
	I5CC[5][5] = 0.5*I4*(2*A[0]*A[2]*B[2][0] + A[0]*A[0]*B[2][2] + A[2]*A[2]*B[0][0]);	// (5,5) = [0,2,0,2]

	I5CC[0][3] = 0.5*I4*(A[0]*A[1]*B[0][0] + A[0]*A[0]*B[0][1] + A[0]*A[0]*B[0][1] + A[0]*A[1]*B[0][0]); // (0,3) = [0,0,0,1];
	I5CC[0][4] = 0.5*I4*(A[0]*A[2]*B[0][1] + A[0]*A[1]*B[0][2] + A[0]*A[1]*B[0][2] + A[0]*A[2]*B[0][1]); // (0,4) = [0,0,1,2];
	I5CC[0][5] = 0.5*I4*(A[0]*A[2]*B[0][0] + A[0]*A[0]*B[0][2] + A[0]*A[0]*B[0][2] + A[0]*A[2]*B[0][0]); // (0,5) = [0,0,0,2];
	I5CC[1][3] = 0.5*I4*(A[1]*A[1]*B[1][0] + A[1]*A[0]*B[1][1] + A[1]*A[0]*B[1][1] + A[1]*A[1]*B[1][0]); // (1,3) = [1,1,0,1];
	I5CC[1][4] = 0.5*I4*(A[1]*A[2]*B[1][1] + A[1]*A[1]*B[1][2] + A[1]*A[1]*B[1][2] + A[1]*A[2]*B[1][1]); // (1,4) = [1,1,1,2];
	I5CC[1][5] = 0.5*I4*(A[1]*A[2]*B[1][0] + A[1]*A[0]*B[1][2] + A[1]*A[0]*B[1][2] + A[1]*A[2]*B[1][0]); // (1,5) = [1,1,0,2];
	I5CC[2][3] = 0.5*I4*(A[2]*A[1]*B[2][0] + A[2]*A[0]*B[2][1] + A[2]*A[0]*B[2][1] + A[2]*A[1]*B[2][0]); // (2,3) = [2,2,0,1];
	I5CC[2][4] = 0.5*I4*(A[2]*A[2]*B[2][1] + A[2]*A[1]*B[2][2] + A[2]*A[1]*B[2][2] + A[2]*A[2]*B[2][1]); // (2,4) = [2,2,1,2];
	I5CC[2][5] = 0.5*I4*(A[2]*A[2]*B[2][0] + A[2]*A[0]*B[2][2] + A[2]*A[0]*B[2][2] + A[2]*A[2]*B[2][0]); // (2,5) = [2,2,0,2];

	I5CC[3][4] = 0.5*I4*(A[0]*A[2]*B[1][1] + A[0]*A[1]*B[1][2] + A[1]*A[1]*B[0][2] + A[1]*A[2]*B[0][1]); // (3,4) = [0,1,1,2];
	I5CC[3][5] = 0.5*I4*(A[0]*A[2]*B[1][0] + A[0]*A[0]*B[1][2] + A[1]*A[0]*B[0][2] + A[1]*A[2]*B[0][0]); // (3,5) = [0,1,0,2];
	I5CC[4][5] = 0.5*I4*(A[1]*A[2]*B[2][0] + A[1]*A[0]*B[2][2] + A[2]*A[0]*B[1][2] + A[2]*A[2]*B[1][0]); // (4,5) = [1,2,0,2];

	// calculate push forward of I
	double Ib[6][6];

	Ib[0][0] = 0.5*(B[0][0]*B[0][0] + B[0][0]*B[0][0]); // (0,0) = [0,0,0,0];
	Ib[1][1] = 0.5*(B[1][1]*B[1][1] + B[1][1]*B[1][1]); // (1,1) = [1,1,1,1];
	Ib[2][2] = 0.5*(B[2][2]*B[2][2] + B[2][2]*B[2][2]); // (2,2) = [2,2,2,2];

	Ib[0][1] = 0.5*(B[0][1]*B[0][1] + B[0][1]*B[0][1]); // (0,1) = [0,0,1,1];
	Ib[0][2] = 0.5*(B[0][2]*B[0][2] + B[0][2]*B[0][2]); // (0,2) = [0,0,2,2];
	Ib[1][2] = 0.5*(B[1][2]*B[1][2] + B[1][2]*B[1][2]); // (1,2) = [1,1,2,2];

	Ib[3][3] = 0.5*(B[0][0]*B[1][1] + B[0][1]*B[1][0]); // (3,3) = [0,1,0,1];
	Ib[4][4] = 0.5*(B[1][1]*B[2][2] + B[1][2]*B[2][1]); // (4,4) = [1,2,1,2];
	Ib[5][5] = 0.5*(B[0][0]*B[2][2] + B[0][2]*B[2][0]); // (5,5) = [0,2,0,2];

	Ib[0][3] = 0.5*(B[0][0]*B[0][1] + B[0][1]*B[0][0]); // (0,3) = [0,0,0,1];
	Ib[0][4] = 0.5*(B[0][1]*B[0][2] + B[0][2]*B[0][1]); // (0,4) = [0,0,1,2];
	Ib[0][5] = 0.5*(B[0][0]*B[0][2] + B[0][2]*B[0][0]); // (0,5) = [0,0,0,2];
	Ib[1][3] = 0.5*(B[1][0]*B[1][1] + B[1][1]*B[1][0]); // (1,3) = [1,1,0,1];
	Ib[1][4] = 0.5*(B[1][1]*B[1][2] + B[1][2]*B[1][1]); // (1,4) = [1,1,1,2];
	Ib[1][5] = 0.5*(B[1][0]*B[1][2] + B[1][2]*B[1][0]); // (1,5) = [1,1,0,2];
	Ib[2][3] = 0.5*(B[2][0]*B[2][1] + B[2][1]*B[2][0]); // (2,3) = [2,2,0,1];
	Ib[2][4] = 0.5*(B[2][1]*B[2][2] + B[2][2]*B[2][1]); // (2,4) = [2,2,1,2];
	Ib[2][5] = 0.5*(B[2][0]*B[2][2] + B[2][2]*B[2][0]); // (2,5) = [2,2,0,2];

	Ib[3][4] = 0.5*(B[0][1]*B[1][2] + B[0][2]*B[1][1]); // (3,4) = [0,1,1,2];
	Ib[3][5] = 0.5*(B[0][0]*B[1][2] + B[0][2]*B[1][0]); // (3,5) = [0,1,0,2];
	Ib[4][5] = 0.5*(B[1][0]*B[2][2] + B[1][2]*B[2][0]); // (4,5) = [1,2,0,2];

	// calculate push forward of dW/dCdC:C
	double WCCC[3][3];
	for (i=0; i<3; ++i)
		for (j=0; j<3; ++j)
		{
			WCCC[i][j]  = (W11*I1 + W12*I1*I1 + W2*I1 + 2*W12*I2 + 2*W22*I1*I2 + W14*I4 + W24*I1*I4 + 2*W15*I5 + 2*W25*I1*I5)*B[i][j];
			WCCC[i][j] -= (W12*I1 + 2*W22*I2 + W2 + W24*I4 + 2*W25*I5)*B2[i][j];
			WCCC[i][j] += (W14*I1 + 2*W24*I2 + W44*I4 + 2*W45*I5)*I4*A[i]*A[j];
			WCCC[i][j] += (W15*I1 + 2*W25*I2 + W45*I4 + 2*W55*I5 + W5)*I5C[i][j];
		}

	// this look-up-table is used to convert 4th-order tensors to the Voigt notation
	const int LUT[6][2] = {{0,0}, {1,1}, {2,2}, {0,1}, {1,2}, {0,2}};

	// calculate push-forward of dW2/dCdC
	double W2CC[6][6];
	for (m=0; m<6; ++m)
	{
		i = LUT[m][0];
		j = LUT[m][1];
		for (n=m; n<6; ++n)
		{
			k = LUT[n][0];
			l = LUT[n][1];

			W2CC[m][n] = (W11 + 2.0*W12*I1 + W2 + W22*I1*I1)*B[i][j]*B[k][l];
			W2CC[m][n] += -(W12+W22*I1)*(B[i][j]*B2[k][l] + B2[i][j]*B[k][l]) + W22*B2[i][j]*B2[k][l]- W2*Ib[m][n];
			W2CC[m][n] += (W14 + W24*I1)*I4*(B[i][j]*A[k]*A[l] + A[i]*A[j]*B[k][l]);
			W2CC[m][n] += (W15 + W25*I1)*(B[i][j]*I5C[k][l] + I5C[i][j]*B[k][l]);
			W2CC[m][n] += -W24*I4*(B2[i][j]*A[k]*A[l] + A[i]*A[j]*B2[k][l]);
			W2CC[m][n] += W44*I4*I4*A[i]*A[j]*A[k]*A[l];
			W2CC[m][n] += W45*I4*(A[i]*A[j]*I5C[k][l] + I5C[i][j]*A[k]*A[l]);
			W2CC[m][n] += W55*I5C[i][j]*I5C[k][l] + W5*I5CC[m][n];
		}
	}

	// let's put it all together
	for (m=0; m<6; ++m)
	{
		i = LUT[m][0];
		j = LUT[m][1];
		for (n=m; n<6; ++n)
		{
			k = LUT[n][0];
			l = LUT[n][1];

			D[m][n] += 4*Ji*W2CC[m][n] - (4.0/3.0*Ji)*(WCCC[i][j]*ID[k][l] + ID[i][j]*WCCC[k][l]) + (4.0/9.0*Ji)*(CW2CCC)*IxI[m][n];
		}
	}

	// set symmetric components
	D[1][0] = D[0][1]; D[2][0] = D[0][2]; D[3][0] = D[0][3]; D[4][0] = D[0][4]; D[5][0] = D[0][5];
	D[2][1] = D[1][2]; D[3][1] = D[1][3]; D[4][1] = D[1][4]; D[5][1] = D[1][5];
	D[3][2] = D[2][3]; D[4][2] = D[2][4]; D[5][2] = D[2][5];
	D[4][3] = D[3][4]; D[5][3] = D[3][5];
	D[5][4] = D[4][5];
}
