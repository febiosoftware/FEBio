// FEArrudaBoyce.cpp: implementation of the FEMArrudaBoyce class.
//
// After Kalliske & Rothert, Eng Computations 14(2)(1997):216-232
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FEArrudaBoyce.h"

// register the material with the framework
REGISTER_MATERIAL(FEArrudaBoyce, "Arruda-Boyce");

// define the material parameters
BEGIN_PARAMETER_LIST(FEArrudaBoyce, FEIncompressibleMaterial)
	ADD_PARAMETER(c1, FE_PARAM_DOUBLE, "c1");
	ADD_PARAMETER(c2, FE_PARAM_DOUBLE, "c2");
END_PARAMETER_LIST();

//////////////////////////////////////////////////////////////////////
// FEArrudaBoyce
//////////////////////////////////////////////////////////////////////

// rename material parameters
#define mu c1
#define N c2

void FEArrudaBoyce::Init()
{
	FEIncompressibleMaterial::Init();

	if (m_K <= 0.0) throw MaterialError("Invalid value for k");
}

mat3ds FEArrudaBoyce::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

        // Check the value for N is >0
        if (c2 <= 0.0) throw MaterialError("Invalid value for N");

	const double third = 1.0/3.0;

	// deformation gradient
	mat3d &F = pt.F;
	double J = pt.J;
	double Ji = 1.0/J;
	double Jm13 = pow(J, -1.0/3.0);
	double Jm23 = Jm13*Jm13;
	double twoJi = 2.0*Ji;

	// left Cauchy-Green tensor and its square
	double B[3][3];

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

	// Invariants of B (= invariants of C)
	// Note that these are the invariants of Btilde, not of B!
	double I1 = B[0][0]+B[1][1]+B[2][2];

	// T = a1.B; called tau_tilde by Kalliske
	double T[3][3];
	double IoN=I1/N;
	double IoN2=IoN*IoN;
	double a1 = 2.0*mu*(0.5+0.1*IoN+11.0*IoN2/350.0+19.0*IoN2*IoN/1750.0+519.0*IoN2*IoN2/134750.0);

	T[0][0] = a1*B[0][0];
	T[0][1] = a1*B[0][1];
	T[0][2] = a1*B[0][2];

	T[1][1] = a1*B[1][1];
	T[1][2] = a1*B[1][2];

	T[2][2] = a1*B[2][2];

	// trT = tr(T)/3
	double trT = (T[0][0] + T[1][1] + T[2][2])*third;

	// get pressure at material point
	double p = pt.avgp;

	// calculate stress: s = pI + (1/J)dev[T]
	mat3ds s;

	s.xx() = p + Ji*(T[0][0] - trT);
	s.yy() = p + Ji*(T[1][1] - trT);
	s.zz() = p + Ji*(T[2][2] - trT);
	s.xy() = Ji*T[0][1];
	s.yz() = Ji*T[1][2];
	s.xz() = Ji*T[0][2];

	return s;
}

tens4ds FEArrudaBoyce::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	const double third = 1.0 / 3.0;
        const double fn=4.0/9.0;

	// deformation gradient
	mat3d &F = pt.F;
	double J = pt.J;
	double Jm13 = pow(J, -1.0/3.0);
	double Jm23 = Jm13*Jm13;
	double Ji = 1.0/J;

	// calculate deviatoric left Cauchy-Green tensor: B = F*Ft
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

      	double I1 = B[0][0]+B[1][1]+B[2][2];  // first invariant of B

        // calculate tau_tilde (called T here)
	double T[3][3];
        double IoN=I1/N;
        double IoN2=IoN*IoN;
	double a1 = 2.0*mu*(0.5+0.1*IoN+11.0*IoN2/350.0+19.0*IoN2*IoN/1750.0+519.0*IoN2*IoN2/134750.0);

	T[0][0] = a1*B[0][0];
	T[0][1] = a1*B[0][1];
	T[0][2] = a1*B[0][2];

	T[1][0] = a1*B[1][0];
	T[1][1] = a1*B[1][1];
	T[1][2] = a1*B[1][2];

	T[2][0] = a1*B[2][0];
	T[2][1] = a1*B[2][1];
	T[2][2] = a1*B[2][2];

	// trT = tr(T)/3
	double trT = (T[0][0] + T[1][1] + T[2][2])*third;

        // Only the deviatoric part of T is needed from now on...
        T[0][0]-=trT;
        T[1][1]-=trT;
        T[2][2]-=trT;

	// calculate a_tilde (called A here)
	// only the three diagonal components need to be found now
	double b1 = 4.0*mu*(0.1/N+22.0*IoN/(350.0*N)+57.0*IoN2/(1750.0*N)+2076.0*IoN*IoN2/(134750.0*N));

	double At0000=b1*B[0][0]*B[0][0];
	double At1111=b1*B[1][1]*B[1][1];
	double At2222=b1*B[2][2]*B[2][2];
	double trA=third*(At0000+At1111+At2222); // volumetric part, to be subtracted

	// calculate the constitutive tensor ...

	double D[6][6] = {0};

	// D[0][0] = c(0,0,0,0)
	D[0][0] = Ji*(fn*trT-4.0*third*T[0][0]+At0000-trA);

	// D[1][1] = c(1,1,1,1)
	D[1][1] = Ji*(fn*trT-4.0*third*T[1][1]+At1111-trA);

	// D[2][2] = c(2,2,2,2)
	D[2][2] = Ji*(fn*trT-4.0*third*T[2][2]+At2222-trA);


	// D[0][1] = D[1][0] = c(0,0,1,1)
	D[0][1] = Ji*(-0.5*fn*trT-2.0*third*(T[0][0]+T[1][1]));
        D[0][1] += Ji*b1*B[0][0]*B[1][1];

	// D[1][2] = D[2][1] = c(1,1,2,2)
	D[1][2] = Ji*(-0.5*fn*trT-2.0*third*(T[1][1]+T[2][2]));
        D[1][2] += Ji*b1*B[1][1]*B[2][2];

	// D[0][2] = D[2][0] = c(0,0,2,2)
	D[0][2] = Ji*(-0.5*fn*trT-2.0*third*(T[0][0]+T[2][2]));
        D[0][2] += Ji*b1*B[0][0]*B[2][2];


	// D[3][3] = 0.5*(c(0,1,0,1) + c(0,1,1,0))
        D[3][3] = Ji*third*trT;
        D[3][3] += Ji*0.5*b1*(B[0][1]*B[0][1]+B[0][1]*B[1][0]);

	// D[4][4] = 0.5*(c(1,2,1,2) + c(1,2,2,1))
        D[4][4]=Ji*third*trT;
        D[4][4] += Ji*0.5*b1*(B[1][2]*B[1][2]+B[1][2]*B[2][1]);

	// D[5][5] = 0.5*(c(0,2,0,2) + c(0,2,2,0))
        D[5][5] = Ji*third*trT;
        D[5][5] += Ji*0.5*b1*(B[0][2]*B[0][2]+B[0][2]*B[2][0]);


	// D[0][3] = 0.5*(c(0,0,0,1) + c(0,0,1,0))
	D[0][3] = -Ji*third*(T[0][1]+T[1][0]);
        D[0][3] += Ji*0.5*b1*(B[0][0]*B[0][1]+B[0][0]*B[1][0]);

	// D[0][4] = 0.5*(c(0,0,1,2) + c(0,0,2,1))
	D[0][4] = -Ji*third*(T[1][2]+T[2][1]);
        D[0][4] += Ji*0.5*b1*(B[0][0]*B[1][2]+B[0][0]*B[2][1]);

	// D[0][5] = 0.5*(c(0,0,0,2) + c(0,0,2,0))
	D[0][5] = -Ji*third*(T[0][2]+T[2][0]);
        D[0][5] += Ji*0.5*b1*(B[0][0]*B[0][2]+B[0][0]*B[2][0]);

	// D[1][3] = 0.5*(c(1,1,0,1) + c(1,1,1,0))
	D[1][3] = -Ji*third*(T[0][1]+T[1][0]);
        D[1][3] += Ji*0.5*b1*(B[1][1]*B[0][1]+B[1][1]*B[1][0]);

	// D[1][4] = 0.5*(c(1,1,1,2) + c(1,1,2,1))
	D[1][4] = -Ji*third*(T[1][2]+T[2][1]);
        D[1][4] += Ji*0.5*b1*(B[1][1]*B[1][2]+B[1][1]*B[2][1]);

	// D[1][5] = 0.5*(c(1,1,0,2) + c(1,1,2,0))
	D[1][5] = -Ji*third*(T[0][2]+T[2][0]);
        D[1][5] += Ji*0.5*b1*(B[1][1]*B[0][2]+B[1][1]*B[2][0]);

	// D[2][3] = 0.5*(c(2,2,0,1) + c(2,2,1,0))
	D[2][3] = -Ji*third*(T[0][1]+T[1][0]);
        D[2][3] += Ji*0.5*b1*(B[2][2]*B[0][1]+B[2][2]*B[1][0]);

	// D[2][4] = 0.5*(c(2,2,1,2) + c(2,2,2,1))
	D[2][4] = -Ji*third*(T[1][2]+T[2][1]);
        D[2][4] += Ji*0.5*b1*(B[2][2]*B[1][2]+B[2][2]*B[2][1]);

	// D[2][5] = 0.5*(c(2,2,0,2) + c(2,2,2,0))
	D[2][5] = -Ji*third*(T[0][2]+T[2][0]);
        D[2][5] += Ji*0.5*b1*(B[2][2]*B[0][2]+B[2][2]*B[2][0]);


	// D[3][4] = 0.5*(c(0,1,1,2) + c(0,1,2,1))
        D[3][4] = Ji*0.5*b1*(B[0][1]*B[1][2]+B[0][1]*B[2][1]);

	// D[3][5] = 0.5*(c(0,1,0,2) + c(0,1,2,0))
        D[3][5] = Ji*0.5*b1*(B[0][1]*B[0][2]+B[0][1]*B[2][0]);

	// D[4][5] = 0.5*(c(1,2,0,2) + c(1,2,2,0))
        D[4][5] = Ji*0.5*b1*(B[1][2]*B[0][2]+B[1][2]*B[2][0]);


	// set symmetric components
	D[1][0] = D[0][1]; D[2][0] = D[0][2]; D[3][0] = D[0][3]; D[4][0] = D[0][4]; D[5][0] = D[0][5];
	D[2][1] = D[1][2]; D[3][1] = D[1][3]; D[4][1] = D[1][4]; D[5][1] = D[1][5];
	D[3][2] = D[2][3]; D[4][2] = D[2][4]; D[5][2] = D[2][5];
	D[4][3] = D[3][4]; D[5][3] = D[3][5];
	D[5][4] = D[4][5];

	return tens4ds(D);
}
