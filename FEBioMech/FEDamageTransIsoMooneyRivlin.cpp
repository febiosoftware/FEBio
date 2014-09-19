#include "stdafx.h"
#include "FEDamageTransIsoMooneyRivlin.h"

// define the material parameters
BEGIN_PARAMETER_LIST(FEDamageTransIsoMooneyRivlin, FEUncoupledMaterial)
	ADD_PARAMETER(m_c1, FE_PARAM_DOUBLE, "c1");
	ADD_PARAMETER(m_c2, FE_PARAM_DOUBLE, "c2");
	ADD_PARAMETER(m_c3, FE_PARAM_DOUBLE, "c3");
	ADD_PARAMETER(m_c4, FE_PARAM_DOUBLE, "c4");
	ADD_PARAMETER(m_Mbeta, FE_PARAM_DOUBLE, "Mbeta");
	ADD_PARAMETER(m_Msmin, FE_PARAM_DOUBLE, "Msmin");
	ADD_PARAMETER(m_Msmax, FE_PARAM_DOUBLE, "Msmax");
	ADD_PARAMETER(m_Fbeta, FE_PARAM_DOUBLE, "Fbeta");
	ADD_PARAMETER(m_Fsmin, FE_PARAM_DOUBLE, "Fsmin");
	ADD_PARAMETER(m_Fsmax, FE_PARAM_DOUBLE, "Fsmax");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FEDamageTransIsoMooneyRivlin::FEDamageTransIsoMooneyRivlin(FEModel* pfem) : FEUncoupledMaterial(pfem)
{

}

//-----------------------------------------------------------------------------
void FEDamageTransIsoMooneyRivlin::Init()
{
	if (m_c4 <= 0) throw MaterialError("c4 must be > 0");

}

//-----------------------------------------------------------------------------
mat3ds FEDamageTransIsoMooneyRivlin::DevStress(FEMaterialPoint& mp)
{
	// matrix stress
	mat3ds sm = MatrixStress(mp);

	// matrix reduction factor
	double gm = MatrixDamage(mp);

	// fiber stress
	mat3ds sf = FiberStress(mp);

	// fiber reduction factor
	double gf = FiberDamage(mp);

	return sm*gm + sf*gf;
}

//-----------------------------------------------------------------------------
//! Calculate the deviatoric stress
mat3ds FEDamageTransIsoMooneyRivlin::MatrixStress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// determinant of deformation gradient
	double J = pt.m_J;

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
	// W = C1*(I1 - 3) + C2*(I2 - 3)
	//
	// Wi = dW/dIi
	double W1 = m_c1;
	double W2 = m_c2;
	// ---

	// calculate T = F*dW/dC*Ft
	// T = F*dW/dC*Ft
	mat3ds T = B*(W1 + W2*I1) - B2*W2;

	return T.dev()*(2.0/J);
}

//-----------------------------------------------------------------------------
// Calculate the fiber stress
mat3ds FEDamageTransIsoMooneyRivlin::FiberStress(FEMaterialPoint &mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// deformation gradient
	mat3d& F = pt.m_F;
	double J = pt.m_J;
	double Jm13 = pow(J, -1.0/3.0);

	// get the initial fiber direction
	vec3d a0;
	a0.x = pt.m_Q[0][0];
	a0.y = pt.m_Q[1][0];
	a0.z = pt.m_Q[2][0];

	// calculate the current material axis lam*a = F*a0;
	vec3d a = F*a0;

	// normalize material axis and store fiber stretch
	double lam, lamd;
	lam = a.unit();
	lamd = lam*Jm13; // i.e. lambda tilde

	// invariant I4
	double I4 = lamd*lamd;

	// strain energy derivative
	double W4 = (I4 - 1)*m_c3*exp(m_c4*(I4-1)*(I4-1));

	// calculate dyad of a: AxA = (a x a)
	mat3ds AxA = dyad(a);

	// calculate FdWf/dCFt = I4*W4*(a x a)
	mat3ds T = AxA*(W4*I4);

	// return stress
	return T.dev()*(2.0/J);
}

//-----------------------------------------------------------------------------
tens4ds FEDamageTransIsoMooneyRivlin::DevTangent(FEMaterialPoint& mp)
{
	// matrix tangent
	tens4ds Cm = MatrixTangent(mp);

	// matrix damage
	double gm = MatrixDamage(mp);

	// fiber tangent
	tens4ds Cf = FiberTangent(mp);

	// fiber damage
	double gf = FiberDamage(mp);

	return Cm*gm + Cf*gf;
}

//-----------------------------------------------------------------------------
//! Calculate the deviatoric tangent
tens4ds FEDamageTransIsoMooneyRivlin::MatrixTangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// determinant of deformation gradient
	double J = pt.m_J;
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
	double W1, W2;
	W1 = m_c1;
	W2 = m_c2;
	// ---

	// calculate dWdC:C
	double WC = W1*I1 + 2*W2*I2;

	// calculate C:d2WdCdC:C
	double CWWC = 2*I2*W2;

	// deviatoric cauchy-stress, trs = trace[s]/3
	mat3ds devs = pt.m_s.dev();

	// Identity tensor
	mat3ds I(1,1,1,0,0,0);

	tens4ds IxI = dyad1s(I);
	tens4ds I4  = dyad4s(I);
	tens4ds BxB = dyad1s(B);
	tens4ds B4  = dyad4s(B);

	// d2W/dCdC:C
	mat3ds WCCxC = B*(W2*I1) - B2*W2;

	tens4ds cw = (BxB - B4)*(W2*4.0*Ji) - dyad1s(WCCxC, I)*(4.0/3.0*Ji) + IxI*(4.0/9.0*Ji*CWWC);

	tens4ds c = dyad1s(devs, I)*(-2.0/3.0) + (I4 - IxI/3.0)*(4.0/3.0*Ji*WC) + cw;

	return c;
}


//-----------------------------------------------------------------------------
// Fiber material tangent
//
tens4ds FEDamageTransIsoMooneyRivlin::FiberTangent(FEMaterialPoint &mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// deformation gradient
	mat3d& F = pt.m_F;
	double J = pt.m_J;
	double Jm13 = pow(J, -1.0/3.0);
	double Jm23 = Jm13*Jm13;
	double Ji = 1.0/J;

	// get initial local material axis
	vec3d a0;
	a0.x = pt.m_Q[0][0];
	a0.y = pt.m_Q[1][0];
	a0.z = pt.m_Q[2][0];

	// calculate current local material axis
	vec3d a = F*a0;

	double lam = a.unit();

	// deviatoric stretch
	double lamd = lam*Jm13;

	double I4 = lamd*lamd;

	// strain energy derivative
	double W4 = (I4 - 1)*m_c3*exp(m_c4*(I4-1)*(I4-1));
	double W44 = m_c3*exp(m_c4*(I4 - 1)*(I4 - 1)) + 2*m_c3*m_c4*(I4 - 1)*(I4 - 1)*exp(m_c4*(I4 - 1)*(I4 - 1));

	// calculate dWdC:C
	double WC = W4*I4;

	// calculate C:d2WdCdC:C
	double CWWC = W44*I4*I4;

	mat3dd I(1);	// Identity
	tens4ds IxI = dyad1s(I);
	tens4ds Id4  = dyad4s(I);

	mat3ds AxA = dyad(a);
	tens4ds AxAxAxA = dyad1s(AxA);

	tens4ds cw = AxAxAxA*(4.0*Ji*W44*I4*I4) - dyad1s(I, AxA)*(4.0/3.0*Ji*W44*I4*I4);

	tens4ds c = (Id4 - IxI/3.0)*(4.0/3.0*Ji*WC) + IxI*(4.0/9.0*Ji*CWWC) + cw;

	// see if we need to add the stress
	FETIMRDamageMaterialPoint& dp = *mp.ExtractData<FETIMRDamageMaterialPoint>();
	if (dp.m_FEtrial > dp.m_FEmax)
	{
		mat3ds devs = pt.m_s.dev();
		double dg = FiberDamageDerive(mp);
		c += dyad1s(devs)*(J*dg/dp.m_FEtrial);
	}

	return c;
}

//-----------------------------------------------------------------------------
double FEDamageTransIsoMooneyRivlin::DevStrainEnergyDensity(FEMaterialPoint& mp)
{
	// matrix sed
	double sedm = MatrixStrainEnergyDensity(mp);
    
	// matrix reduction factor
	double gm = MatrixDamage(mp);
    
	// fiber sed
	double sedf = FiberStrainEnergyDensity(mp);
    
	// fiber reduction factor
	double gf = FiberDamage(mp);
    
	return sedm*gm + sedf*gf;
}

//-----------------------------------------------------------------------------
//! Calculate the matrix strain energy density
double FEDamageTransIsoMooneyRivlin::MatrixStrainEnergyDensity(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
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
	// W = C1*(I1 - 3) + C2*(I2 - 3)
	//
    double sed = m_c1*(I1 - 3) + m_c2*(I2 - 3);
    
	return sed;
}

//-----------------------------------------------------------------------------
// Calculate the fiber strain energy density
double FEDamageTransIsoMooneyRivlin::FiberStrainEnergyDensity(FEMaterialPoint &mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
	// deformation gradient
	mat3d& F = pt.m_F;
	double J = pt.m_J;
	double Jm13 = pow(J, -1.0/3.0);
    
	// get the initial fiber direction
	vec3d a0;
	a0.x = pt.m_Q[0][0];
	a0.y = pt.m_Q[1][0];
	a0.z = pt.m_Q[2][0];
    
	// calculate the current material axis lam*a = F*a0;
	vec3d a = F*a0;
    
	// normalize material axis and store fiber stretch
	double lam, lamd;
	lam = a.unit();
	lamd = lam*Jm13; // i.e. lambda tilde
    
	// invariant I4
	double I4 = lamd*lamd;
    
	// strain energy derivative
	double sed = 0.5*m_c3/m_c4*(exp(m_c4*(I4-1)*(I4-1))-1);
    
	// return sed
	return sed;
}

//-----------------------------------------------------------------------------
// Calculate damage reduction factor for matrix
double FEDamageTransIsoMooneyRivlin::MatrixDamage(FEMaterialPoint &mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// calculate right Cauchy-Green tensor
	mat3ds C = pt.DevRightCauchyGreen();
	mat3ds C2 = C*C;

	// Invariants
	double I1 = C.tr();
	double I2 = 0.5*(I1*I1 - C2.tr());

	// strain-energy value
	double SEF = m_c1*(I1 - 3) + m_c2*(I2 - 3);

	// get the damage material point data
	FETIMRDamageMaterialPoint& dp = *mp.ExtractData<FETIMRDamageMaterialPoint>();

	// calculate trial-damage parameter
	dp.m_MEtrial = sqrt(2.0*fabs(SEF));

	// calculate damage parameter
	double Es = max(dp.m_MEtrial, dp.m_MEmax);

	// calculate reduction parameter
	double g = 1.0;
	if (Es < m_Msmin) g = 1.0;
	else if (Es > m_Msmax) g = 0.0;
	else 
	{
		double F = (Es - m_Msmin)/(m_Msmin - m_Msmax);
		g = 1.0 - (1.0 - m_Mbeta + m_Mbeta*F*F)*(F*F);
	}

	dp.m_Dm = 1-g;
	return g;
}

//-----------------------------------------------------------------------------
// Calculate damage reduction factor for matrix
double FEDamageTransIsoMooneyRivlin::MatrixDamageDerive(FEMaterialPoint &mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// calculate right Cauchy-Green tensor
	mat3ds C = pt.DevRightCauchyGreen();
	mat3ds C2 = C*C;

	// Invariants
	double I1 = C.tr();
	double I2 = 0.5*(I1*I1 - C2.tr());

	// strain-energy value
	double SEF = m_c1*(I1 - 3) + m_c2*(I2 - 3);

	// get the damage material point data
	FETIMRDamageMaterialPoint& dp = *mp.ExtractData<FETIMRDamageMaterialPoint>();

	// calculate trial-damage parameter
	dp.m_MEtrial = sqrt(2.0*fabs(SEF));

	// calculate damage parameter
	double Es = max(dp.m_MEtrial, dp.m_MEmax);

	// calculate reduction parameter
	double dg = 0.0;
	if (Es < m_Msmin) dg = 0.0;
	else if (Es > m_Msmax) dg = 0.0;
	else 
	{
		double h = m_Msmax - m_Msmin;
		double F = (Es - m_Msmin)/h;
		dg = -2.0*F/h*(1 - m_Mbeta + m_Mbeta*F*F)-(F*F)*(2.0*m_Mbeta*F/h);
	}

	return dg;
}


//-----------------------------------------------------------------------------
// Calculate damage reduction factor for fibers
double FEDamageTransIsoMooneyRivlin::FiberDamage(FEMaterialPoint &mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// deformation gradient
	mat3d& F = pt.m_F;
	double J = pt.m_J;
	double Jm13 = pow(J, -1.0/3.0);

	// get the initial fiber direction
	vec3d a0;
	a0.x = pt.m_Q[0][0];
	a0.y = pt.m_Q[1][0];
	a0.z = pt.m_Q[2][0];

	// calculate the current material axis lam*a = F*a0;
	vec3d a = F*a0;

	// normalize material axis and store fiber stretch
	double lam, lamd;
	lam = a.unit();
	lamd = lam*Jm13; // i.e. lambda tilde

	// invariant I4
	double I4 = lamd*lamd;

	// strain energy value
	double SEF = 0.5*m_c3/m_c4*(exp(m_c4*(I4-1)*(I4-1))-1);

	// get the damage material point data
	FETIMRDamageMaterialPoint& dp = *mp.ExtractData<FETIMRDamageMaterialPoint>();

	// calculate trial-damage parameter
	dp.m_FEtrial = sqrt(2.0*fabs(SEF));

	// calculate damage parameter
	double Es = max(dp.m_FEtrial, dp.m_FEmax);

	// calculate reduction parameter
	double g = 1.0;
	if (Es < m_Fsmin) g = 1.0;
	else if (Es > m_Fsmax) g = 0.0;
	else 
	{
		double F = (Es - m_Fsmin)/(m_Fsmin - m_Fsmax);
		g = 1.0 - (1.0 - m_Fbeta + m_Fbeta*F*F)*(F*F);
	}

	dp.m_Df = 1-g;
	return g;
}


//-----------------------------------------------------------------------------
// Calculate damage reduction factor for fibers
double FEDamageTransIsoMooneyRivlin::FiberDamageDerive(FEMaterialPoint &mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// deformation gradient
	mat3d& F = pt.m_F;
	double J = pt.m_J;
	double Jm13 = pow(J, -1.0/3.0);

	// get the initial fiber direction
	vec3d a0;
	a0.x = pt.m_Q[0][0];
	a0.y = pt.m_Q[1][0];
	a0.z = pt.m_Q[2][0];

	// calculate the current material axis lam*a = F*a0;
	vec3d a = F*a0;

	// normalize material axis and store fiber stretch
	double lam, lamd;
	lam = a.unit();
	lamd = lam*Jm13; // i.e. lambda tilde

	// invariant I4
	double I4 = lamd*lamd;

	// strain energy value
	double SEF = 0.5*m_c3/m_c4*(exp(m_c4*(I4-1)*(I4-1))-1);

	// get the damage material point data
	FETIMRDamageMaterialPoint& dp = *mp.ExtractData<FETIMRDamageMaterialPoint>();

	// calculate trial-damage parameter
	dp.m_FEtrial = sqrt(2.0*fabs(SEF));

	// calculate damage parameter
	double Es = max(dp.m_FEtrial, dp.m_FEmax);

	// calculate reduction parameter
	double dg = 0.0;
	if (Es < m_Fsmin) dg = 0.0;
	else if (Es > m_Fsmax) dg = 0.0;
	else 
	{
		double h = m_Fsmin - m_Fsmax;
		double F = (Es - m_Fsmin)/h;
		dg = -2.0*F/h*(1 - m_Fbeta + m_Fbeta*F*F)-(F*F)*(2.0*m_Fbeta*F/h);
	}

	return dg;
}
