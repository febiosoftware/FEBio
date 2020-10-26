/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2020 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/
#include "stdafx.h"
#include "FEContinuousElasticDamage.h"
#include <FECore/FEModel.h>

// Macauley Bracket
#define MB(a) ((a)>0.0?(a):0.0)
#define SIGN(a) ((a)>0.0?1.0:0.0)

//=========================================================================================

class FEFiberDamagePoint : public FEMaterialPoint
{
public:
	FEFiberDamagePoint(FEMaterialPoint* pm) : FEMaterialPoint(pm)
	{
		m_D = 0.0;
		m_psi_f0_ini = 0.0;
		m_psi_f0 = 0.0;
		m_psi_f0_prev = 0.0;
		m_bt_ini = 0.0;
		m_bt = 0.0;
		m_bt_prev = 0.0;
		m_gamma = 0.0;
		m_gamma_prev = 0.0;
		m_psf_c = 0.0;
		m_init = false;
	}

	void Init() override
	{
		m_D = 0.0;
		m_psi_f0_ini = m_psf_c;// we set this to the initial value
		m_psi_f0 = m_psf_c;		// we set this to the initial value
		m_psi_f0_prev = m_psf_c; // we set this to the initial value
		m_bt_ini = 0.0;
		m_bt = 0.0;
		m_bt_prev = 0.0;
		m_gamma = 0.0;
		m_gamma_prev = 0.0;
		m_init = false;
	}

	void Update(const FETimeInfo& timeInfo) override
	{
		m_gamma_prev = m_gamma;
		m_bt_prev = m_bt;
		m_psi_f0_prev = m_psi_f0;
	}

public:
	bool	m_init;	// initialization flag
	double	m_D;		// accumulated damage

	double	m_psi_f0_ini, m_psf_c;
	double	m_psi_f0, m_psi_f0_prev;

	double	m_bt_ini;
	double	m_bt, m_bt_prev;

	double	m_gamma, m_gamma_prev;
};

//=========================================================================================
double FEDamageInterface::Damage(FEMaterialPoint& mp)
{
	FEFiberDamagePoint& damagePoint = *mp.ExtractData<FEFiberDamagePoint>();
	return damagePoint.m_D;
}

//=========================================================================================

BEGIN_FECORE_CLASS(FEContinuousElasticDamage, FEElasticMaterial)
	ADD_PARAMETER(m_c1, FE_RANGE_GREATER(0.0), "c1");
	ADD_PARAMETER(m_k, FE_RANGE_GREATER(0.0), "k");
	ADD_PARAMETER(m_fiber, "fiber");
	ADD_PARAMETER(m_fib.m_a1, FE_RANGE_GREATER_OR_EQUAL(0.0), "a1");
	ADD_PARAMETER(m_fib.m_a2, FE_RANGE_GREATER(1.0), "a2");
	ADD_PARAMETER(m_fib.m_kappa, FE_RANGE_CLOSED(0.0, 2.0/3.0), "kappa");
	ADD_PARAMETER(m_fib.m_tinit, FE_RANGE_GREATER_OR_EQUAL(0.0), "t0");
	ADD_PARAMETER(m_fib.m_Dmax, FE_RANGE_CLOSED(0.0, 1.0), "Dmax");
	ADD_PARAMETER(m_fib.m_beta_s, FE_RANGE_GREATER(0.0), "beta_s");
	ADD_PARAMETER(m_fib.m_gamma_max, FE_RANGE_GREATER(0.0), "gamma_max");
END_FECORE_CLASS();

FEContinuousElasticDamage::FEContinuousElasticDamage(FEModel* fem) : FEElasticMaterial(fem), m_fib(fem)
{
	m_c1 = 0.0;
	m_k = 0.0;
	m_fiber = vec3d(1, 0, 0);
}

FEMaterialPoint* FEContinuousElasticDamage::CreateMaterialPointData()
{
	return new FEFiberDamagePoint(new FEElasticMaterialPoint);
}

double FEContinuousElasticDamage::StrainEnergyDensity(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	FEFiberDamagePoint& damagePoint = *mp.ExtractData<FEFiberDamagePoint>();

	mat3ds C = pt.RightCauchyGreen();
	mat3ds C2 = C.sqr();
	double I1 = C.tr();
	double J = pt.m_J;
	double J23 = pow(J, 2. / 3.);
	double lnJ = log(J);

	double I1_tilde = I1 / J23;

	// matrix strain-energy
	double Psi_iso = m_c1 * (I1_tilde - 3.0);
	double Psi_vol = 0.5*m_k*(lnJ*lnJ);
	double Psi_matrix = Psi_iso + Psi_vol;

	// fiber strain energy
	vec3d a0 = m_fiber(mp);
	double Psi_fiber = m_fib.FiberStrainEnergyDensity(mp, a0);

	// add it all up
	double Psi = Psi_matrix + Psi_fiber;

	return Psi;
}

mat3ds FEContinuousElasticDamage::Stress(FEMaterialPoint& mp)
{
	// calculate the matrix stress
	mat3ds s_matrix = MatrixStress(mp);

	// calculate the fiber stress
	vec3d a0 = m_fiber(mp);
	mat3ds s_fiber = m_fib.FiberStress(mp, a0);

	// add it all up
	mat3ds s = s_matrix + s_fiber;

	// all done
	return s;
}

mat3ds FEContinuousElasticDamage::MatrixStress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	double J = pt.m_J;
	double Jm53 = pow(J, -5. / 3.);

	mat3ds b = pt.LeftCauchyGreen();
	double I1 = b.tr();

	double p = m_k * log(J) / J;

	mat3dd I(1.0);
	mat3ds s_iso = (b - I*(I1/3.0))*(2.0*m_c1 * Jm53);
	mat3ds s_vol = I * p;

	mat3ds s = s_iso + s_vol;
	
	return s;
}

tens4ds FEContinuousElasticDamage::Tangent(FEMaterialPoint& mp)
{
	// calculate matrix tangent
	tens4ds c_matrix = MatrixTangent(mp);

	// calculate fiber tangent
	vec3d a0 = m_fiber(mp);
	tens4ds c_fiber = m_fib.FiberTangent(mp, a0);

	// add it all up
	tens4ds c = c_matrix + c_fiber;

	// all done
	return c;
}

//-----------------------------------------------------------------------------
//! calculate spatial tangent stiffness at material point, using secant method
mat3ds FEContinuousElasticDamage::SecantStress(FEMaterialPoint& mp)
{
	// extract the deformation gradient
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	mat3d F = pt.m_F;
	double J = pt.m_J;
	mat3ds E = pt.Strain();
	mat3dd I(1);
	mat3d FiT = F.transinv();

	// calculate the 2nd P-K stress at the current deformation gradient
	double W = StrainEnergyDensity(mp);

	// create deformation gradient increment
	double eps = 1e-9;
	vec3d e[3];
	e[0] = vec3d(1, 0, 0); e[1] = vec3d(0, 1, 0); e[2] = vec3d(0, 0, 1);
	mat3ds S(0.0);
	for (int k = 0; k < 3; ++k) {
		for (int l = k; l < 3; ++l) {
			// evaluate incremental stress
			mat3d dF = FiT * ((e[k] & e[l]))*(eps*0.5);
			mat3d F1 = F + dF;
			pt.m_F = F1;
			pt.m_J = pt.m_F.det();

			double dW = StrainEnergyDensity(mp) - W;

			// evaluate the secant modulus
			S(k, l) = 2.0 * dW / eps;
		}
	}

	// push from material to spatial frame
	mat3ds s = pt.push_forward(S);

	// restore values
	pt.m_F = F;
	pt.m_J = J;

	// return secant stress
	return s;
}

tens4ds FEContinuousElasticDamage::MatrixTangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	double J = pt.m_J;
	double Jm53 = pow(J, -5.0 / 3.0);

	mat3ds b = pt.LeftCauchyGreen();
	double I1 = b.tr();

	mat3ds I(mat3dd(1.0));

	double p = m_k * log(J) / J;
	double dp = m_k * (1.0 - log(J))/(J*J);

	tens4ds I4 = dyad4s(I);
	tens4ds IxI = dyad1s(I);

	tens4ds c_iso = (I4*I1 - dyad1s(I, b) + IxI*(I1 / 3.0))*(4.0*m_c1*Jm53 / 3.0);

	tens4ds c_vol = IxI * (dp*J) + (IxI - I4*2.0)*p;

	tens4ds c = c_iso + c_vol;

	return c;
}

//=================================================================================================
BEGIN_FECORE_CLASS(FEDamageFiberPower, FEElasticFiberMaterial)
	ADD_PARAMETER(m_a1, FE_RANGE_GREATER_OR_EQUAL(0.0), "a1");
	ADD_PARAMETER(m_a2, FE_RANGE_GREATER(1.0), "a2");
	ADD_PARAMETER(m_kappa, FE_RANGE_CLOSED(0.0, 2.0/3.0), "kappa");
	ADD_PARAMETER(m_tinit, FE_RANGE_GREATER_OR_EQUAL(0.0), "t0");
	ADD_PARAMETER(m_Dmax, FE_RANGE_CLOSED(0.0, 1.0), "Dmax");
	ADD_PARAMETER(m_beta_s, FE_RANGE_GREATER(0.0), "beta_s");
	ADD_PARAMETER(m_gamma_max, FE_RANGE_GREATER(0.0), "gamma_max");
END_FECORE_CLASS();

FEDamageFiberPower::FEDamageFiberPower(FEModel* fem) : FEElasticFiberMaterial(fem)
{
	m_a1 = 0.0;
	m_a2 = 0.0;
	m_kappa = 0.0;

	m_tinit = 1e9;	// large value so, no damage accumulation by default
	m_Dmax = 1.0;

	m_beta_s = 0.0;
	m_gamma_max = 0.0;

	// Looks like these are hard-coded
	m_r_s = 0.99;
	m_r_inf = 0.99;
}

FEMaterialPoint* FEDamageFiberPower::CreateMaterialPointData()
{
	FEFiberDamagePoint* mp = new FEFiberDamagePoint(new FEElasticMaterialPoint);
	mp->m_psf_c = 2.0;
	return mp;
}

//! Strain energy density
double FEDamageFiberPower::FiberStrainEnergyDensity(FEMaterialPoint& mp, const vec3d& a0)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	FEFiberDamagePoint& damagePoint = *mp.ExtractData<FEFiberDamagePoint>();

	mat3ds C = pt.RightCauchyGreen();
	mat3ds C2 = C.sqr();
	double I1 = C.tr();
	double J = pt.m_J;
	double J23 = pow(J, 2. / 3.);
	double lnJ = log(J);

	double I1_tilde = I1 / J23;

	// fiber strain energy
	double I4 = a0 * (C*a0);
	double I5 = a0 * (C2*a0);
	double K3 = I1 * I4 - I5;

	double D = damagePoint.m_D;
	double Psi0 = m_kappa * I1 + (1.0 - 3.0*m_kappa / 2.0)*K3;
	double P = (1.0 - D)*Psi0 - 2.0;
	double Psi_fiber = m_a1 * pow(P, m_a2);

	return Psi_fiber;
}

// calculate stress in fiber direction a0
mat3ds FEDamageFiberPower::FiberStress(FEMaterialPoint& mp, const vec3d& a0)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	FEFiberDamagePoint& damagePoint = *mp.ExtractData<FEFiberDamagePoint>();

	mat3d F = pt.m_F;
	mat3ds C = pt.RightCauchyGreen();
	mat3ds C2 = C.sqr();
	mat3ds b = pt.LeftCauchyGreen();
	double I1 = C.tr();
	double J = pt.m_J;

	vec3d a = F * a0;
	a.unit();

	double I4 = a0 * (C *a0);
	double I5 = a0 * (C2*a0);
	double K3 = I1 * I4 - I5;

	double g1 = 1.0 - 3.0*m_kappa / 2.0;

	// get internal variables
	double D = damagePoint.m_D;
	double bt_prev = damagePoint.m_bt_prev;
	double psi_f0_prev = damagePoint.m_psi_f0_prev;
	double gamma_prev = damagePoint.m_gamma_prev;

	// get current simulation time.
	double t = GetFEModel()->GetTime().currentTime;

	// (i) compute trans-iso strain energy
	double psi_f0 = m_kappa * I1 + g1 * K3;

	// (ii) check initial damage state
	double eps = 1e-9; // NOTE: should be machine eps
	if (t >= (m_tinit - eps))
	{
		// (b) compute bt
		double bt = bt_prev + MB(psi_f0 - psi_f0_prev);

		// init damage 
		if (damagePoint.m_init == false)
		{
			damagePoint.m_psi_f0_ini = psi_f0;
			damagePoint.m_bt_ini = bt;
			damagePoint.m_init = true;
		}

		// (iii) calculate max damage saturation value 
		// trial criterion
		double phi_trial = MB(psi_f0 - damagePoint.m_psi_f0_ini) - gamma_prev;

		// check algorithmic saturation criterion
		double gamma = 0;
		if (phi_trial > eps) gamma = MB(psi_f0 - damagePoint.m_psi_f0_ini);
		else gamma = gamma_prev;
		assert(gamma >= gamma_prev);

		// compute damage saturation value
		double Ds = m_Dmax * (1.0 - exp(log(1.0 - m_r_inf)*gamma / m_gamma_max));

		// (iv) compute internal variable
		double beta = MB(bt - damagePoint.m_bt_ini);

		// (v) evaluate damage function
		D = Ds * (1.0 - exp(log(1.0 - m_r_s)*beta / m_beta_s));

		// update internal variables
		damagePoint.m_bt = bt;
		damagePoint.m_psi_f0 = psi_f0;
		damagePoint.m_gamma = gamma;
		damagePoint.m_D = D;
	}

	double P = (1.0 - D)*(psi_f0)-2.0;
	if (P < 0.0) P = 0.0;
	double dm = m_a1 * m_a2*pow(P, m_a2 - 1.0);

	mat3dd I(1.0);
	mat3ds m = dyad(a);
	mat3ds aob = dyads(a, b*a);
	mat3ds ts = b * I4 + m * (I1*I4) - aob * I4;
	mat3ds s = (b*m_kappa + ts * g1)*(2 * dm*(1.0 - D) / J);

	return s;
}

// Spatial tangent
tens4ds FEDamageFiberPower::FiberTangent(FEMaterialPoint& mp, const vec3d& a0)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	FEFiberDamagePoint& damagePoint = *mp.ExtractData<FEFiberDamagePoint>();

	mat3d F = pt.m_F;
	mat3ds C = pt.RightCauchyGreen();
	mat3ds C2 = C.sqr();
	mat3ds b = pt.LeftCauchyGreen();
	double I1 = C.tr();
	double J = pt.m_J;

	vec3d a = F * a0;
	a.unit();

	double I4 = a0 * (C *a0);
	double I5 = a0 * (C2*a0);
	double K3 = I1 * I4 - I5;

	double k = 1.0 - 3.0*m_kappa / 2.0;

	double D = damagePoint.m_D;

	double Psi0 = m_kappa * I1 + k * K3;
	double P = (1.0 - D)*(Psi0)-2.0;
	if (P < 0) P = 0.0;
	double dm = m_a1 * m_a2*pow(P, m_a2 - 1.0);
	double d2m = m_a1 * m_a2 * (m_a2 - 1.0)*pow(P, m_a2 - 2.0);

	mat3ds m = dyad(a);
	mat3ds aBa = dyads(a, b*a);

	mat3ds t = b * I4 + m * (I1*I4) - aBa * I4;
	mat3ds s0 = (b * m_kappa + t * k)*(2.0 / J);

	tens4ds aba = dyad4s(a, b, a)*(I4*2.0);
	tens4ds bom = dyad1s(b, m);
	tens4ds bot = dyad1s(b, t);
	tens4ds txt = dyad1s(t);
	tens4ds bxb = dyad1s(b);
	tens4ds sxs = dyad1s(s0);

	// elastic stiffness 
	tens4ds c0 = (bom*I4 - aba)*(4.0*k / J);
	tens4ds ceD = sxs * (J*d2m*(1.0 - D)) + c0 * ((1.0 - D)*dm);

	tens4ds c = ceD;

	// damage stiffness
	double bt = damagePoint.m_bt;
	double gamma = damagePoint.m_gamma;
	double Ds = m_Dmax * (1.0 - exp(log(1.0 - m_r_inf)*gamma / m_gamma_max));
	double beta = MB(bt - damagePoint.m_bt_ini);

	double psi0_prev = damagePoint.m_psi_f0_prev;
	double psi0_ini = damagePoint.m_psi_f0_ini;
	double gamma_prev = damagePoint.m_gamma_prev;

	double dD_dbeta = -Ds * (log(1 - m_r_s) / m_beta_s)*exp(log(1 - m_r_s)*beta / m_beta_s);
	double dDs_dgamma = -m_Dmax * (log(1 - m_r_inf) / m_gamma_max)*exp(log(1 - m_r_inf)*gamma / m_gamma_max);
	double dD_dDs = 1.0 - exp(log(1 - m_r_s)*beta / m_beta_s);
	double dbeta_dpsi0 = 0.25*(SIGN(bt - damagePoint.m_bt_ini) + 1)*(SIGN(Psi0 - psi0_prev) + 1);
	double dgamma_dpsi0 = 0.5*(SIGN(Psi0 - psi0_ini) + 1);

	double meps = 1e-9; // NOTE: should be machine eps
	if (Psi0 - psi0_prev > meps)
	{
		tens4ds Cd = sxs * ((dm + d2m * (1 - D)*Psi0)*dD_dbeta*dbeta_dpsi0);
		c -= Cd;
	}

	double phi_trial = MB(Psi0 - damagePoint.m_psi_f0_ini) - gamma_prev;
	if (phi_trial > meps)
	{
		tens4ds Cd = sxs * ((dm + d2m * (1 - D)*Psi0)*dD_dDs*dDs_dgamma*dgamma_dpsi0);
		c -= Cd;
	}

	return c;
}


//=================================================================================================
BEGIN_FECORE_CLASS(FEDamageFiberExponential, FEElasticFiberMaterial)
	ADD_PARAMETER(m_k1, FE_RANGE_GREATER_OR_EQUAL(0.0), "k1");
	ADD_PARAMETER(m_k2, FE_RANGE_GREATER(1.0), "k2");
	ADD_PARAMETER(m_kappa, FE_RANGE_CLOSED(0.0, 1.0/3.0), "kappa");
	ADD_PARAMETER(m_tinit, FE_RANGE_GREATER_OR_EQUAL(0.0), "t0");
	ADD_PARAMETER(m_Dmax, FE_RANGE_CLOSED(0.0, 1.0), "Dmax");
	ADD_PARAMETER(m_beta_s, FE_RANGE_GREATER(0.0), "beta_s");
	ADD_PARAMETER(m_gamma_max, FE_RANGE_GREATER(0.0), "gamma_max");
END_FECORE_CLASS();

FEDamageFiberExponential::FEDamageFiberExponential(FEModel* fem) : FEElasticFiberMaterial(fem)
{
	m_k1 = 0.0;
	m_k2 = 0.0;
	m_kappa = 0.0;

	m_tinit = 1e9;	// large value so, no damage accumulation by default
	m_Dmax = 1.0;

	m_beta_s = 0.0;
	m_gamma_max = 0.0;

	// Looks like these are hard-coded
	m_r_s = 0.99;
	m_r_inf = 0.99;
}

FEMaterialPoint* FEDamageFiberExponential::CreateMaterialPointData()
{
	FEFiberDamagePoint* mp = new FEFiberDamagePoint(new FEElasticMaterialPoint);
	mp->m_psf_c = 1.0;
	return mp;
}

//! Strain energy density
double FEDamageFiberExponential::FiberStrainEnergyDensity(FEMaterialPoint& mp, const vec3d& a0)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	FEFiberDamagePoint& damagePoint = *mp.ExtractData<FEFiberDamagePoint>();

	mat3ds C = pt.RightCauchyGreen();
	mat3ds C2 = C.sqr();
	double I1 = C.tr();
	double J = pt.m_J;
	double J23 = pow(J, 2. / 3.);
	double lnJ = log(J);

	double I1_tilde = I1 / J23;

	// fiber strain energy
	double I4 = a0 * (C*a0);

	double D = damagePoint.m_D;
	double Psi0 = m_kappa * I1 + (1.0 - 3.0*m_kappa)*I4;
	double P = (1.0 - D)*Psi0 - 2.0;
	if (P < 0) P = 0.0;
	double Psi_fiber = 0.5*m_k1*(exp(m_k2*P*P) - 1) / m_k2;

	return Psi_fiber;
}

// calculate stress in fiber direction a0
mat3ds FEDamageFiberExponential::FiberStress(FEMaterialPoint& mp, const vec3d& a0)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	FEFiberDamagePoint& damagePoint = *mp.ExtractData<FEFiberDamagePoint>();

	mat3d F = pt.m_F;
	mat3ds C = pt.RightCauchyGreen();
	mat3ds C2 = C.sqr();
	mat3ds b = pt.LeftCauchyGreen();
	double I1 = C.tr();
	double J = pt.m_J;

	vec3d a = F * a0;
	a.unit();

	double I4 = a0 * (C *a0);

	double g1 = 1.0 - 3.0*m_kappa;

	// get internal variables
	double D = damagePoint.m_D;
	double bt_prev = damagePoint.m_bt_prev;
	double psi_f0_prev = damagePoint.m_psi_f0_prev;
	double gamma_prev = damagePoint.m_gamma_prev;

	// get current simulation time.
	double t = GetFEModel()->GetTime().currentTime;

	// (i) compute trans-iso strain energy
	double psi_f0 = m_kappa * I1 + g1 * I4;

	// (ii) check initial damage state
	double eps = 1e-9; // NOTE: should be machine eps
	if (t >= (m_tinit - eps))
	{
		// (b) compute bt
		double bt = bt_prev + MB(psi_f0 - psi_f0_prev);

		// init damage 
		if (damagePoint.m_init == false)
		{
			damagePoint.m_psi_f0_ini = psi_f0;
			damagePoint.m_bt_ini = bt;
			damagePoint.m_init = true;
		}

		// (iii) calculate max damage saturation value 
		// trial criterion
		double phi_trial = MB(psi_f0 - damagePoint.m_psi_f0_ini) - gamma_prev;

		// check algorithmic saturation criterion
		double gamma = 0;
		if (phi_trial > eps) gamma = MB(psi_f0 - damagePoint.m_psi_f0_ini);
		else gamma = gamma_prev;
		assert(gamma >= gamma_prev);

		// compute damage saturation value
		double Ds = m_Dmax * (1.0 - exp(log(1.0 - m_r_inf)*gamma / m_gamma_max));

		// (iv) compute internal variable
		double beta = MB(bt - damagePoint.m_bt_ini);

		// (v) evaluate damage function
		D = Ds * (1.0 - exp(log(1.0 - m_r_s)*beta / m_beta_s));

		// update internal variables
		damagePoint.m_bt = bt;
		damagePoint.m_psi_f0 = psi_f0;
		damagePoint.m_gamma = gamma;
		damagePoint.m_D = D;
	}

	double P = (1.0 - D)*(psi_f0)-1.0;
	if (P < 0.0) P = 0.0;
	double dm = m_k1 * P*exp(m_k2*P*P);

	mat3dd I(1.0);
	mat3ds m = dyad(a);
	mat3ds s = (b*m_kappa + m * (I4*g1))*(2 * dm*(1.0 - D) / J);

	return s;
}

// Spatial tangent
tens4ds FEDamageFiberExponential::FiberTangent(FEMaterialPoint& mp, const vec3d& a0)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	FEFiberDamagePoint& damagePoint = *mp.ExtractData<FEFiberDamagePoint>();

	mat3d F = pt.m_F;
	mat3ds C = pt.RightCauchyGreen();
	mat3ds C2 = C.sqr();
	mat3ds b = pt.LeftCauchyGreen();
	double I1 = C.tr();
	double J = pt.m_J;

	vec3d a = F * a0;
	a.unit();

	double I4 = a0 * (C *a0);

	double k = 1.0 - 3.0*m_kappa;

	double D = damagePoint.m_D;

	double Psi0 = m_kappa * I1 + k * I4;
	double P = (1.0 - D)*(Psi0)-2.0;
	if (P < 0) P = 0.0;
	double dm = m_k1 * P*exp(m_k2*P*P);
	double d2m = m_k1*(1.0 + 2.0*m_k2* P*P)* exp(m_k2*P*P);

	mat3ds m = dyad(a);
	mat3ds s0 = (b * m_kappa + m *(I4 * k))*(2.0 / J);

	tens4ds sxs = dyad1s(s0);

	// elastic stiffness 
	tens4ds c0(0.0);
	tens4ds ceD = sxs * (J*d2m*(1.0 - D)) + c0 * ((1.0 - D)*dm);

	tens4ds c = ceD;

	// damage stiffness
	double bt = damagePoint.m_bt;
	double gamma = damagePoint.m_gamma;
	double Ds = m_Dmax * (1.0 - exp(log(1.0 - m_r_inf)*gamma / m_gamma_max));
	double beta = MB(bt - damagePoint.m_bt_ini);

	double psi0_prev = damagePoint.m_psi_f0_prev;
	double psi0_ini = damagePoint.m_psi_f0_ini;
	double gamma_prev = damagePoint.m_gamma_prev;

	double dD_dbeta = -Ds * (log(1 - m_r_s) / m_beta_s)*exp(log(1 - m_r_s)*beta / m_beta_s);
	double dDs_dgamma = -m_Dmax * (log(1 - m_r_inf) / m_gamma_max)*exp(log(1 - m_r_inf)*gamma / m_gamma_max);
	double dD_dDs = 1.0 - exp(log(1 - m_r_s)*beta / m_beta_s);
	double dbeta_dpsi0 = 0.25*(SIGN(bt - damagePoint.m_bt_ini) + 1)*(SIGN(Psi0 - psi0_prev) + 1);
	double dgamma_dpsi0 = 0.5*(SIGN(Psi0 - psi0_ini) + 1);

	double meps = 1e-9; // NOTE: should be machine eps
	if (Psi0 - psi0_prev > meps)
	{
		tens4ds Cd = sxs * ((dm + d2m * (1 - D)*Psi0)*dD_dbeta*dbeta_dpsi0);
		c -= Cd;
	}

	double phi_trial = MB(Psi0 - damagePoint.m_psi_f0_ini) - gamma_prev;
	if (phi_trial > meps)
	{
		tens4ds Cd = sxs * ((dm + d2m * (1 - D)*Psi0)*dD_dDs*dDs_dgamma*dgamma_dpsi0);
		c -= Cd;
	}

	return c;
}
