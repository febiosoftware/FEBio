#include "stdafx.h"
#include "FEDamageNeoHookean.h"

// register the material with the framework
REGISTER_MATERIAL(FEDamageNeoHookean, "damage neo-Hookean");

// define the material parameters
BEGIN_PARAMETER_LIST(FEDamageNeoHookean, FEElasticMaterial)
	ADD_PARAMETER(m_E, FE_PARAM_DOUBLE, "E");
	ADD_PARAMETER(m_v, FE_PARAM_DOUBLE, "v");
	ADD_PARAMETER(m_alpha, FE_PARAM_DOUBLE, "a");
	ADD_PARAMETER(m_beta , FE_PARAM_DOUBLE, "b");
END_PARAMETER_LIST();


//-----------------------------------------------------------------------------
// Constructor
FEDamageNeoHookean::FEDamageNeoHookean(void)
{
	m_E = 0;
	m_v = 0;

	m_alpha = 0.014;
	m_beta = 0.34;
}

//-----------------------------------------------------------------------------
// Initialization routine and parameter checking
void FEDamageNeoHookean::Init()
{
	if (m_E <= 0) throw MaterialError("Invalid value for E");
	if (!INRANGE(m_v, -1.0, 0.5)) throw MaterialError("Invalid value for v");
	if (!INRANGE(m_beta, 0.0, 1.0)) throw MaterialError("Invalid value for b: must be in range [0,1]");
	if (m_alpha < 0) throw MaterialError("Invalid value of a: must be a non-negative number");

	// calculate Lame parameters
	m_lam = m_v*m_E/((1+m_v)*(1-2*m_v));
	m_mu  = 0.5*m_E/(1+m_v);
}

//-----------------------------------------------------------------------------
//! Calculate the stress. This happens in two phases. First, we calculate 
//! the stress for the undamaged material. Second, we update the damage
//! parameter and correct the stress accordingly.
mat3ds FEDamageNeoHookean::Stress(FEMaterialPoint& mp)
{
	// --- A. Calculate neo-Hookean stress ----
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	mat3d &F = pt.F;
	double detF = pt.J;
	double detFi = 1.0/detF;
	double lndetF = log(detF);

	// calculate left Cauchy-Green tensor
	mat3ds b = pt.LeftCauchyGreen();

	// Identity
	mat3dd I(1);

	// calculate stress
	mat3ds s = (b - I)*(m_mu*detFi) + I*(m_lam*lndetF*detFi);

	// --- B. Calculate the damage reduction factor ---
	double g = Damage(mp);

	return s*g;
}

//-----------------------------------------------------------------------------
// Calculate damage reduction factor 
double FEDamageNeoHookean::Damage(FEMaterialPoint &mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// get the deformation gradient
	mat3d &F = pt.F;

	// calculate right Cauchy-Green tensor
	mat3ds C = pt.RightCauchyGreen();

	// Invariants
	double I1 = C.tr();
	double J = pt.J;

	// strain-energy value
	double lnJ = log(J);
	double SEF = 0.5*m_mu*(I1 - 3) - m_mu*lnJ + 0.5*m_lam*(lnJ*lnJ);

	// get the damage material point data
	FEDamageMaterialPoint& dp = *mp.ExtractData<FEDamageMaterialPoint>();

	// calculate trial-damage parameter
	dp.m_Etrial = sqrt(2.0*fabs(SEF));

	// calculate damage parameter
	double Es = max(dp.m_Etrial, dp.m_Emax);

	// calculate reduction parameter
	double g = 1.0;
	if (fabs(Es) > 1e-12) g = m_beta + (1.0 - m_beta)*(1.0 - exp(-Es/m_alpha))/(Es/m_alpha);
	else g = 1.0 - 0.5*(1.0 - m_beta)/m_alpha*Es;

	dp.m_D = g;
	return g;
}

//-----------------------------------------------------------------------------
// Calculate tangent. I'm not sure if the tangent needs to be modified for the 
// damage model For now, I don't modify it.
tens4ds FEDamageNeoHookean::Tangent(FEMaterialPoint& mp)
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
	
	double D[6][6] = {0};
	D[0][0] = lam1+2.*mu1; D[0][1] = lam1       ; D[0][2] = lam1       ;
	D[1][0] = lam1       ; D[1][1] = lam1+2.*mu1; D[1][2] = lam1       ;
	D[2][0] = lam1       ; D[2][1] = lam1       ; D[2][2] = lam1+2.*mu1;
	D[3][3] = mu1;
	D[4][4] = mu1;
	D[5][5] = mu1;

	double g = Damage(mp);

	return tens4ds(D)*g;
}
