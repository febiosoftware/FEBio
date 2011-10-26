#include "FEVonMisesPlasticity.h"

//-----------------------------------------------------------------------------
// register the material with the framework
REGISTER_MATERIAL(FEVonMisesPlasticity, "von-Mises plastic");

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_PARAMETER_LIST(FEVonMisesPlasticity, FEElasticMaterial)
	ADD_PARAMETER(m_K, FE_PARAM_DOUBLE, "K");
	ADD_PARAMETER(m_G, FE_PARAM_DOUBLE, "G");
	ADD_PARAMETER(m_Y, FE_PARAM_DOUBLE, "Y");
END_PARAMETER_LIST();


//-----------------------------------------------------------------------------
FEVonMisesPlasticity::FEVonMisesPlasticity(void)
{
	m_K = m_G = m_Y = 0;
}

//-----------------------------------------------------------------------------
void FEVonMisesPlasticity::Init()
{
	if (m_K <= 0) throw MaterialError("K must be postive number");
	if (m_G <= 0) throw MaterialError("G must be postive number");
	if (m_Y <= 0) throw MaterialError("Y must be postitive number");
}

//-----------------------------------------------------------------------------
mat3ds FEVonMisesPlasticity::Stress(FEMaterialPoint &mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	FEJ2PlasticMaterialPoint& pp = *mp.ExtractData<FEJ2PlasticMaterialPoint>();
	mat3d& F = pt.F;
	
	// get the current strain
	mat3ds e = F.sym() - mat3dd(1.0);

	// calculate strain increment
	mat3ds de = (e - pp.e0);

	// get the trial stress
	mat3ds strial = pp.sn + (de.dev()*(2.0*m_G) + de.iso()*(3.0*m_K));

	double k = m_Y / sqrt(3.0);
	double fac = strial.dev().norm() / (sqrt(2.0)*k);

	mat3ds s;
	if (fac<=1)
	{
		s = strial;
		pp.b = false;
	}
	else
	{
		s = strial.iso() + strial.dev() / fac;
		pp.b = true;
	}

	// store the current strain measure
	pp.e1 = e;

	return s;
}

//-----------------------------------------------------------------------------
tens4ds FEVonMisesPlasticity::Tangent(FEMaterialPoint &mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	FEJ2PlasticMaterialPoint& pp = *mp.ExtractData<FEJ2PlasticMaterialPoint>();

	// lame parameters
	double lam = m_K - m_G*2.0/3.0;
	double mu  = m_G;

	double D[6][6] = {0};
	D[0][0] = lam+2.*mu; D[0][1] = lam      ; D[0][2] = lam      ;
	D[1][0] = lam      ; D[1][1] = lam+2.*mu; D[1][2] = lam      ;
	D[2][0] = lam      ; D[2][1] = lam      ; D[2][2] = lam+2.*mu;
	D[3][3] = mu;
	D[4][4] = mu;
	D[5][5] = mu;
	tens4ds C(D);

	// see if we are in plastic flow mode
	if (pp.b)
	{
		// get the stress
		mat3ds s = pt.s;
		mat3ds n = s.dev()*2.0;

		mat3ds A = C.dot(n);
		double G = n.dotdot(C.dot(n));

		C -= dyad4s(A)/G;
	}

	return C;
}
