#include "stdafx.h"
#include "FEActiveFiberContraction.h"
#include "FEElasticMaterial.h"

//-----------------------------------------------------------------------------
BEGIN_PARAMETER_LIST(FEActiveFiberContraction, FEMaterial);
	ADD_PARAMETER(m_ascl , "ascl");
	ADD_PARAMETER(m_Tmax , "Tmax");
	ADD_PARAMETER(m_ca0  , "ca0");
	ADD_PARAMETER(m_camax, "camax");
	ADD_PARAMETER(m_beta , "beta");
	ADD_PARAMETER(m_l0   , "l0");
	ADD_PARAMETER(m_refl , "refl");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FEActiveFiberContraction::FEActiveFiberContraction(FEModel* pfem) : FEMaterial(pfem)
{
	m_ascl = 0;
	m_Tmax = 1.0;
	m_ca0 = 1.0;
	m_camax = 0.0;
}

//-----------------------------------------------------------------------------
bool FEActiveFiberContraction::Init()
{
	if (FEMaterial::Init() == false) return false;

	// for backward compatibility we set m_camax to m_ca0 if it is not defined
	if (m_camax == 0.0) m_camax = m_ca0;
	if (m_camax <= 0.0) return MaterialError("camax must be larger than zero");

	return true;
}

//-----------------------------------------------------------------------------
mat3ds FEActiveFiberContraction::FiberStress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// get the deformation gradient
	mat3d F = pt.m_F;
	double J = pt.m_J;
	double Jm13 = pow(J, -1.0 / 3.0);

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

	// calculate dyad of a: AxA = (a x a)
	mat3ds AxA = dyad(a);

	// get the activation
	double saf = 0.0;
	if (m_ascl > 0)
	{
		double ctenslm = m_ascl;

		// current sarcomere length
		double strl = m_refl*lamd;

		// sarcomere length change
		double dl = strl - m_l0;

		if (dl >= 0)
		{
			// calcium sensitivity
			double eca50i = (exp(m_beta*dl) - 1);

			// ratio of Camax/Ca0
			double rca = m_camax / m_ca0;

			// active fiber stress
			saf = m_Tmax*(eca50i / (eca50i + rca*rca))*ctenslm;
		}
	}
	return AxA*saf;
}

//-----------------------------------------------------------------------------
tens4ds FEActiveFiberContraction::FiberStiffness(FEMaterialPoint& mp)
{
	/*	if (lcna >= 0)
	{
	double ctenslm = m_plc->Value();

	// current sarcomere length
	double strl = refl*lamd;

	// sarcomere length change
	double dl = strl - l0;

	if (dl >= 0) W44 += J*2*beta*refl*exp(-beta*dl);
	}
	*/
	return tens4ds(0.0);
}
