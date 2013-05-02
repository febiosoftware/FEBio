#include "stdafx.h"
#include "FEUncoupledActiveContraction.h"

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_PARAMETER_LIST(FEUncoupledActiveContraction, FEUncoupledMaterial)
	ADD_PARAMETER(m_Tmax , FE_PARAM_DOUBLE, "Tmax" );
	ADD_PARAMETER(m_ca0  , FE_PARAM_DOUBLE, "ca0"  );
	ADD_PARAMETER(m_camax, FE_PARAM_DOUBLE, "camax");
	ADD_PARAMETER(m_beta , FE_PARAM_DOUBLE, "beta" );
	ADD_PARAMETER(m_l0   , FE_PARAM_DOUBLE, "l0"   );
	ADD_PARAMETER(m_refl , FE_PARAM_DOUBLE, "refl" );
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FEUncoupledActiveContraction::FEUncoupledActiveContraction()
{

}

//-----------------------------------------------------------------------------
void FEUncoupledActiveContraction::Init()
{

}

//-----------------------------------------------------------------------------
mat3ds FEUncoupledActiveContraction::DevStress(FEMaterialPoint &mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// deformation gradient
	mat3d& F = pt.F;
	double J = pt.J;
	double Jm13 = pow(J, -1.0/3.0);

	// get the initial fiber direction
	vec3d a0;
	a0.x = pt.Q[0][0];
	a0.y = pt.Q[1][0];
	a0.z = pt.Q[2][0];

	// calculate the current material axis lam*a = F*a0;
	vec3d a = F*a0;

	// normalize material axis and store fiber stretch
	double lam = a.unit();
	double lamd = lam*Jm13; // i.e. lambda tilde

	// current sarcomere length
	double strl = m_refl*lamd;

	// sarcomere length change
	double dl = strl - m_l0;

	// calculate stress
	mat3ds s(0.0);
	if (dl >= 0)
	{
		// calcium sensitivity
		double eca50i = (exp(m_beta*dl) - 1);

		// ratio of Camax/Ca0
		double rca = m_camax/m_ca0;

		// active fiber stress
		double saf = m_Tmax*(eca50i / ( eca50i + rca*rca ));

		// calculate dyad of a: AxA = (a x a)
		mat3ds AxA = dyad(a);

		// add saf*(a x a)
		s += AxA*saf;
	}

	return s;
}

//-----------------------------------------------------------------------------
tens4ds FEUncoupledActiveContraction::DevTangent(FEMaterialPoint &pt)
{
	// TODO: implement this
	return tens4ds(0.0);
}
