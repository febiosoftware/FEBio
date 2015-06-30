#include "FEThermoNeoHookean.h"
#include "FEHeatTransferMaterial.h"

//-----------------------------------------------------------------------------
BEGIN_PARAMETER_LIST(FEThermoNeoHookean, FEThermalElastic);
	ADD_PARAMETER(m_mu   , FE_PARAM_DOUBLE, "G"    );
	ADD_PARAMETER(m_k    , FE_PARAM_DOUBLE, "K"    );
	ADD_PARAMETER(m_a0   , FE_PARAM_DOUBLE, "a0"   );
	ADD_PARAMETER(m_gamma, FE_PARAM_DOUBLE, "gamma");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! constructor
FEThermoNeoHookean::FEThermoNeoHookean(FEModel* pfem) : FEThermalElastic(pfem)
{

}

//-----------------------------------------------------------------------------
void FEThermoNeoHookean::Init()
{
	if (m_mu <= 0.0) throw MaterialError("G must be positive.");
	if (m_k  <= 0.0) throw MaterialError("K must be positive.");
	if (m_a0 <  0.0) throw MaterialError("a0 must be non-negative.");
	if (m_gamma <= 0.0) throw MaterialError("gamma must be positive.");
}

//-----------------------------------------------------------------------------
//! Cauchy stress
mat3ds FEThermoNeoHookean::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());
	FEHeatMaterialPoint& ht = *(mp.ExtractData<FEHeatMaterialPoint>());

	// get the jacobian
	double J = pt.m_J;
	double Jg = pow(J, m_gamma - 1.0);

	// get the initial and current temperature
	double T = ht.m_T;
	double T0 = ht.m_T0;

	// calculate temperature-dependant properties
	double mu = m_mu*T/T0;
	double k  = m_k*T/T0;

	// evaluate Cauchy stress
	mat3dd I(1.0);
	mat3ds b = pt.LeftCauchyGreen();

	mat3ds s = (b - I)*(mu/J) + I*(k*log(J)/J) - I*(3.0*m_a0*m_k*Jg*(T - T0));

	return s;
}

//-----------------------------------------------------------------------------
//! spatial elasticity tangent
tens4ds FEThermoNeoHookean::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());
	FEHeatMaterialPoint& ht = *(mp.ExtractData<FEHeatMaterialPoint>());

	// get the jacobian
	double J = pt.m_J;
	double Jg = pow(J, m_gamma - 1.0);
	double lnJ = log(J);

	// get the initial and current temperature
	double T = ht.m_T;
	double T0 = ht.m_T0;

	// calculate temperature-dependant properties
	double mu = m_mu*T/T0;
	double k  = m_k*T/T0;

	// evaluate spatial tangent
	mat3dd I(1.0);
	tens4ds IxI = dyad1s(I);
	tens4ds I4  = dyad4s(I);

	tens4ds c = I4*(2.0*(mu - k*lnJ)/J) + IxI*(k/J) - (IxI - I4*2.0)*(3.0*m_a0*m_k*Jg*(T - T0));

	return c;
}

//-----------------------------------------------------------------------------
mat3ds FEThermoNeoHookean::ThermalTangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());
	FEHeatMaterialPoint& ht = *(mp.ExtractData<FEHeatMaterialPoint>());

	// get the jacobian
	double J = pt.m_J;
	double Jg = pow(J, m_gamma - 1.0);

	// get the initial temperature
	double T0 = ht.m_T0;

	// evaluate thermal tangent
	mat3dd I(1.0);
	mat3ds b = pt.LeftCauchyGreen();

	mat3ds dsdT = (b - I)*(m_mu/(T0*J)) + I*(m_k*log(J)/(J*T0)) - I*(3.0*m_a0*m_k*Jg);

	return dsdT;
}
