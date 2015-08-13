#include "FEThermoNeoHookean.h"
#include "FEHeatTransferMaterial.h"

//-----------------------------------------------------------------------------
BEGIN_PARAMETER_LIST(FEThermoNeoHookean, FEThermalElastic);
	ADD_PARAMETER2(m_mu   , FE_PARAM_DOUBLE, FE_RANGE_GREATER(0.0), "G"    );
	ADD_PARAMETER2(m_k    , FE_PARAM_DOUBLE, FE_RANGE_GREATER(0.0), "K"    );
	ADD_PARAMETER2(m_gamma, FE_PARAM_DOUBLE, FE_RANGE_GREATER(0.0), "gamma");
	ADD_PARAMETER2(m_a0   , FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "a0"   );
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! constructor
FEThermoNeoHookean::FEThermoNeoHookean(FEModel* pfem) : FEThermalElastic(pfem)
{

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

	tens4ds c = I4*(2.0*(mu - k*lnJ)/J) + IxI*(k/J) - (IxI*m_gamma - I4*2.0)*(3.0*m_a0*m_k*Jg*(T - T0));

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
