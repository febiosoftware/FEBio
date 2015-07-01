#include "FEThermalConductivity.h"
#include <FEBioMech/FEElasticMaterial.h>
#include "FEHeatTransferMaterial.h"

//-----------------------------------------------------------------------------
BEGIN_PARAMETER_LIST(FEConstReferenceThermalConductivity, FEThermalConductivity);
	ADD_PARAMETER(m_k0, FE_PARAM_DOUBLE, "k0");
	ADD_PARAMETER(m_wt, FE_PARAM_DOUBLE, "wt");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FEConstReferenceThermalConductivity::FEConstReferenceThermalConductivity(FEModel* pfem) : FEThermalConductivity(pfem)
{
	m_k0 = 0;
	m_wt = 0;
}

//-----------------------------------------------------------------------------
void FEConstReferenceThermalConductivity::Init()
{
	if (m_k0 <= 0.0) throw MaterialError("k0 must be positive.");
	if (m_wt <  0.0) throw MaterialError("wt must be non-negative.");
}

//-----------------------------------------------------------------------------
// evaluate conductivity
mat3ds FEConstReferenceThermalConductivity::Conductivity(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());
	FEHeatMaterialPoint& ht = *(mp.ExtractData<FEHeatMaterialPoint>());

	mat3ds b = pt.RightCauchyGreen();
	double J = pt.m_J;
	double T = ht.m_T;
	double T0 = ht.m_T0;

	double k = m_k0*(1 - m_wt*(T - T0))/J;

	return b*k;
}

//-----------------------------------------------------------------------------
// evaluate the derivative of the conductivity with respect to C
tens4ds FEConstReferenceThermalConductivity::Tangent_Conductivity_Strain(FEMaterialPoint& mp)
{
	tens4ds T(0.0);
	return T;
}

//-----------------------------------------------------------------------------
// evaluate the derivative of the conductivity with respect to temperature
mat3ds FEConstReferenceThermalConductivity::Tangent_Conductivity_Temperature(FEMaterialPoint& mp)
{
	// TODO: Implement this
	mat3ds C(0.0);
	return C;
}

//-----------------------------------------------------------------------------
BEGIN_PARAMETER_LIST(FEConstThermalConductivity, FEThermalConductivity);
	ADD_PARAMETER(m_k0, FE_PARAM_DOUBLE, "k0");
	ADD_PARAMETER(m_wt, FE_PARAM_DOUBLE, "wt");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FEConstThermalConductivity::FEConstThermalConductivity(FEModel* pfem) : FEThermalConductivity(pfem)
{
	m_k0 = 0;
	m_wt = 0;
}

//-----------------------------------------------------------------------------
void FEConstThermalConductivity::Init()
{
	if (m_k0 <= 0.0) throw MaterialError("k0 must be positive.");
	if (m_wt <  0.0) throw MaterialError("wt must be non-negative.");
}

//-----------------------------------------------------------------------------
// evaluate conductivity
mat3ds FEConstThermalConductivity::Conductivity(FEMaterialPoint& mp)
{
	FEHeatMaterialPoint& ht = *(mp.ExtractData<FEHeatMaterialPoint>());

	// evaluate the temperature-dependant conductivity
	double T = ht.m_T;
	double T0 = ht.m_T0;
	double k = m_k0*(1 - m_wt*(T - T0));

	mat3dd I(1.0);
	return I*k;
}

//-----------------------------------------------------------------------------
// evaluate the derivative of the conductivity with respect to C
tens4ds FEConstThermalConductivity::Tangent_Conductivity_Strain(FEMaterialPoint& mp)
{
	FEHeatMaterialPoint& ht = *(mp.ExtractData<FEHeatMaterialPoint>());

	// evaluate the temperature-dependant conductivity
	double T = ht.m_T;
	double T0 = ht.m_T0;
	double k = m_k0*(1 - m_wt*(T - T0));

	mat3dd I(1.0);
	tens4ds IxI = dyad1s(I);
	tens4ds I4 = dyad4s(I);

	tens4ds D = (IxI - I4*2.0)*k;
	return D;
}

//-----------------------------------------------------------------------------
// evaluate the derivative of the conductivity with respect to temperature
mat3ds FEConstThermalConductivity::Tangent_Conductivity_Temperature(FEMaterialPoint& mp)
{
	// TODO: Implement this
	mat3ds C(0.0);
	return C;
}
