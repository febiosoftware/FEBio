#include "stdafx.h"
#include "tens4d.h"
#include "FEBulkElasticity.h"

// register the material with the framework
REGISTER_MATERIAL(FEBulkElasticity, "bulk elasticity");

// define the material parameters
BEGIN_PARAMETER_LIST(FEBulkElasticity, FEElasticMaterial)
ADD_PARAMETER(m_K, FE_PARAM_DOUBLE, "k");
END_PARAMETER_LIST();

//////////////////////////////////////////////////////////////////////
// Provide bulk modulus, e.g., to enforce incompressibility
//////////////////////////////////////////////////////////////////////

void FEBulkElasticity::Init()
{
	FEElasticMaterial::Init();
	
	if (m_unstable) throw MaterialError("This material is unstable when used alone.  Combine it in a solid mixture with a stable material.");
	if (m_K <= 0) throw MaterialError("k must be positive");
}

mat3ds FEBulkElasticity::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	double J = pt.J;
	mat3dd I(1); // Identity

	// calculate stress
	mat3ds s = m_K*log(J)/J*I;

	return s;
}

tens4ds FEBulkElasticity::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	double J = pt.J;
	mat3dd I(1); // Identity
	
	// calculate elasticity tensor
	tens4ds c = m_K/J*(dyad1s(I) - 2*log(J)*dyad4s(I));
	
	return c;
}
