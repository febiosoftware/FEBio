#pragma once
#include <FEBioMech/FEElasticMaterial.h>

//-----------------------------------------------------------------------------
//! This class is the base class of elastic materials that are used in a thermo-elastic
//! analysis. These materials can be dependent on the temperature and must also implement
//! the thermal tangent, i.e. the derivative of stress with respect to the temperature.
class FEThermalElastic : public FEElasticMaterial
{
public:
	FEThermalElastic(FEModel* pfem) : FEElasticMaterial(pfem) {}

	//! spatial thermal tangent, i.e. derivative of Caucy stress wrt temperature
	virtual mat3ds ThermalTangent(FEMaterialPoint& mp) = 0;
};
