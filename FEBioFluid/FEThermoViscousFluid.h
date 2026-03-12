#pragma once
#include "FEViscousFluid.h"
#include "febiofluid_api.h"

class FEBIOFLUID_API FEThermoViscousFluid : public FEViscousFluid
{
public:
	FEThermoViscousFluid(FEModel* fem) : FEViscousFluid(fem) {}

	//! tangent of stress with respect to temperature
	virtual mat3ds Tangent_Temperature(FEMaterialPoint& mp) = 0;
};
