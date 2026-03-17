#pragma once
#include <FEBioFluid/FEViscousFluid.h>
#include <FECore/FEFunction1D.h>
#include <FECore/FEMaterial.h>
#include "febiothermofluid_api.h"

class FEBIOTHERMOFLUID_API FEFluidHeatFlux : public FEMaterial
{
public:
    FEFluidHeatFlux(FEModel* fem);

	//! tangent of stress with respect to temperature
	virtual vec3d HeatFlux(FEMaterialPoint& mp) = 0;
    
    virtual double Conductivity(FEMaterialPoint& mp) = 0;
    
    virtual double Tangent_Conductivity_Temperature(FEMaterialPoint& mp) = 0;
    
    virtual double Tangent_Conductivity_Strain(FEMaterialPoint& mp) = 0;
};
