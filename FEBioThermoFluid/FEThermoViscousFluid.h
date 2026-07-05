#pragma once

#include <FEBioFluid/FEViscousFluid.h>
#include <FECore/FEFunction1D.h>
#include "febiothermofluid_api.h"

class FEBIOTHERMOFLUID_API FEThermoViscousFluid : public FEViscousFluid
{
public:
    FEThermoViscousFluid(FEModel* fem);

	//! tangent of stress with respect to temperature
	virtual mat3ds Tangent_Temperature(FEMaterialPoint& mp) = 0;
    
public:
    double  m_Tr;   //!< referential temperature
};
