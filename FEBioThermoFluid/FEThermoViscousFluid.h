#pragma once
#include <FEBioFluid/FEViscousFluid.h>
#include <FECore/FEFunction1D.h>
#include "febiothermofluid_api.h"

class FEBIOTHERMOFLUID_API FEThermoViscousFluid : public FEViscousFluid
{
public:
    enum { MAX_PROPS = 6 };
    
public:
    FEThermoViscousFluid(FEModel* fem);

    //! initialization
    bool Init() override;
    
    //! Serialization
    void Serialize(DumpStream& ar) override;
    
	//! tangent of stress with respect to temperature
	virtual mat3ds Tangent_Temperature(FEMaterialPoint& mp) = 0;
    
    double ThermoProp(FEMaterialPoint& mp, const int iprop);
    double TangentThermoProp(FEMaterialPoint& mp, const int iprop);

public:
    double  m_Tr;   //!< referential temperature
    FEFunction1D*   m_prop[MAX_PROPS];
};
