#pragma once
#include <FECore/FEMaterial.h>
#include "febiomix_api.h"

//-----------------------------------------------------------------------------
//! Base class for active momentum supply.
//! These materials need to define the momentum supply and its tangents.
//!
class FEBIOMIX_API FEActiveMomentumSupply : public FEMaterial
{
public:
    FEActiveMomentumSupply(FEModel* pfem) : FEMaterial(pfem) {}
    virtual ~FEActiveMomentumSupply(){}
    
    //! active momentum supply
    virtual vec3d ActiveSupply(FEMaterialPoint& pt) = 0;
    
    //! tangent of active momentum supply with respect to strain
    virtual vec3d Tangent_ActiveSupply_Strain(FEMaterialPoint& mp) = 0;
    
    //! tangent of hydraulic permeability with respect to concentration
    vec3d Tangent_ActiveSupply_Concentration(FEMaterialPoint& mp, const int isol);
};
