#pragma once
#include <FECore/FEMaterial.h>
#include <FECore/tens4d.h>

//-----------------------------------------------------------------------------
//! Base class for the viscous part of the fluid response.
//! These materials provide the viscous stress and its tangents.
//!
class FEViscousFluid : public FEMaterial
{
public:
    FEViscousFluid(FEModel* pfem) : FEMaterial(pfem) {}
    virtual ~FEViscousFluid() {}
    
    //! viscous stress
    virtual mat3ds Stress(FEMaterialPoint& pt) = 0;
    
    //! tangent of stress with respect to strain J
    virtual mat3ds Tangent_Strain(FEMaterialPoint& mp) = 0;
    
    //! tangent of stress with respect to rate of deformation tensor D
    virtual tens4ds Tangent_RateOfDeformation(FEMaterialPoint& mp) = 0;
    
    //! dynamic viscosity
    virtual double ShearViscosity(FEMaterialPoint& mp) = 0;
    
    //! bulk viscosity
    virtual double BulkViscosity(FEMaterialPoint& mp) = 0;
    
};
