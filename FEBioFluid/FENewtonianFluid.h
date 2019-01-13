#pragma once
#include "FEViscousFluid.h"

//-----------------------------------------------------------------------------
// This class evaluates the viscous stress in a Newtonian fluid

class FEBIOFLUID_API FENewtonianFluid :	public FEViscousFluid
{
public:
    //! constructor
    FENewtonianFluid(FEModel* pfem);
    
    //! viscous stress
    mat3ds Stress(FEMaterialPoint& pt) override;
    
    //! tangent of stress with respect to strain J
    mat3ds Tangent_Strain(FEMaterialPoint& mp) override;
    
    //! tangent of stress with respect to rate of deformation tensor D
    tens4ds Tangent_RateOfDeformation(FEMaterialPoint& mp) override;
    
    //! dynamic viscosity
    double ShearViscosity(FEMaterialPoint& mp) override;
    
    //! bulk viscosity
    double BulkViscosity(FEMaterialPoint& mp) override;
    
public:
    double	m_kappa;	//!< bulk viscosity
    double	m_mu;		//!< shear viscosity
    
    // declare parameter list
    DECLARE_FECORE_CLASS();
};
