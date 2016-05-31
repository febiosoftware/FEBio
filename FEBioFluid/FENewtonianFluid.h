#pragma once
#include "FEViscousFluid.h"

//-----------------------------------------------------------------------------
// This class evaluates the viscous stress in a Newtonian fluid

class FENewtonianFluid :	public FEViscousFluid
{
public:
    //! constructor
    FENewtonianFluid(FEModel* pfem);
    
    //! viscous stress
    mat3ds Stress(FEMaterialPoint& pt);
    
    //! tangent of stress with respect to strain J
    mat3ds Tangent_Strain(FEMaterialPoint& mp);
    
    //! tangent of stress with respect to rate of deformation tensor D
    tens4ds Tangent_RateOfDeformation(FEMaterialPoint& mp);
    
    //! dynamic viscosity
    double ShearViscosity(FEMaterialPoint& mp);
    
    //! bulk viscosity
    double BulkViscosity(FEMaterialPoint& mp);
    
public:
    double	m_kappa;	//!< bulk viscosity
    double	m_mu;		//!< shear viscosity
    
    // declare parameter list
    DECLARE_PARAMETER_LIST();
};
