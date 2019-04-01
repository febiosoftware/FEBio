#pragma once
#include "FEViscousFluid.h"

//-----------------------------------------------------------------------------
// This class evaluates the viscous stress in a Carreau-Yasuda power-law fluid

class FEBIOFLUID_API FECarreauYasudaFluid :	public FEViscousFluid
{
public:
    //! constructor
    FECarreauYasudaFluid(FEModel* pfem);
    
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
    double	m_mu0;		//!< shear viscosity at zero shear rate
    double	m_mui;		//!< shear viscosity at infinite shear rate
    double  m_lam;      //!< time constant
    double  m_n;        //!< power-law index
    double  m_a;        //!< exponent
    
    // declare parameter list
    DECLARE_FECORE_CLASS();
};
