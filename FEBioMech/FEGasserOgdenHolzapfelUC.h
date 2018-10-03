#pragma once
#include "FEUncoupledMaterial.h"

class FEGasserOgdenHolzapfelUC : public FEUncoupledMaterial
{
public:
    double	m_c;			// neo-Hookean c coefficient
    double	m_k1,m_k2;		// fiber material constants
    double	m_kappa;		// structure coefficient
    double	m_g;			// fiber angle
    
public:
    FEGasserOgdenHolzapfelUC(FEModel* pfem) : FEUncoupledMaterial(pfem) {}
    
    //! calculate deviatoric stress at material point
    mat3ds DevStress(FEMaterialPoint& pt) override;
    
    //! calculate deviatoric tangent stiffness at material point
    tens4ds DevTangent(FEMaterialPoint& pt) override;
    
    //! calculate deviatoric strain energy density at material point
    double DevStrainEnergyDensity(FEMaterialPoint& pt) override;
    
    // declare parameter list
    DECLARE_FECORE_CLASS();
};
