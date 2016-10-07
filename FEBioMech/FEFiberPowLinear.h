#pragma once
#include "FEElasticFiberMaterial.h"

//-----------------------------------------------------------------------------
//! Material class for single fiber, tension only
//! Power law - linear

class FEFiberPowLinear : public FEElasticFiberMaterial
{
public:
    FEFiberPowLinear(FEModel* pfem);
    
    //! Cauchy stress
    mat3ds Stress(FEMaterialPoint& mp);
    
    // Spatial tangent
    tens4ds Tangent(FEMaterialPoint& mp);
    
    //! Strain energy density
    double StrainEnergyDensity(FEMaterialPoint& mp);
    
    // declare the parameter list
    DECLARE_PARAMETER_LIST();
    
public:
    double	m_E;		// fiber modulus
    double  m_lam0;     // stretch ratio at end of toe region
    double  m_beta;     // power law exponent in toe region
};
