#pragma once
#include "FEElasticMaterial.h"

//-----------------------------------------------------------------------------
//! Material class that implements osmotic pressure using a virial expansion.
//
class FEOsmoticVirialExpansion : public FEElasticMaterial
{
public:
	FEOsmoticVirialExpansion(FEModel* pfem);
    
    //! Returns the Cauchy stress
    mat3ds Stress(FEMaterialPoint& mp) override;
    
    //! Returs the spatial tangent
    tens4ds Tangent(FEMaterialPoint& mp) override;
    
    // declare the parameter list
    DECLARE_FECORE_CLASS();
    
public:
    double	m_phiwr;	//!< fluid volume fraction in reference configuration
    double	m_cr;		//!< concentration in reference configuration
    double  m_c1;       //!< first virial coefficient
    double  m_c2;       //!< second virial coefficient
    double  m_c3;       //!< third virial coefficient
};
