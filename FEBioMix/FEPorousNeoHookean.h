#pragma once
#include <FEBioMech/FEElasticMaterial.h>
#include <FECore/FEModelParam.h>
#include "febiomix_api.h"

//-----------------------------------------------------------------------------
//! Porous Neo Hookean material
//! This compressible material produces infinite stress as the porosity approaches
//! zero (i.e., when J = 1 - porosity).  This behavior helps to prevent pore collape.

//! Implementation of a porous neo-Hookean hyperelastic material.
class FEBIOMIX_API FEPorousNeoHookean : public FEElasticMaterial
{
public:
    FEPorousNeoHookean(FEModel* pfem) : FEElasticMaterial(pfem) { m_phisr = 1; }
    
public:
    double  m_E;        //!< Young's modulus
    double  m_lam;      //!< first Lamé coefficient
    double  m_mu;       //!< second Lamé coefficient
    double  m_phiwr;    //!< referential porosity
    double  m_phisr;    //!< referential solidity
    
public:
    //! calculate stress at material point
    mat3ds Stress(FEMaterialPoint& pt) override;
    
    //! calculate tangent stiffness at material point
    tens4ds Tangent(FEMaterialPoint& pt) override;
    
    //! calculate strain energy density at material point
    double StrainEnergyDensity(FEMaterialPoint& pt) override;
    
    //! initialize
    bool Init() override;
    
    // declare the parameter list
    DECLARE_FECORE_CLASS();
};
