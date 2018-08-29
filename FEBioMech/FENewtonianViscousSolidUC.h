//
//  FENewtonianViscousSolidUC.hpp
//  FEBioMech
//
//  Created by Gerard Ateshian on 6/14/17.
//  Copyright © 2017 febio.org. All rights reserved.
//

#ifndef FENewtonianViscousSolidUC_hpp
#define FENewtonianViscousSolidUC_hpp

#include "FEUncoupledMaterial.h"

class FENewtonianViscousSolidUC : public FEUncoupledMaterial
{
public:
    FENewtonianViscousSolidUC(FEModel* pfem) : FEUncoupledMaterial(pfem) {}
    
public:
    double	m_kappa;	//!< bulk viscosity
    double	m_mu;       //!< shear viscosity
    
public:
    //! calculate stress at material point
    mat3ds DevStress(FEMaterialPoint& pt) override;
    
    //! calculate tangent stiffness at material point
    tens4ds DevTangent(FEMaterialPoint& pt) override;
    
    //! calculate strain energy density at material point
    double DevStrainEnergyDensity(FEMaterialPoint& pt) override;
    
    // declare the parameter list
    DECLARE_PARAMETER_LIST();
};

#endif /* FENewtonianViscousSolidUC_hpp */
