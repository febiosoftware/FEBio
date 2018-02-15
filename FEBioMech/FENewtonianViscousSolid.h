//
//  FENewtonianViscousSolid.hpp
//  FEBioMech
//
//  Created by Gerard Ateshian on 4/7/16.
//  Copyright © 2016 febio.org. All rights reserved.
//

#ifndef FENewtonianViscousSolid_hpp
#define FENewtonianViscousSolid_hpp

#include "FEElasticMaterial.h"

class FENewtonianViscousSolid : public FEElasticMaterial
{
public:
    FENewtonianViscousSolid(FEModel* pfem) : FEElasticMaterial(pfem) {}
    
public:
    double	m_kappa;	//!< bulk viscosity
    double	m_mu;       //!< shear viscosity
    
public:
    //! calculate stress at material point
    virtual mat3ds Stress(FEMaterialPoint& pt) override;
    
    //! calculate tangent stiffness at material point
    virtual tens4ds Tangent(FEMaterialPoint& pt) override;
    
    //! calculate strain energy density at material point
    virtual double StrainEnergyDensity(FEMaterialPoint& pt) override;
    
    // declare the parameter list
    DECLARE_PARAMETER_LIST();
};

#endif /* FENewtonianViscousSolid_hpp */
