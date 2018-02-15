//
//  FENewtonianViscousSolidUC.cpp
//  FEBioMech
//
//  Created by Gerard Ateshian on 6/14/17.
//  Copyright © 2017 febio.org. All rights reserved.
//

#include "FENewtonianViscousSolidUC.h"
#include <FECore/FEModel.h>

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_PARAMETER_LIST(FENewtonianViscousSolidUC, FEUncoupledMaterial)
ADD_PARAMETER2(m_kappa, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(      0.0), "kappa");
ADD_PARAMETER2(m_mu   , FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(      0.0), "mu"   );
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
mat3ds FENewtonianViscousSolidUC::DevStress(FEMaterialPoint& mp)
{
    FEElasticMaterialPoint& pe = *mp.ExtractData<FEElasticMaterialPoint>();
    
    mat3ds D = pe.RateOfDeformation();
    
    // Identity
    mat3dd I(1);
    
    // calculate stress
    mat3ds s = I*(D.tr()*(m_kappa - 2*m_mu/3)) + D*(2*m_mu);
    
    // determinant of deformation gradient
    double J = pe.m_J;
    
    return s.dev()*(2.0/J);
}

//-----------------------------------------------------------------------------
tens4ds FENewtonianViscousSolidUC::DevTangent(FEMaterialPoint& mp)
{
    FEElasticMaterialPoint& pe = *mp.ExtractData<FEElasticMaterialPoint>();
    
    mat3dd I(1);
    
    double dt = GetFEModel()->GetTime().timeIncrement;
    tens4ds Cv = (dyad1s(I, I)*(m_kappa - 2 * m_mu / 3) + dyad4s(I, I)*(2 * m_mu)) / (2 * dt);
    
    return Cv;
}

//-----------------------------------------------------------------------------
double FENewtonianViscousSolidUC::DevStrainEnergyDensity(FEMaterialPoint& mp)
{
    return 0;
}

