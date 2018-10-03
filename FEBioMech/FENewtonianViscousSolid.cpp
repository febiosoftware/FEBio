//
//  FENewtonianViscousSolid.cpp
//  FEBioMech
//
//  Created by Gerard Ateshian on 4/7/16.
//  Copyright © 2016 febio.org. All rights reserved.
//

#include "FENewtonianViscousSolid.h"
#include <FECore/FEModel.h>

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FENewtonianViscousSolid, FEElasticMaterial)
	ADD_PARAMETER(m_kappa, FE_RANGE_GREATER_OR_EQUAL(      0.0), "kappa");
	ADD_PARAMETER(m_mu   , FE_RANGE_GREATER_OR_EQUAL(      0.0), "mu"   );
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
mat3ds FENewtonianViscousSolid::Stress(FEMaterialPoint& mp)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    mat3ds D = pt.RateOfDeformation();
    
    // Identity
    mat3dd I(1);
    
    // calculate stress
    mat3ds s = I*(D.tr()*(m_kappa - 2*m_mu/3)) + D*(2*m_mu);
    
    return s;
}

//-----------------------------------------------------------------------------
tens4ds FENewtonianViscousSolid::Tangent(FEMaterialPoint& mp)
{
    FEElasticMaterialPoint& pe = *mp.ExtractData<FEElasticMaterialPoint>();
    
    mat3dd I(1);
    
	double dt = GetFEModel()->GetTime().timeIncrement;
	tens4ds Cv = (dyad1s(I, I)*(m_kappa - 2 * m_mu / 3) + dyad4s(I, I)*(2 * m_mu)) / (2 * dt);
    
    return Cv;
}

//-----------------------------------------------------------------------------
double FENewtonianViscousSolid::StrainEnergyDensity(FEMaterialPoint& mp)
{
    return 0;
}

