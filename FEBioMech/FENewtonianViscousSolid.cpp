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
BEGIN_PARAMETER_LIST(FENewtonianViscousSolid, FEElasticMaterial)
ADD_PARAMETER2(m_kappa, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(      0.0), "kappa");
ADD_PARAMETER2(m_mu   , FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(      0.0), "mu"   );
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
mat3ds FENewtonianViscousSolid::Stress(FEMaterialPoint& mp)
{
    FEViscousMaterialPoint& pt = *mp.ExtractData<FEViscousMaterialPoint>();
    
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
    FEViscousMaterialPoint& pv = *mp.ExtractData<FEViscousMaterialPoint>();
    
    mat3ds A = (pe.m_F.inverse()*pv.m_Fp).sym();
    
    mat3dd I(1);
    
	double dt = GetFEModel()->GetTime().timeIncrement;
	tens4ds Cv = (dyad1s(I, A)*(m_kappa - 2 * m_mu / 3) + dyad4s(I, A)*(2 * m_mu)) / (2 * dt);
    
    return Cv;
}

//-----------------------------------------------------------------------------
double FENewtonianViscousSolid::StrainEnergyDensity(FEMaterialPoint& mp)
{
    return 0;
}

