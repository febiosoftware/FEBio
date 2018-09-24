//
//  FEReactionRateExpSED.cpp
//  FEBioMix
//
//  Created by Gerard Ateshian on 2/3/15.
//  Copyright (c) 2015 febio.org. All rights reserved.
//

#include "FEReactionRateExpSED.h"
#include "FEBioMech/FERemodelingElasticMaterial.h"

// Material parameters for the FEMultiphasic material
BEGIN_PARAMETER_LIST(FEReactionRateExpSED, FEMaterial)
	ADD_PARAMETER(m_B   , "B");
	ADD_PARAMETER(m_Psi0, FE_RANGE_NOT_EQUAL(0.0), "Psi0");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! reaction rate at material point
double FEReactionRateExpSED::ReactionRate(FEMaterialPoint& pt)
{
    double phir = m_pReact->m_pMP->SolidReferentialVolumeFraction(pt);
    
    FERemodelingMaterialPoint& rpt = *(pt.ExtractData<FERemodelingMaterialPoint>());
    FEElasticMaterialPoint& et = *pt.ExtractData<FEElasticMaterialPoint>();
    double J = et.m_J;
    double sed = rpt.m_sed;
    double zhat = m_B*exp(sed/m_Psi0)/(J-phir);
    return zhat;
}

//-----------------------------------------------------------------------------
//! tangent of reaction rate with strain at material point
mat3ds FEReactionRateExpSED::Tangent_ReactionRate_Strain(FEMaterialPoint& pt)
{
    double phir = m_pReact->m_pMP->SolidReferentialVolumeFraction(pt);
    
    FEElasticMaterialPoint& et = *pt.ExtractData<FEElasticMaterialPoint>();
    FEBiphasicMaterialPoint& bt = *pt.ExtractData<FEBiphasicMaterialPoint>();
    double J = et.m_J;
    double p = bt.m_pa;
    double zhat = ReactionRate(pt);
    mat3dd I(1);
    mat3ds dzhatde = (I/(phir-J) + (et.m_s+I*p)/m_Psi0)*zhat;
    return dzhatde;
}

//-----------------------------------------------------------------------------
//! tangent of reaction rate with effective fluid pressure at material point
double FEReactionRateExpSED::Tangent_ReactionRate_Pressure(FEMaterialPoint& pt)
{
    return 0;
}

