/*
 *  FEReactionRateHuiskes.cpp
 *  FEBioXCode
 *
 *  Created by Gerard Ateshian on 5/15/13.
 *  Copyright 2013 Columbia University. All rights reserved.
 *
 */

#include "FEReactionRateHuiskes.h"
#include "FEBioMech/FERemodelingElasticMaterial.h"

// Material parameters for the FEMultiphasic material
BEGIN_PARAMETER_LIST(FEReactionRateHuiskes, FEMaterial)
	ADD_PARAMETER(m_B, FE_PARAM_DOUBLE, "B");
	ADD_PARAMETER(m_psi0, FE_PARAM_DOUBLE, "psi0");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
void FEReactionRateHuiskes::Init()
{
	FEMaterial::Init();
	
	if (m_psi0 < 0) throw MaterialError("k must be positive");

}

//-----------------------------------------------------------------------------
//! reaction rate at material point
double FEReactionRateHuiskes::ReactionRate(FEMaterialPoint& pt)
{
	double rhor = m_pReact->m_pMP->SolidReferentialApparentDensity(pt);
    double phir = m_pReact->m_pMP->SolidReferentialVolumeFraction(pt);
	
    FERemodelingMaterialPoint& rpt = *(pt.ExtractData<FERemodelingMaterialPoint>());
	FEElasticMaterialPoint& et = *pt.ExtractData<FEElasticMaterialPoint>();
    double J = et.m_J;
	double sed = rpt.m_sed;
	double zhat = m_B*(sed/rhor - m_psi0)/(J-phir);
	return zhat;
}

//-----------------------------------------------------------------------------
//! tangent of reaction rate with strain at material point
mat3ds FEReactionRateHuiskes::Tangent_ReactionRate_Strain(FEMaterialPoint& pt)
{
	double rhor = m_pReact->m_pMP->SolidReferentialApparentDensity(pt);
    double phir = m_pReact->m_pMP->SolidReferentialVolumeFraction(pt);
	
	FEElasticMaterialPoint& et = *pt.ExtractData<FEElasticMaterialPoint>();
    FEBiphasicMaterialPoint& bt = *pt.ExtractData<FEBiphasicMaterialPoint>();
    double J = et.m_J;
    double p = bt.m_pa;
    double zhat = ReactionRate(pt);
    mat3dd I(1);
    mat3ds dzhatde = (I*(-zhat) + (et.m_s+I*p)*(m_B/rhor))/(J-phir);
	return dzhatde;
}

//-----------------------------------------------------------------------------
//! tangent of reaction rate with effective fluid pressure at material point
double FEReactionRateHuiskes::Tangent_ReactionRate_Pressure(FEMaterialPoint& pt)
{
	return 0;
}

