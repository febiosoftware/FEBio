/*
 *  FEReactionRateHuiskes.cpp
 *  FEBioXCode
 *
 *  Created by Gerard Ateshian on 5/15/13.
 *  Copyright 2013 Columbia University. All rights reserved.
 *
 */

#include "stdafx.h"
#include "FEReactionRateHuiskes.h"

// register the material with the framework
REGISTER_MATERIAL(FEReactionRateHuiskes, "Huiskes");

// Material parameters for the FEMultiphasic material
BEGIN_PARAMETER_LIST(FEReactionRateHuiskes, FEMaterial)
ADD_PARAMETER(m_B, FE_PARAM_DOUBLE, "B");
ADD_PARAMETER(m_k, FE_PARAM_DOUBLE, "k");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
void FEReactionRateHuiskes::Init()
{
	FEMaterial::Init();
	
	if (m_k < 0) throw MaterialError("k must be positive");
	
}

//-----------------------------------------------------------------------------
//! reaction rate at material point
double FEReactionRateHuiskes::ReactionRate(FEMaterialPoint& pt)
{
	double rhor = m_pReact->m_pMP->SolidReferentialApparentDensity(pt);
    double phir = m_pReact->m_pMP->SolidReferentialVolumeFraction(pt);
	
	FEElasticMaterialPoint& et = *pt.ExtractData<FEElasticMaterialPoint>();
    double J = et.J;
	double sed = et.sed;
	double zhat = m_B*(sed/rhor - m_k)/(J-phir);
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
    double J = et.J;
    double p = bt.m_pa;
    double zhat = ReactionRate(pt);
    mat3dd I(1);
    mat3ds dzhatde = (I*(-zhat) + (et.s+I*p)*(m_B/rhor))/(J-phir);
	return dzhatde;
}

//-----------------------------------------------------------------------------
//! tangent of reaction rate with effective fluid pressure at material point
double FEReactionRateHuiskes::Tangent_ReactionRate_Pressure(FEMaterialPoint& pt)
{
	return 0;
}

