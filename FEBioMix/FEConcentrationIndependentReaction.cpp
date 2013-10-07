//
//  FEConcentrationIndependentReaction.cpp
//  FEBioXCode4
//
//  Created by Gerard Ateshian on 7/21/13.
//  Copyright (c) 2013 Columbia University. All rights reserved.
//

#include "stdafx.h"
#include "FEConcentrationIndependentReaction.h"

//-----------------------------------------------------------------------------
//! molar supply at material point
double FEConcentrationIndependentReaction::ReactionSupply(FEMaterialPoint& pt)
{
	FESolutesMaterialPoint& spt = *pt.ExtractData<FESolutesMaterialPoint>();
	
	// get reaction rate
	double kF = m_pFwd->ReactionRate(pt);
	
	// evaluate the reaction molar supply
	double zhat = kF;
	
	// contribution of solid-bound molecules
	const int nsol = (int)spt.m_ca.size();
	const int nsbm = (int)spt.m_sbmr.size();
	for (int i=0; i<nsbm; ++i) {
		int vR = m_vR[nsol+i];
		if (vR > 0) {
			double c = m_pMP->SBMConcentration(pt, i);
			zhat *= pow(c, vR);
		}
	}
	
	return zhat;
}

//-----------------------------------------------------------------------------
//! tangent of molar supply with strain at material point
mat3ds FEConcentrationIndependentReaction::Tangent_ReactionSupply_Strain(FEMaterialPoint& pt)
{
	FEElasticMaterialPoint& ept = *pt.ExtractData<FEElasticMaterialPoint>();
	FEBiphasicMaterialPoint& bpt = *pt.ExtractData<FEBiphasicMaterialPoint>();
	
	const int nsol = m_nsol;
	const int nsbm = (int)m_v.size() - nsol;
	double J = ept.m_J;
	double phi0 = bpt.m_phi0;
    
	double kF = m_pFwd->ReactionRate(pt);
	mat3ds dkFde = m_pFwd->Tangent_ReactionRate_Strain(pt);
	double zhat = ReactionSupply(pt);
	mat3ds dzhatde = mat3dd(0);
	if (kF > 0) {
		dzhatde += dkFde/kF;
	}
	mat3ds I = mat3dd(1);
	for (int isbm = 0; isbm<nsbm; ++isbm)
		dzhatde -= I*(m_vR[nsol+isbm]/(J-phi0));
	
	dzhatde *= zhat;
	
	return dzhatde;
}

//-----------------------------------------------------------------------------
//! tangent of molar supply with effective pressure at material point
double FEConcentrationIndependentReaction::Tangent_ReactionSupply_Pressure(FEMaterialPoint& pt)
{
	double kF = m_pFwd->ReactionRate(pt);
	double dkFdp = m_pFwd->Tangent_ReactionRate_Pressure(pt);
	double zhat = ReactionSupply(pt);
	double dzhatdp = 0;
	if (kF > 0) {
		dzhatdp = dkFdp*zhat/kF;
	}
	return dzhatdp;
}

//-----------------------------------------------------------------------------
//! tangent of molar supply with effective concentration at material point
double FEConcentrationIndependentReaction::Tangent_ReactionSupply_Concentration(FEMaterialPoint& pt, const int sol)
{
	return 0;
}
