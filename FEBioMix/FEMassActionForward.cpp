/*
 *  FEMassActionForward.cpp
 *  FEBioXCode
 *
 *  Created by Gerard Ateshian on 3/6/13.
 *  Copyright 2013 Columbia University. All rights reserved.
 *
 */

#include "FEMassActionForward.h"

//-----------------------------------------------------------------------------
//! molar supply at material point
double FEMassActionForward::ReactionSupply(FEMaterialPoint& pt)
{
	FESolutesMaterialPoint& spt = *pt.ExtractData<FESolutesMaterialPoint>();
	
	// get reaction rate
	double kF = m_pFwd->ReactionRate(pt);
	
	// evaluate the reaction molar supply
	double zhat = kF;
	
	// start with contribution from solutes
	const int nsol = (int)spt.m_ca.size();
	for (int i=0; i<nsol; ++i) {
		int vR = m_vR[i];
		if (vR > 0) {
			double c = spt.m_ca[i];
			zhat *= pow(c, vR);
		}
	}
	
	// add contribution of solid-bound molecules
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
mat3ds FEMassActionForward::Tangent_ReactionSupply_Strain(FEMaterialPoint& pt)
{
	FEElasticMaterialPoint& ept = *pt.ExtractData<FEElasticMaterialPoint>();
	FEBiphasicMaterialPoint& bpt = *pt.ExtractData<FEBiphasicMaterialPoint>();
	FESolutesMaterialPoint& spt = *pt.ExtractData<FESolutesMaterialPoint>();
	
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
	for (int isol=0; isol<nsol; ++isol)
			dzhatde += I*(m_vR[isol]*spt.m_dkdJ[isol]/spt.m_k[isol]);
	for (int isbm = 0; isbm<nsbm; ++isbm)
		dzhatde -= I*(m_vR[nsol+isbm]/(J-phi0));
	
	dzhatde *= zhat;
	
	return dzhatde;
}

//-----------------------------------------------------------------------------
//! tangent of molar supply with effective pressure at material point
double FEMassActionForward::Tangent_ReactionSupply_Pressure(FEMaterialPoint& pt)
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
double FEMassActionForward::Tangent_ReactionSupply_Concentration(FEMaterialPoint& pt, const int sol)
{
	const int nsol = m_nsol;
	
	// if the derivative is taken with respect to a solid-bound molecule, return 0
	if (sol >= nsol) {
		return 0;
	}
	
	FESolutesMaterialPoint& spt = *pt.ExtractData<FESolutesMaterialPoint>();
	double zhat = ReactionSupply(pt);
	double dzhatdc = 0;
	for (int isol=0; isol<nsol; ++isol) {
		dzhatdc += m_vR[isol]*spt.m_dkdc[isol][sol]/spt.m_k[isol];
		if ((isol == sol) && (spt.m_c[sol] > 0))
			dzhatdc += m_vR[isol]/spt.m_c[sol];
	}
	
	dzhatdc *= zhat;
	
	return dzhatdc;
}
