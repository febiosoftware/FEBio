/*
 *  FEMichaelisMenten.cpp
 *  FEBioXCode
 *
 *  Created by Gerard Ateshian on 3/8/13.
 *  Copyright 2013 Columbia University. All rights reserved.
 *
 */

#include "FEMichaelisMenten.h"

// define the material parameters
BEGIN_PARAMETER_LIST(FEMichaelisMenten, FEChemicalReaction)
	ADD_PARAMETER(m_Km, FE_PARAM_DOUBLE, "Km");
	ADD_PARAMETER(m_c0, FE_PARAM_DOUBLE, "c0");
END_PARAMETER_LIST();

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

//-----------------------------------------------------------------------------
//! data initialization and checking
bool FEMichaelisMenten::Init()
{
    // Initialize base class
    if (FEChemicalReaction::Init() == false) return false;
    
	// there is only one reactant and one product in a Michaelis-Menten reaction
	if (m_solR.size() + m_sbmR.size() > 1)
		return MaterialError("Provide only one vR for this reaction");
	if (m_solP.size() + m_sbmP.size() > 1)
		return MaterialError("Provide only one vP for this reaction");

	if (m_c0 < 0) return MaterialError("c0 must be positive");
	
	const int ntot = (int)m_v.size();
	for (int itot=0; itot<ntot; itot++) {
		if (m_vR[itot] > 0) m_Rid = itot;
		if (m_vP[itot] > 0) m_Pid = itot;
	}
	
	if (m_Rid == -1) return MaterialError("Provide vR for the reactant");
	
	// check if reactant is a solute or a solid-bound molecule
	if (m_Rid >= m_nsol) m_Rtype = true;
	
	return true;
}

//-----------------------------------------------------------------------------
//! molar supply at material point
double FEMichaelisMenten::ReactionSupply(FEMaterialPoint& pt)
{
	FESolutesMaterialPoint& spt = *pt.ExtractData<FESolutesMaterialPoint>();

	// get reaction rate
	double Vmax = m_pFwd->ReactionRate(pt);
	double c;
	if (m_Rtype) {
		c = m_pMP->SBMConcentration(pt, m_Rid);
	}
	else {
		c = spt.m_ca[m_Rid];
	}

	double zhat = 0;
	if (c > m_c0) zhat = Vmax*c/(m_Km + c);
	
	return zhat;
}

//-----------------------------------------------------------------------------
//! tangent of molar supply with strain at material point
mat3ds FEMichaelisMenten::Tangent_ReactionSupply_Strain(FEMaterialPoint& pt)
{
	FEElasticMaterialPoint& ept = *pt.ExtractData<FEElasticMaterialPoint>();
	FEBiphasicMaterialPoint& bpt = *pt.ExtractData<FEBiphasicMaterialPoint>();
	FESolutesMaterialPoint& spt = *pt.ExtractData<FESolutesMaterialPoint>();
	
	double c;
	double dcdJ;
	if (m_Rtype) {
		c = m_pMP->SBMConcentration(pt, m_Rid);
		double J = ept.m_J;
		double phi0 = bpt.m_phi0;
		dcdJ = -c/(J-phi0);
	}
	else {
		c = spt.m_ca[m_Rid];
		dcdJ = spt.m_dkdJ[m_Rid]*spt.m_c[m_Rid];
	}
	
	double dzhatdJ = 0;
	if (c > m_c0) {
        double Vmax = m_pFwd->ReactionRate(pt);
        dzhatdJ = dcdJ*m_Km*Vmax/SQR(m_Km + c);
    }
	
	return mat3dd(1)*dzhatdJ;
}

//-----------------------------------------------------------------------------
//! tangent of molar supply with effective pressure at material point
double FEMichaelisMenten::Tangent_ReactionSupply_Pressure(FEMaterialPoint& pt)
{
	return 0;
}

//-----------------------------------------------------------------------------
//! tangent of molar supply with effective concentration at material point
double FEMichaelisMenten::Tangent_ReactionSupply_Concentration(FEMaterialPoint& pt, const int sol)
{
	FESolutesMaterialPoint& spt = *pt.ExtractData<FESolutesMaterialPoint>();
	
	if (m_Rtype) {
        return 0;
	}
	else if (m_Rid != sol)
        return 0;
	
    double c = spt.m_ca[m_Rid];
	double dzhatdc = 0;
	if (c > m_c0) {
        double Vmax = m_pFwd->ReactionRate(pt);
        dzhatdc = m_Km*Vmax/SQR(m_Km + c)*(spt.m_k[m_Rid] + spt.m_dkdc[m_Rid][m_Rid]*spt.m_c[m_Rid]);
    }
	
	return dzhatdc;
}
