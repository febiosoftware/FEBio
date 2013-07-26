/*
 *  FEMichaelisMenten.cpp
 *  FEBioXCode
 *
 *  Created by Gerard Ateshian on 3/8/13.
 *  Copyright 2013 Columbia University. All rights reserved.
 *
 */

#include "FEMichaelisMenten.h"

// register the material with the framework
REGISTER_MATERIAL(FEMichaelisMenten, "Michaelis-Menten");

// define the material parameters
BEGIN_PARAMETER_LIST(FEMichaelisMenten, FEChemicalReaction)
	ADD_PARAMETER(m_Km, FE_PARAM_DOUBLE, "Km");
	ADD_PARAMETER(m_c0, FE_PARAM_DOUBLE, "c0");
END_PARAMETER_LIST();


//-----------------------------------------------------------------------------
//! data initialization and checking
void FEMichaelisMenten::Init()
{
    // Initialize base class
    FEChemicalReaction::Init();
    
	// there is only one reactant and one product in a Michaelis-Menten reaction
	if (m_solR.size() + m_sbmR.size() > 1)
		throw MaterialError("Provide only one vR for this reaction");
	if (m_solP.size() + m_sbmP.size() > 1)
		throw MaterialError("Provide only one vP for this reaction");

	if (m_c0 < 0) throw MaterialError("c0 must be positive");
	
	const int ntot = (int)m_v.size();
	for (int itot=0; itot<ntot; itot++) {
		if (m_vR[itot] > 0) m_Rid = itot;
		if (m_vP[itot] > 0) m_Pid = itot;
	}
	
	if (m_Rid == -1) throw MaterialError("Provide vR for the reactant");
	
	// check if reactant is a solute or a solid-bound molecule
	if (m_Rid >= m_nsol) m_Rtype = true;
	
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
	
	double zhat = ReactionSupply(pt);
	double dzhatdJ = 0;
	if (c > m_c0) dzhatdJ = dcdJ*m_Km*zhat/(m_Km + c);
	
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
	
	double c;
	double dcdc;
	if (m_Rtype) {
		c = m_pMP->SBMConcentration(pt, m_Rid);
		dcdc = (m_Rid == sol) ? 1 : 0;
	}
	else {
		c = spt.m_ca[m_Rid];
		dcdc = (m_Rid == sol) ? 1 : 0;
	}
	
	double zhat = ReactionSupply(pt);
	double dzhatdc = 0;
	if (c > m_c0) dzhatdc = dcdc*m_Km*zhat/(m_Km + c);
	
	return dzhatdc;
}
