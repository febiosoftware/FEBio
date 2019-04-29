/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, The Trustees of Columbia University in 
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/



#include "stdafx.h"
#include "FEMichaelisMenten.h"
#include <FECore/log.h>

// define the material parameters
BEGIN_FECORE_CLASS(FEMichaelisMenten, FEChemicalReaction)
	ADD_PARAMETER(m_Km, "Km");
	ADD_PARAMETER(m_c0, "c0");
END_FECORE_CLASS();

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
	if (m_solR.size() + m_sbmR.size() > 1) {
		feLogError("Provide only one vR for this reaction");
		return false;
	}

	if (m_solP.size() + m_sbmP.size() > 1) {
		feLogError("Provide only one vP for this reaction");
		return false;
	}

	if (m_c0 < 0) {
		feLogError("c0 must be positive");
		return false;
	}
	
	const int ntot = (int)m_v.size();
	for (int itot=0; itot<ntot; itot++) {
		if (m_vR[itot] > 0) m_Rid = itot;
		if (m_vP[itot] > 0) m_Pid = itot;
	}
	
	if (m_Rid == -1) {
		feLogError("Provide vR for the reactant");
		return false;
	}
	
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
