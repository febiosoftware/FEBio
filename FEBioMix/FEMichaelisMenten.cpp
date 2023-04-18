/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
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
#include "FESoluteInterface.h"
#include "FEMultiphasic.h"
#include <FECore/log.h>

// define the material parameters
BEGIN_FECORE_CLASS(FEMichaelisMenten, FEChemicalReaction)
	ADD_PARAMETER(m_Km, "Km");
	ADD_PARAMETER(m_c0, "c0");

	// set material properties
	ADD_PROPERTY(m_pFwd, "forward_rate", FEProperty::Optional);

END_FECORE_CLASS();

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

//-----------------------------------------------------------------------------
FEMichaelisMenten::FEMichaelisMenten(FEModel* pfem) : FEChemicalReaction(pfem) 
{ 
	m_Rid = m_Pid = -1; 
	m_Km = m_c0 = 0; 
	m_Rtype = false; 
}

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
	// get reaction rate
	double Vmax = m_pFwd->ReactionRate(pt);
	double c = 0.0;
	if (m_Rtype) {
		c = m_psm->SBMConcentration(pt, m_Rid);
	}
	else {
		c = m_psm->GetActualSoluteConcentration(pt, m_Rid);
	}

	double zhat = 0;
	if (c > m_c0) zhat = Vmax*c/(m_Km + c);
	
	return zhat;
}

//-----------------------------------------------------------------------------
//! tangent of molar supply with strain at material point
mat3ds FEMichaelisMenten::Tangent_ReactionSupply_Strain(FEMaterialPoint& pt)
{
	double ca = 0.0;
	double dcdJ = 0.0;
	if (m_Rtype) {
		FEBiphasicInterface* pbm = dynamic_cast<FEBiphasicInterface*>(GetAncestor()); assert(pbm);
		FEElasticMaterialPoint& ept = *pt.ExtractData<FEElasticMaterialPoint>();
		ca = m_psm->SBMConcentration(pt, m_Rid);
		double J = ept.m_J;
        double phi0 = pbm->GetReferentialSolidVolumeFraction(pt);
		dcdJ = -ca/(J-phi0);
	}
	else {
		ca = m_psm->GetActualSoluteConcentration(pt, m_Rid);
		double c = m_psm->GetEffectiveSoluteConcentration(pt, m_Rid);
		double dkdJ = m_psm->dkdJ(pt, m_Rid);
		dcdJ = dkdJ*c;
	}
	
	double dzhatdJ = 0;
	if (ca > m_c0) {
        double Vmax = m_pFwd->ReactionRate(pt);
        dzhatdJ = dcdJ*m_Km*Vmax/SQR(m_Km + ca);
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
	if (m_Rtype || (m_Rid != sol)) return 0;

	double dzhatdc = 0;
	double ca = m_psm->GetActualSoluteConcentration(pt, m_Rid);
	if (ca > m_c0) {
        double Vmax = m_pFwd->ReactionRate(pt);
		double k = m_psm->GetPartitionCoefficient(pt, m_Rid);
		double c = m_psm->GetEffectiveSoluteConcentration(pt, m_Rid);
		double dkdc = m_psm->dkdc(pt, m_Rid, m_Rid);
		dzhatdc = m_Km*Vmax/SQR(m_Km + ca)*(k + dkdc*c);
    }
	
	return dzhatdc;
}
