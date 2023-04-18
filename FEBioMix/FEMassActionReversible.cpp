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
#include "FEMassActionReversible.h"
#include "FESoluteInterface.h"
#include "FEBiphasic.h"

BEGIN_FECORE_CLASS(FEMassActionReversible, FEChemicalReaction)
    // set material properties
    ADD_PROPERTY(m_pFwd, "forward_rate", FEProperty::Optional);
    ADD_PROPERTY(m_pRev, "reverse_rate", FEProperty::Optional);
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEMassActionReversible::FEMassActionReversible(FEModel* pfem) : FEChemicalReaction(pfem)
{
}

//-----------------------------------------------------------------------------
//! molar supply at material point
double FEMassActionReversible::FwdReactionSupply(FEMaterialPoint& pt)
{
    // get forward reaction rate
    double k = m_pFwd->ReactionRate(pt);
    
    // evaluate the reaction molar supply
    double zhat = k;
    
    // start with contribution from solutes
    int nsol = (int) m_psm->Solutes();
    for (int i=0; i<nsol; ++i) {
        int vR = m_vR[i];
        if (vR > 0) {
            double c = m_psm->GetActualSoluteConcentration(pt, i);
            zhat *= pow(c, vR);
        }
    }

    // add contribution of solid-bound molecules
    const int nsbm = (int) m_psm->SBMs();
    for (int i=0; i<nsbm; ++i) {
        int vR = m_vR[nsol+i];
        if (vR > 0) {
            double c = m_psm->SBMConcentration(pt, i);
            zhat *= pow(c, vR);
        }
    }
    
    return zhat;
}

//-----------------------------------------------------------------------------
//! molar supply at material point
double FEMassActionReversible::RevReactionSupply(FEMaterialPoint& pt)
{
    // get forward reaction rate
    double k = m_pRev->ReactionRate(pt);
    
    // evaluate the reaction molar supply
    double zhat = k;
    
    // start with contribution from solutes
    int nsol = m_psm->Solutes();
    for (int i=0; i<nsol; ++i) {
        int vP = m_vP[i];
        if (vP > 0) {
            double c = m_psm->GetActualSoluteConcentration(pt, i);
            zhat *= pow(c, vP);
        }
    }

    // add contribution of solid-bound molecules
    const int nsbm = m_psm->SBMs();
    for (int i=0; i<nsbm; ++i) {
        int vP = m_vP[nsol+i];
        if (vP > 0) {
            double c = m_psm->SBMConcentration(pt, i);
            zhat *= pow(c, vP);
        }
    }
    
    return zhat;
}

//-----------------------------------------------------------------------------
//! molar supply at material point
double FEMassActionReversible::ReactionSupply(FEMaterialPoint& pt)
{
	double zhatF = FwdReactionSupply(pt);
	double zhatR = RevReactionSupply(pt);
	return zhatF - zhatR;
}

//-----------------------------------------------------------------------------
//! tangent of molar supply with strain at material point
mat3ds FEMassActionReversible::Tangent_ReactionSupply_Strain(FEMaterialPoint& pt)
{
    FEBiphasicInterface* pbm = dynamic_cast<FEBiphasicInterface*>(GetAncestor());
	FEElasticMaterialPoint& ept = *pt.ExtractData<FEElasticMaterialPoint>();
	
	const int nsol = m_nsol;
	const int nsbm = (int)m_v.size() - nsol;
    double J = ept.m_J;
    double phi0 = pbm->GetReferentialSolidVolumeFraction(pt);
	
	// forward reaction
	double kF = m_pFwd->ReactionRate(pt);
	mat3ds dkFde = m_pFwd->Tangent_ReactionRate_Strain(pt);
	double zhatF = FwdReactionSupply(pt);
	mat3ds dzhatFde = mat3dd(0);
	if (kF > 0) {
		dzhatFde += dkFde/kF;
	}
	mat3ds I = mat3dd(1);
	for (int isol=0; isol<nsol; ++isol)
    {
        double dkdJ = m_psm->dkdJ(pt, isol);
        double k = m_psm->GetPartitionCoefficient(pt, isol);
        dzhatFde += I*(m_vR[isol]*dkdJ/k);
    }
	for (int isbm = 0; isbm<nsbm; ++isbm)
		dzhatFde += I*(m_vR[nsol+isbm]/(J-phi0));
	
	dzhatFde *= zhatF;
	
	// reverse reaction
	double kR = m_pRev->ReactionRate(pt);
	mat3ds dkRde = m_pRev->Tangent_ReactionRate_Strain(pt);
	double zhatR = RevReactionSupply(pt);
	mat3ds dzhatRde = mat3dd(0);
	if (kR > 0) {
		dzhatRde += dkRde/kR;
	}
	for (int isol=0; isol<nsol; ++isol)
    {
        double dkdJ = m_psm->dkdJ(pt, isol);
        double k = m_psm->GetPartitionCoefficient(pt, isol);
        dzhatRde += I*(m_vP[isol]*dkdJ/k);
    }
	for (int isbm = 0; isbm<nsbm; ++isbm)
		dzhatRde -= I*(m_vP[nsol+isbm]/(J-phi0));
	
	dzhatRde *= zhatR;
	
	return dzhatFde - dzhatRde;
}

//-----------------------------------------------------------------------------
//! tangent of molar supply with effective pressure at material point
double FEMassActionReversible::Tangent_ReactionSupply_Pressure(FEMaterialPoint& pt)
{
	// forward reaction
	double kF = m_pFwd->ReactionRate(pt);
	double dzhatFdp = 0;
	if (kF > 0) {
		double dkFdp = m_pFwd->Tangent_ReactionRate_Pressure(pt);
		double zhatF = FwdReactionSupply(pt);
		dzhatFdp = dkFdp*zhatF/kF;
	}
	
	// reverse reaction
	double kR = m_pRev->ReactionRate(pt);
	double dzhatRdp = 0;
	if (kR > 0) {
		double dkRdp = m_pRev->Tangent_ReactionRate_Pressure(pt);
		double zhatR = RevReactionSupply(pt);
		dzhatRdp = dkRdp*zhatR/kR;
	}
	
	return dzhatFdp - dzhatRdp;
}

//-----------------------------------------------------------------------------
//! tangent of molar supply with effective concentration at material point
double FEMassActionReversible::Tangent_ReactionSupply_Concentration(FEMaterialPoint& pt, const int sol)
{
    const int nsol = m_nsol;
    
    // if the derivative is taken with respect to a solid-bound molecule, return 0
    if (sol >= nsol) {
        return 0;
    }
    
    // forward reaction
    double zhatF = FwdReactionSupply(pt);
    double dzhatFdc = 0;
    for (int isol=0; isol<nsol; ++isol) {
        double dkdc = m_psm->dkdc(pt, isol, sol);
        double k = m_psm->GetPartitionCoefficient(pt, isol);
        double c = m_psm->GetEffectiveSoluteConcentration(pt, sol);
        dzhatFdc += m_vR[isol]*dkdc/k;
        if ((isol == sol) && (c > 0))
            dzhatFdc += m_vR[isol]/c;
    }
    
    dzhatFdc *= zhatF;
    
    // reverse reaction
    double zhatR = RevReactionSupply(pt);
    double dzhatRdc = 0;
    for (int isol=0; isol<nsol; ++isol) {
        double dkdc = m_psm->dkdc(pt, isol, sol);
        double k = m_psm->GetPartitionCoefficient(pt, isol);
        double c = m_psm->GetEffectiveSoluteConcentration(pt, sol);

        dzhatRdc += m_vP[isol]*dkdc/k;
        if ((isol == sol) && (c > 0))
            dzhatRdc += m_vP[isol]/c;
    }
    
    dzhatRdc *= zhatR;
    
    return dzhatFdc - dzhatRdc;
}
