//
//  FEMembraneMassActionReversible.cpp
//  FEBioMix
//
//  Created by Gerard Ateshian on 3/4/18.
//  Copyright © 2018 febio.org. All rights reserved.
//

#include "FEMembraneMassActionReversible.h"
//-----------------------------------------------------------------------------
//! molar supply at material point
double FEMembraneMassActionReversible::FwdReactionSupply(FEMaterialPoint& pt)
{
    FESolutesMaterialPoint& spt = *pt.ExtractData<FESolutesMaterialPoint>();
    
    // get forward reaction rate
    double k = m_pFwd->ReactionRate(pt);
    
    // evaluate the reaction molar supply
    double zhat = k;
    
    // start with contribution from solutes
    const int nsol = (int)spt.m_c.size();
    for (int i=0; i<nsol; ++i) {
        int vR = m_vR[i];
        if (vR > 0) {
            double c = spt.m_c[i];
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
    
    // add contribution from internal and external solutes
    const int nse = (int)spt.m_ce.size();
    for (int i=0; i<nse; ++i) {
        int vRe = m_vRe[spt.m_ide[i]];
        if (vRe > 0) {
            // evaluate nodal effective concentrations
            double c = spt.m_ce[i];
            zhat *= pow(c, vRe);
        }
    }
    const int nsi = (int)spt.m_ci.size();
    for (int i=0; i<nsi; ++i) {
        int vRi = m_vRi[spt.m_idi[i]];
        if (vRi > 0) {
            // evaluate nodal effective concentrations
            double c = spt.m_ci[i];
            zhat *= pow(c, vRi);
        }
    }

    return zhat;
}

//-----------------------------------------------------------------------------
//! molar supply at material point
double FEMembraneMassActionReversible::RevReactionSupply(FEMaterialPoint& pt)
{
    FESolutesMaterialPoint& spt = *pt.ExtractData<FESolutesMaterialPoint>();
    
    // get forward reaction rate
    double k = m_pRev->ReactionRate(pt);
    
    // evaluate the reaction molar supply
    double zhat = k;
    
    // start with contribution from solutes
    const int nsol = (int)spt.m_c.size();
    for (int i=0; i<nsol; ++i) {
        int vP = m_vP[i];
        if (vP > 0) {
            double c = spt.m_c[i];
            zhat *= pow(c, vP);
        }
    }
    
    // add contribution of solid-bound molecules
    const int nsbm = (int)spt.m_sbmr.size();
    for (int i=0; i<nsbm; ++i) {
        int vP = m_vP[nsol+i];
        if (vP > 0) {
            double c = m_pMP->SBMConcentration(pt, i);
            zhat *= pow(c, vP);
        }
    }
    
    // add contribution from internal and external solutes
    const int nse = (int)spt.m_ce.size();
    for (int i=0; i<nse; ++i) {
        int vPe = m_vPe[spt.m_ide[i]];
        if (vPe > 0) {
            // evaluate nodal effective concentrations
            double c = spt.m_ce[i];
            zhat *= pow(c, vPe);
        }
    }
    const int nsi = (int)spt.m_ci.size();
    for (int i=0; i<nsi; ++i) {
        int vPi = m_vPi[spt.m_idi[i]];
        if (vPi > 0) {
            // evaluate nodal effective concentrations
            double c = spt.m_ci[i];
            zhat *= pow(c, vPi);
        }
    }

    return zhat;
}

//-----------------------------------------------------------------------------
//! molar supply at material point
double FEMembraneMassActionReversible::ReactionSupply(FEMaterialPoint& pt)
{
    double zhatF = FwdReactionSupply(pt);
    double zhatR = RevReactionSupply(pt);
    return zhatF - zhatR;
}

//-----------------------------------------------------------------------------
//! tangent of molar supply with strain at material point
double FEMembraneMassActionReversible::Tangent_ReactionSupply_Strain(FEMaterialPoint& pt)
{
    // forward reaction
    double kF = m_pFwd->ReactionRate(pt);
    double dkFde = m_pFwd->Tangent_ReactionRate_Strain(pt);
    double zhatF = FwdReactionSupply(pt);
    double dzhatFde = 0;
    if (kF > 0) dzhatFde = dkFde*(zhatF/kF);
    
    // reverse reaction
    double kR = m_pRev->ReactionRate(pt);
    double dkRde = m_pRev->Tangent_ReactionRate_Strain(pt);
    double zhatR = RevReactionSupply(pt);
    double dzhatRde = 0;
    if (kR > 0) dzhatRde = dkRde*(zhatR/kR);
    
    return dzhatFde - dzhatRde;
}

//-----------------------------------------------------------------------------
//! tangent of molar supply with effective pressure at material point
double FEMembraneMassActionReversible::Tangent_ReactionSupply_Pressure(FEMaterialPoint& pt)
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
//! tangent of molar supply with effective pressure at material point
double FEMembraneMassActionReversible::Tangent_ReactionSupply_Pi(FEMaterialPoint& pt)
{
    return 0;
}

//-----------------------------------------------------------------------------
//! tangent of molar supply with effective pressure at material point
double FEMembraneMassActionReversible::Tangent_ReactionSupply_Pe(FEMaterialPoint& pt)
{
    return 0;
}

//-----------------------------------------------------------------------------
//! tangent of molar supply with effective concentration at material point
double FEMembraneMassActionReversible::Tangent_ReactionSupply_Concentration(FEMaterialPoint& pt, const int sol)
{
    const int nsol = m_nsol;
    
    // if the derivative is taken with respect to a solid-bound molecule, return 0
    if (sol >= nsol) {
        return 0;
    }
    
    FESolutesMaterialPoint& spt = *pt.ExtractData<FESolutesMaterialPoint>();
    
    // forward reaction
    double zhatF = FwdReactionSupply(pt);
    double dzhatFdc = 0;
    if ((zhatF > 0) && (spt.m_c[sol] > 0)) dzhatFdc = m_vR[sol]*zhatF/spt.m_c[sol];
    
    // reverse reaction
    double zhatR = RevReactionSupply(pt);
    double dzhatRdc = 0;
    if ((zhatR > 0) && (spt.m_c[sol] > 0)) dzhatRdc = m_vP[sol]*zhatR/spt.m_c[sol];
    
    return dzhatFdc - dzhatRdc;
}

//-----------------------------------------------------------------------------
//! tangent of molar supply with effective concentration at material point
double FEMembraneMassActionReversible::Tangent_ReactionSupply_Ci(FEMaterialPoint& pt, const int sol)
{
    FESolutesMaterialPoint& spt = *pt.ExtractData<FESolutesMaterialPoint>();
    
    // forward reaction
    double zhatF = FwdReactionSupply(pt);
    double kF = m_pFwd->ReactionRate(pt);
    double dkFdci = m_pFwd->Tangent_ReactionRate_Ci(pt, sol);
    double dzhatFdc = 0;
    if (kF != 0) dzhatFdc = dkFdci/kF*zhatF;
    if ((zhatF > 0) && (spt.m_ci[sol] > 0)) dzhatFdc = m_vRi[spt.m_idi[sol]]*zhatF/spt.m_ci[sol];
    
    // reverse reaction
    double zhatR = RevReactionSupply(pt);
    double kR = m_pRev->ReactionRate(pt);
    double dkRdci = m_pRev->Tangent_ReactionRate_Ci(pt, sol);
    double dzhatRdc = 0;
    if (kR != 0) dzhatRdc = dkRdci/kR*zhatR;
    if ((zhatR > 0) && (spt.m_ci[sol] > 0)) dzhatRdc += m_vPi[spt.m_idi[sol]]*zhatR/spt.m_ci[sol];
    
    return dzhatFdc - dzhatRdc;
}

//-----------------------------------------------------------------------------
//! tangent of molar supply with effective concentration at material point
double FEMembraneMassActionReversible::Tangent_ReactionSupply_Ce(FEMaterialPoint& pt, const int sol)
{
    FESolutesMaterialPoint& spt = *pt.ExtractData<FESolutesMaterialPoint>();
    
    // forward reaction
    double zhatF = FwdReactionSupply(pt);
    double kF = m_pFwd->ReactionRate(pt);
    double dkFdce = m_pFwd->Tangent_ReactionRate_Ce(pt, sol);
    double dzhatFdc = 0;
    if (kF != 0) dzhatFdc = dkFdce/kF*zhatF;
    if ((zhatF > 0) && (spt.m_ce[sol] > 0)) dzhatFdc = m_vRe[spt.m_ide[sol]]*zhatF/spt.m_ce[sol];
    
    // reverse reaction
    double zhatR = RevReactionSupply(pt);
    double kR = m_pRev->ReactionRate(pt);
    double dkRdce = m_pRev->Tangent_ReactionRate_Ce(pt, sol);
    double dzhatRdc = 0;
    if (kR != 0) dzhatRdc = dkRdce/kR*zhatR;
    if ((zhatR > 0) && (spt.m_ce[sol] > 0)) dzhatRdc += m_vPe[spt.m_ide[sol]]*zhatR/spt.m_ce[sol];
    
    return dzhatFdc - dzhatRdc;
}
