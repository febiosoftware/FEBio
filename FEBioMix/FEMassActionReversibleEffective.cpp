//
//  FEMassActionReversibleEffective.cpp
//  FEBioMix
//
//  Created by Gerard Ateshian on 8/23/16.
//  Copyright © 2016 febio.org. All rights reserved.
//

#include "FEMassActionReversibleEffective.h"


//-----------------------------------------------------------------------------
//! molar supply at material point
double FEMassActionReversibleEffective::FwdReactionSupply(FEMaterialPoint& pt)
{
    FESolutesMaterialPoint& spt = *pt.ExtractData<FESolutesMaterialPoint>();
    
    // get forward reaction rate
    double k = m_pFwd->ReactionRate(pt);
    
    // evaluate the reaction molar supply
    double zhat = k;
    
    // start with contribution from solutes
    const int nsol = spt.m_c.size();
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
    
    return zhat;
}

//-----------------------------------------------------------------------------
//! molar supply at material point
double FEMassActionReversibleEffective::RevReactionSupply(FEMaterialPoint& pt)
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
    
    return zhat;
}

//-----------------------------------------------------------------------------
//! molar supply at material point
double FEMassActionReversibleEffective::ReactionSupply(FEMaterialPoint& pt)
{
    double zhatF = FwdReactionSupply(pt);
    double zhatR = RevReactionSupply(pt);
    return zhatF - zhatR;
}

//-----------------------------------------------------------------------------
//! tangent of molar supply with strain at material point
mat3ds FEMassActionReversibleEffective::Tangent_ReactionSupply_Strain(FEMaterialPoint& pt)
{
    FEElasticMaterialPoint& ept = *pt.ExtractData<FEElasticMaterialPoint>();
    FEBiphasicMaterialPoint& bpt = *pt.ExtractData<FEBiphasicMaterialPoint>();
    FESolutesMaterialPoint& spt = *pt.ExtractData<FESolutesMaterialPoint>();
    
    const int nsol = m_nsol;
    const int nsbm = (int)m_v.size() - nsol;
    double J = ept.m_J;
    double phi0 = bpt.m_phi0;
    
    // forward reaction
    double kF = m_pFwd->ReactionRate(pt);
    mat3ds dkFde = m_pFwd->Tangent_ReactionRate_Strain(pt);
    double zhatF = FwdReactionSupply(pt);
    mat3ds dzhatFde = mat3dd(0);
    if (kF > 0) dzhatFde = dkFde*(zhatF/kF);
    
    // reverse reaction
    double kR = m_pRev->ReactionRate(pt);
    mat3ds dkRde = m_pRev->Tangent_ReactionRate_Strain(pt);
    double zhatR = RevReactionSupply(pt);
    mat3ds dzhatRde = mat3dd(0);
    if (kR > 0) dzhatRde = dkRde*(zhatR/kR);
    
    return dzhatFde - dzhatRde;
}

//-----------------------------------------------------------------------------
//! tangent of molar supply with effective pressure at material point
double FEMassActionReversibleEffective::Tangent_ReactionSupply_Pressure(FEMaterialPoint& pt)
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
double FEMassActionReversibleEffective::Tangent_ReactionSupply_Concentration(FEMaterialPoint& pt, const int sol)
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
    if (zhatF > 0) dzhatFdc = m_vR[sol]*zhatF/spt.m_c[sol];
    
    // reverse reaction
    double zhatR = RevReactionSupply(pt);
    double dzhatRdc = 0;
    if (zhatR > 0) dzhatRdc = m_vP[sol]*zhatR/spt.m_c[sol];

    return dzhatFdc - dzhatRdc;
}
