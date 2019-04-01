#include "stdafx.h"
#include "FEMembraneMassActionForward.h"
#include "FECore/DOFS.h"

//-----------------------------------------------------------------------------
//! molar supply at material point
double FEMembraneMassActionForward::ReactionSupply(FEMaterialPoint& pt)
{
    FESolutesMaterialPoint& spt = *pt.ExtractData<FESolutesMaterialPoint>();
    
    // get reaction rate
    double kF = m_pFwd->ReactionRate(pt);
    
    // evaluate the reaction molar supply
    double zhat = kF;
    
    // start with contribution from membrane solutes
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
//! tangent of molar supply with strain at material point
double FEMembraneMassActionForward::Tangent_ReactionSupply_Strain(FEMaterialPoint& pt)
{
    double kF = m_pFwd->ReactionRate(pt);
    double dkFde = m_pFwd->Tangent_ReactionRate_Strain(pt);
    double zhat = ReactionSupply(pt);
    double dzhatde = 0;
    if (kF > 0) dzhatde = dkFde*(zhat/kF);
    
    return dzhatde;
}

//-----------------------------------------------------------------------------
//! tangent of molar supply with effective pressure at material point
double FEMembraneMassActionForward::Tangent_ReactionSupply_Pressure(FEMaterialPoint& pt)
{
    double kF = m_pFwd->ReactionRate(pt);
    double dkFdp = m_pFwd->Tangent_ReactionRate_Pressure(pt);
    double zhat = ReactionSupply(pt);
    double dzhatdp = 0;
    if (kF > 0) dzhatdp = dkFdp*zhat/kF;
    return dzhatdp;
}

//-----------------------------------------------------------------------------
//! tangent of molar supply with effective pressure at material point
double FEMembraneMassActionForward::Tangent_ReactionSupply_Pi(FEMaterialPoint& pt)
{
    double kF = m_pFwd->ReactionRate(pt);
    double dkFdp = m_pFwd->Tangent_ReactionRate_Pi(pt);
    double zhat = ReactionSupply(pt);
    double dzhatdp = 0;
    if (kF > 0) dzhatdp = dkFdp*zhat/kF;
    return dzhatdp;
}

//-----------------------------------------------------------------------------
//! tangent of molar supply with effective pressure at material point
double FEMembraneMassActionForward::Tangent_ReactionSupply_Pe(FEMaterialPoint& pt)
{
    double kF = m_pFwd->ReactionRate(pt);
    double dkFdp = m_pFwd->Tangent_ReactionRate_Pe(pt);
    double zhat = ReactionSupply(pt);
    double dzhatdp = 0;
    if (kF > 0) dzhatdp = dkFdp*zhat/kF;
    return dzhatdp;
}

//-----------------------------------------------------------------------------
//! tangent of molar supply with effective concentration at material point
double FEMembraneMassActionForward::Tangent_ReactionSupply_Concentration(FEMaterialPoint& pt, const int sol)
{
    const int nsol = m_nsol;
    
    // if the derivative is taken with respect to a solid-bound molecule, return 0
    if (sol >= nsol) return 0;
    
    FESolutesMaterialPoint& spt = *pt.ExtractData<FESolutesMaterialPoint>();
    double zhat = ReactionSupply(pt);
    double dzhatdc = 0;
    if ((zhat > 0) && (spt.m_c[sol] > 0)) dzhatdc = m_vR[sol]/spt.m_c[sol]*zhat;
    
    return dzhatdc;
}

//-----------------------------------------------------------------------------
//! tangent of molar supply with effective concentration at material point
double FEMembraneMassActionForward::Tangent_ReactionSupply_Ci(FEMaterialPoint& pt, const int sol)
{
    double zhat = ReactionSupply(pt);
    double kF = m_pFwd->ReactionRate(pt);
    double dkFdci = m_pFwd->Tangent_ReactionRate_Ci(pt, sol);
    double dzhatdc = 0;
    FESolutesMaterialPoint& spt = *pt.ExtractData<FESolutesMaterialPoint>();
    if (kF != 0) dzhatdc = dkFdci/kF*zhat;
    if ((zhat > 0) && (spt.m_ci[sol] > 0)) dzhatdc += m_vRi[spt.m_idi[sol]]/spt.m_ci[sol]*zhat;
    
    return dzhatdc;
}

//-----------------------------------------------------------------------------
//! tangent of molar supply with effective concentration at material point
double FEMembraneMassActionForward::Tangent_ReactionSupply_Ce(FEMaterialPoint& pt, const int sol)
{
    double zhat = ReactionSupply(pt);
    double kF = m_pFwd->ReactionRate(pt);
    double dkFdce = m_pFwd->Tangent_ReactionRate_Ce(pt, sol);
    double dzhatdc = 0;
    FESolutesMaterialPoint& spt = *pt.ExtractData<FESolutesMaterialPoint>();
    if (kF != 0) dzhatdc = dkFdce/kF*zhat;
    if ((zhat > 0) && (spt.m_ce[sol] > 0)) dzhatdc += m_vRe[spt.m_ide[sol]]/spt.m_ce[sol]*zhat;
    
    return dzhatdc;
}
