#include "stdafx.h"
#include "FEMassActionForwardEffective.h"

//-----------------------------------------------------------------------------
//! molar supply at material point
double FEMassActionForwardEffective::ReactionSupply(FEMaterialPoint& pt)
{
    FESolutesMaterialPoint& spt = *pt.ExtractData<FESolutesMaterialPoint>();
    
    // get reaction rate
    double kF = m_pFwd->ReactionRate(pt);
    
    // evaluate the reaction molar supply
    double zhat = kF;
    
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
    
    return zhat;
}

//-----------------------------------------------------------------------------
//! tangent of molar supply with strain at material point
mat3ds FEMassActionForwardEffective::Tangent_ReactionSupply_Strain(FEMaterialPoint& pt)
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
    if (kF > 0) dzhatde = dkFde*(zhat/kF);
    
    return dzhatde;
}

//-----------------------------------------------------------------------------
//! tangent of molar supply with effective pressure at material point
double FEMassActionForwardEffective::Tangent_ReactionSupply_Pressure(FEMaterialPoint& pt)
{
    double kF = m_pFwd->ReactionRate(pt);
    double dkFdp = m_pFwd->Tangent_ReactionRate_Pressure(pt);
    double zhat = ReactionSupply(pt);
    double dzhatdp = 0;
    if (kF > 0) dzhatdp = dkFdp*zhat/kF;
    return dzhatdp;
}

//-----------------------------------------------------------------------------
//! tangent of molar supply with effective concentration at material point
double FEMassActionForwardEffective::Tangent_ReactionSupply_Concentration(FEMaterialPoint& pt, const int sol)
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
