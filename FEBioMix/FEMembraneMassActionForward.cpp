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
#include "FEMembraneMassActionForward.h"
#include "FESoluteInterface.h"

BEGIN_FECORE_CLASS(FEMembraneMassActionForward, FEMembraneReaction)
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEMembraneMassActionForward::FEMembraneMassActionForward(FEModel* pfem) : FEMembraneReaction(pfem)
{
    // set material properties
    ADD_PROPERTY(m_pFwd, "forward_rate", FEProperty::Optional);
}

//-----------------------------------------------------------------------------
//! molar supply at material point
double FEMembraneMassActionForward::ReactionSupply(FEMaterialPoint& pt)
{
    // get reaction rate
    double kF = m_pFwd->ReactionRate(pt);
    
    // evaluate the reaction molar supply
    double zhat = kF;
    
    // start with contribution from membrane solutes
    const int nsol = (int) m_psm->Solutes();
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
            double c = m_psm->SBMArealConcentration(pt, i);
            zhat *= pow(c, vR);
        }
    }
    
    // add contribution from internal and external solutes
    const int nse = (int) m_psm->SolutesExternal(pt);
    for (int i=0; i<nse; ++i) {
        int vRe = m_vRe[m_psm->GetSoluteIDExternal(pt,i)];
        if (vRe > 0) {
            // evaluate nodal effective concentrations
            double c = m_psm->GetEffectiveSoluteConcentrationExternal(pt,i);
            zhat *= pow(c, vRe);
        }
    }
    const int nsi = (int) m_psm->SolutesInternal(pt);
    for (int i=0; i<nsi; ++i) {
        int vRi = m_vRi[m_vRi[m_psm->GetSoluteIDInternal(pt,i)]];
        if (vRi > 0) {
            // evaluate nodal effective concentrations
            double c = m_psm->GetEffectiveSoluteConcentrationInternal(pt,i);
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
    
    double zhat = ReactionSupply(pt);
    double dzhatdc = 0;
    double c = m_psm->GetActualSoluteConcentration(pt, sol);
    if ((zhat > 0) && (c > 0)) dzhatdc = m_vR[sol]/c*zhat;
    
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
    if (kF != 0) dzhatdc = dkFdci/kF*zhat;
    double ci = m_psm->GetEffectiveSoluteConcentrationInternal(pt, sol);
    int IDi = m_psm->GetSoluteIDInternal(pt, sol);
    if ((zhat > 0) && (ci > 0)) dzhatdc += m_vRi[IDi]/ci*zhat;
    
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
    if (kF != 0) dzhatdc = dkFdce/kF*zhat;
    double ce = m_psm->GetEffectiveSoluteConcentrationExternal(pt, sol);
    int IDe = m_psm->GetSoluteIDExternal(pt, sol);
    if ((zhat > 0) && (ce > 0)) dzhatdc += m_vRe[IDe]/ce*zhat;
    
    return dzhatdc;
}
