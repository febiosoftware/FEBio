/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2020 University of Utah, The Trustees of Columbia University in 
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
#include "FEMassActionForward.h"

//-----------------------------------------------------------------------------
//! molar supply at material point
double FEMassActionForward::ReactionSupply(FEMaterialPoint& pt)
{
    FESolutesMaterialPoint& spt = *pt.ExtractData<FESolutesMaterialPoint>();
    FEFluidSolutesMaterialPoint& fspt = *pt.ExtractData<FEFluidSolutesMaterialPoint>();
    FESolutesMaterial::Point& smpt = *pt.ExtractData<FESolutesMaterial::Point>();
    FEMultiphasicFSIMaterialPoint& mfpt = *pt.ExtractData<FEMultiphasicFSIMaterialPoint>();
    
    // get reaction rate
    double kF = m_pFwd->ReactionRate(pt);
    
    // evaluate the reaction molar supply
    double zhat = kF;
    
    // start with contribution from solutes
    int nsol = 0;
    if (m_pMP)
        nsol = (int)spt.m_ca.size();
    else if (m_pFS)
        nsol = (int)fspt.m_ca.size();
    else if (m_pSM)
        nsol = (int)smpt.m_ca.size();
    else if (m_pMF)
        nsol = (int)mfpt.m_ca.size();
    for (int i=0; i<nsol; ++i) {
        int vR = m_vR[i];
        if (vR > 0) {
            double c = 0;
            if (m_pMP)
                c = spt.m_ca[i];
            else if (m_pFS)
                c = fspt.m_ca[i];
            else if (m_pSM)
                c = smpt.m_ca[i];
            else if (m_pMF)
                c = mfpt.m_ca[i];
            zhat *= pow(c, vR);
        }
    }
    
    // add contribution of solid-bound molecules
    if (m_pMP)
    {
        const int nsbm = (int)spt.m_sbmr.size();
        for (int i=0; i<nsbm; ++i) {
            int vR = m_vR[nsol+i];
            if (vR > 0) {
                double c = m_pMP->SBMConcentration(pt, i);
                zhat *= pow(c, vR);
            }
        }
    }
    
    return zhat;
}

//-----------------------------------------------------------------------------
//! tangent of molar supply with strain at material point
mat3ds FEMassActionForward::Tangent_ReactionSupply_Strain(FEMaterialPoint& pt)
{
	FEElasticMaterialPoint& ept = *pt.ExtractData<FEElasticMaterialPoint>();
	FEBiphasicMaterialPoint& bpt = *pt.ExtractData<FEBiphasicMaterialPoint>();
	FESolutesMaterialPoint& spt = *pt.ExtractData<FESolutesMaterialPoint>();
    FEMultiphasicFSIMaterialPoint& mfpt = *pt.ExtractData<FEMultiphasicFSIMaterialPoint>();
    FEBiphasicFSIMaterialPoint& bfpt = *pt.ExtractData<FEBiphasicFSIMaterialPoint>();
	
	const int nsol = m_nsol;
	const int nsbm = (int)m_v.size() - nsol;
    
    if (m_pMF)
    {
        double J = ept.m_J;
        double phi0 = bfpt.m_phi0;
        
        double kF = m_pFwd->ReactionRate(pt);
        mat3ds dkFde = m_pFwd->Tangent_ReactionRate_Strain(pt);
        double zhat = ReactionSupply(pt);
        mat3ds dzhatde = mat3dd(0);
        if (kF > 0) {
            dzhatde += dkFde/kF;
        }
        mat3ds I = mat3dd(1);
        for (int isol=0; isol<nsol; ++isol)
            dzhatde += I*(m_vR[isol]*mfpt.m_dkdJ[isol]/mfpt.m_k[isol]);
        for (int isbm = 0; isbm<nsbm; ++isbm)
            dzhatde -= I*(m_vR[nsol+isbm]/(J-phi0));
        
        dzhatde *= zhat;
        
        return dzhatde;
    }
    else
    {
        double J = ept.m_J;
        double phi0 = bpt.m_phi0;
        
        double kF = m_pFwd->ReactionRate(pt);
        mat3ds dkFde = m_pFwd->Tangent_ReactionRate_Strain(pt);
        double zhat = ReactionSupply(pt);
        mat3ds dzhatde = mat3dd(0);
        if (kF > 0) {
            dzhatde += dkFde/kF;
        }
        mat3ds I = mat3dd(1);
        for (int isol=0; isol<nsol; ++isol)
            dzhatde += I*(m_vR[isol]*spt.m_dkdJ[isol]/spt.m_k[isol]);
        for (int isbm = 0; isbm<nsbm; ++isbm)
            dzhatde -= I*(m_vR[nsol+isbm]/(J-phi0));
        
        dzhatde *= zhat;
        
        return dzhatde;
    }
	
}

//-----------------------------------------------------------------------------
//! tangent of molar supply with effective pressure at material point
double FEMassActionForward::Tangent_ReactionSupply_Pressure(FEMaterialPoint& pt)
{
	double kF = m_pFwd->ReactionRate(pt);
	double dkFdp = m_pFwd->Tangent_ReactionRate_Pressure(pt);
	double zhat = ReactionSupply(pt);
	double dzhatdp = 0;
	if (kF > 0) {
		dzhatdp = dkFdp*zhat/kF;
	}
	return dzhatdp;
}

//-----------------------------------------------------------------------------
//! tangent of molar supply with effective concentration at material point
double FEMassActionForward::Tangent_ReactionSupply_Concentration(FEMaterialPoint& pt, const int sol)
{
    const int nsol = m_nsol;
    
    // if the derivative is taken with respect to a solid-bound molecule, return 0
    if (sol >= nsol) {
        return 0;
    }
    
    FESolutesMaterialPoint& spt = *pt.ExtractData<FESolutesMaterialPoint>();
    FEFluidSolutesMaterialPoint& fspt = *pt.ExtractData<FEFluidSolutesMaterialPoint>();
    FESolutesMaterial::Point& smpt = *pt.ExtractData<FESolutesMaterial::Point>();
    FEMultiphasicFSIMaterialPoint& mfpt = *pt.ExtractData<FEMultiphasicFSIMaterialPoint>();
    double zhat = ReactionSupply(pt);
    double dzhatdc = 0;
    for (int isol=0; isol<nsol; ++isol) {
        if (m_pMP)
        {
            dzhatdc += m_vR[isol]*spt.m_dkdc[isol][sol]/spt.m_k[isol];
            if ((isol == sol) && (spt.m_c[sol] > 0))
                dzhatdc += m_vR[isol]/spt.m_c[sol];
        }
        else if (m_pFS)
        {
            dzhatdc += m_vR[isol]*fspt.m_dkdc[isol][sol]/fspt.m_k[isol];
            if ((isol == sol) && (fspt.m_c[sol] > 0))
                dzhatdc += m_vR[isol]/fspt.m_c[sol];
        }
        else if (m_pSM)
        {
            dzhatdc += m_vR[isol]*smpt.m_dkdc[isol][sol]/smpt.m_k[isol];
            if ((isol == sol) && (smpt.m_c[sol] > 0))
                dzhatdc += m_vR[isol]/smpt.m_c[sol];
        }
        else if (m_pMF)
        {
            dzhatdc += m_vR[isol]*mfpt.m_dkdc[isol][sol]/mfpt.m_k[isol];
            if ((isol == sol) && (mfpt.m_c[sol] > 0))
                dzhatdc += m_vR[isol]/mfpt.m_c[sol];
        }
    }
    
    dzhatdc *= zhat;
    
    return dzhatdc;
}
