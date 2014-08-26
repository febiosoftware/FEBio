//
//  FEReactionRateNims.cpp
//  FEBioMix
//
//  Created by Gerard Ateshian on 8/23/14.
//  Copyright (c) 2014 febio.org. All rights reserved.
//

#include "FEReactionRateNims.h"
#include "FECore/DOFS.h"

// Material parameters for the FEMultiphasic material
BEGIN_PARAMETER_LIST(FEReactionRateNims, FEMaterial)
ADD_PARAMETER(m_sol, FE_PARAM_INT, "sol");
ADD_PARAMETER(m_k0, FE_PARAM_DOUBLE, "k0");
ADD_PARAMETER(m_kc, FE_PARAM_DOUBLE, "kc");
ADD_PARAMETER(m_kr, FE_PARAM_DOUBLE, "kr");
ADD_PARAMETER(m_cc, FE_PARAM_DOUBLE, "cc");
ADD_PARAMETER(m_cr, FE_PARAM_DOUBLE, "cr");
ADD_PARAMETER(m_trel, FE_PARAM_DOUBLE, "trel");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
void FEReactionRateNims::Init()
{
	FEMaterial::Init();
	
	if (m_cc <= 0) throw MaterialError("cc must be strictly positive");
	if (m_cr <= 0) throw MaterialError("cr must be strictly positive");
	if (m_trel < 0) throw MaterialError("trel must be positive");

    // do only once
    if (m_lid == -1) {
        // get number of DOFS
        DOFS& fedofs = *DOFS::GetInstance();
        int MAX_CDOFS = fedofs.GetCDOFS();
        // check validity of sol
        if (m_sol < 1 || m_sol > MAX_CDOFS)
            throw MaterialError("sol value outside of valid range for solutes");
        
        // convert global sol value to local id
        FEMultiphasic* pmp = m_pReact->m_pMP;
        int gid = m_sol - 1;
        m_lid = -1;
        for (int isol=0; isol<pmp->Solutes(); ++isol)
            if (pmp->GetSolute(isol)->GetSoluteID() == gid) {
                m_lid = isol;
                break;
            }
        
        // check validity of local id
        if (m_lid == -1)
            throw MaterialError("sol does not match any solute in multiphasic material");
    }
}

//-----------------------------------------------------------------------------
//! reaction rate at material point
double FEReactionRateNims::ReactionRate(FEMaterialPoint& pt)
{
    // get the time
	double t = FEMaterialPoint::time;
    
    FESolutesMaterialPoint& spt = *pt.ExtractData<FESolutesMaterialPoint>();
    double c = spt.m_ca[m_lid];
    double cmax = max(c,spt.m_crd[m_cmax]);
    
    double k = m_k0;
    
    // if we are past the release time and got exposed to the solute
    if ((m_trel > 0) && (t >= m_trel)) {
        if (cmax < m_cr) k += (m_kr - m_k0)*cmax/m_cr;
        else k = m_kr;
    }
    // otherwise
    else {
        // evaluate reaction rate
        if (cmax < m_cc) k += (m_kc - m_k0)*cmax/m_cc;
        else k = m_kc;
    }

	return k;
}

//-----------------------------------------------------------------------------
//! tangent of reaction rate with strain at material point
mat3ds FEReactionRateNims::Tangent_ReactionRate_Strain(FEMaterialPoint& pt)
{
	return mat3dd(0);
}

//-----------------------------------------------------------------------------
//! tangent of reaction rate with effective fluid pressure at material point
double FEReactionRateNims::Tangent_ReactionRate_Pressure(FEMaterialPoint& pt)
{
	return 0;
}

//-----------------------------------------------------------------------------
//! reset, initialize and update chemical reaction data in the FESolutesMaterialPoint
void FEReactionRateNims::ResetElementData(FEMaterialPoint& mp)
{
    // store the solute maximum concentration in the optional
    // chemical reaction data vector m_crd in the solutes material point
    FESolutesMaterialPoint& spt = *mp.ExtractData<FESolutesMaterialPoint>();
    spt.m_crd.push_back(0);
    m_cmax = (int)spt.m_crd.size() - 1;
}

void FEReactionRateNims::InitializeElementData(FEMaterialPoint& mp)
{
    FESolutesMaterialPoint& pt = *mp.ExtractData<FESolutesMaterialPoint>();
    double c = pt.m_ca[m_lid];
    double cmax = pt.m_crd[m_cmax];
    if (c > cmax) pt.m_crd[m_cmax] = c;
}

void FEReactionRateNims::UpdateElementData(FEMaterialPoint& mp)
{
}
