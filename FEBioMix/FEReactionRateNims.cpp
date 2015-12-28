//
//  FEReactionRateNims.cpp
//  FEBioMix
//
//  Created by Gerard Ateshian on 8/23/14.
//  Copyright (c) 2014 febio.org. All rights reserved.
//

#include "FEReactionRateNims.h"
#include "FECore/FEModel.h"

// Material parameters for the FEMultiphasic material
BEGIN_PARAMETER_LIST(FEReactionRateNims, FEMaterial)
	ADD_PARAMETER(m_sol, FE_PARAM_INT, "sol");
	ADD_PARAMETER(m_k0, FE_PARAM_DOUBLE, "k0");
	ADD_PARAMETER(m_kc, FE_PARAM_DOUBLE, "kc");
	ADD_PARAMETER(m_kr, FE_PARAM_DOUBLE, "kr");
	ADD_PARAMETER2(m_cc  , FE_PARAM_DOUBLE, FE_RANGE_GREATER         (0.0), "cc");
	ADD_PARAMETER2(m_cr  , FE_PARAM_DOUBLE, FE_RANGE_GREATER         (0.0), "cr");
	ADD_PARAMETER2(m_trel, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "trel");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
bool FEReactionRateNims::Init()
{
	if (FEMaterial::Init() == false) return false;
	
    // do only once
    if (m_lid == -1) {
        // get number of DOFS
        DOFS& fedofs = GetFEModel()->GetDOFS();
        int MAX_CDOFS = fedofs.GetVariableSize("concentration");
        // check validity of sol
        if (m_sol < 1 || m_sol > MAX_CDOFS)
            return MaterialError("sol value outside of valid range for solutes");
        
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
            return MaterialError("sol does not match any solute in multiphasic material");
    }

	return true;
}

#ifndef max
#define max(a,b) ((a)>(b)?(a):(b))
#endif

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
