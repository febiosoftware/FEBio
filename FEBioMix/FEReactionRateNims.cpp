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
//! Create a shallow copy of the material point data
FEMaterialPoint* FENimsMaterialPoint::Copy()
{
	FENimsMaterialPoint* pt = new FENimsMaterialPoint(*this);
	if (m_pt) pt->m_pt = m_pt->Copy();
	return pt;
}

//-----------------------------------------------------------------------------
//! Initializes material point data.
void FENimsMaterialPoint::Init(bool bflag)
{
	if (bflag)
	{
		// initialize data
        m_lid = -1;
		m_cmax = 0;
	}
	else
	{
        FESolutesMaterialPoint& pt = *m_pt->ExtractData<FESolutesMaterialPoint>();
        if (m_lid != -1) {
            double c = pt.m_ca[m_lid];
            if (c > m_cmax) m_cmax = c;
        }
	}
    
	// don't forget to intialize the nested data
	if (m_pt) m_pt->Init(bflag);
}

//-----------------------------------------------------------------------------
//! Serialize data to the archive
void FENimsMaterialPoint::ShallowCopy(DumpStream& dmp, bool bsave)
{
	if (m_pt) m_pt->ShallowCopy(dmp, bsave);
    
	if (bsave)
	{
		dmp << m_lid << m_cmax;
	}
	else
	{
		dmp >> m_lid >> m_cmax;
	}
}

//-----------------------------------------------------------------------------
//! Serialize data to the archive
void FENimsMaterialPoint::Serialize(DumpFile& ar)
{
	if (m_pt) m_pt->Serialize(ar);
    
	if (ar.IsSaving())
	{
		ar << m_lid << m_cmax;
	}
	else
	{
		ar >> m_lid >> m_cmax;
	}
}

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
    FENimsMaterialPoint& npt = *pt.ExtractData<FENimsMaterialPoint>();
    double cmax = npt.m_cmax;
    
    double k = m_k0;
    double eps = 1.e-9;
    
    // if we are past the release time and got exposed to the solute
    if ((t >= m_trel) && (cmax > eps)) {
        if (cmax < m_cr) k += (m_kr - m_k0)*cmax/m_cc;
        else k = m_kr;

        // use m_lid to turn off checking history of max solute concentration
        npt.m_lid = -1;
    }
    // otherwise
    else {
        // evaluate reaction rate
        double c = spt.m_ca[m_lid];
        if (c < m_cc) k += (m_kc - m_k0)*c/m_cc;
        else k = m_kc;
        
        // use m_lid to turn on checking history of max solute concentration
        npt.m_lid = m_lid;
    }

	return k;
}

//-----------------------------------------------------------------------------
//! tangent of reaction rate with strain at material point
mat3ds FEReactionRateNims::Tangent_ReactionRate_Strain(FEMaterialPoint& pt)
{
    // get the time
	double t = FEMaterialPoint::time;
    
    double dkdJ;
    
    if (t >= m_trel) {
        dkdJ = 0;
    }
    else {
        FESolutesMaterialPoint& spt = *pt.ExtractData<FESolutesMaterialPoint>();
        double c = spt.m_ca[m_lid];
        double dcdJ = spt.m_dkdJ[m_lid]*spt.m_c[m_lid];
        if (c < m_cc) dkdJ = (m_kc - m_k0)*dcdJ/m_cc;
        else dkdJ = 0;
    }
    
	return mat3dd(dkdJ);
}

//-----------------------------------------------------------------------------
//! tangent of reaction rate with effective fluid pressure at material point
double FEReactionRateNims::Tangent_ReactionRate_Pressure(FEMaterialPoint& pt)
{
	return 0;
}

//-----------------------------------------------------------------------------
//! Create material point data for this material
FEMaterialPoint* FEReactionRateNims::CreateMaterialPointData()
{
	return new FENimsMaterialPoint(m_pReact->m_pMP->CreateMaterialPointData());
}
