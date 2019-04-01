#include "stdafx.h"
#include "FEMultiphasicFCD.h"

//-----------------------------------------------------------------------------
FEFCDMaterialPoint::FEFCDMaterialPoint(FEMaterialPoint* ppt) : FESolutesMaterialPoint(ppt)
{
    // initialize to 1 in case user does not specify the value at element level
    // but only at material level.
    m_cFr = 1.0;
}

//-----------------------------------------------------------------------------
void FEFCDMaterialPoint::Init(bool bflag)
{
}

//-----------------------------------------------------------------------------
FEMaterialPoint* FEMultiphasicFCD::CreateMaterialPointData()
{
    return new FEFCDMaterialPoint(new FESolutesMaterialPoint(new FEBiphasicMaterialPoint(m_pSolid->CreateMaterialPointData())));
}

//-----------------------------------------------------------------------------
//! Fixed charge density in current configuration
double FEMultiphasicFCD::FixedChargeDensity(FEMaterialPoint& pt)
{
    FEElasticMaterialPoint& et = *pt.ExtractData<FEElasticMaterialPoint>();
    FEBiphasicMaterialPoint& bt = *pt.ExtractData<FEBiphasicMaterialPoint>();
    FESolutesMaterialPoint& spt = *pt.ExtractData<FESolutesMaterialPoint>();
    FEFCDMaterialPoint& fpt = *pt.ExtractData<FEFCDMaterialPoint>();
    
    // relative volume
    double J = et.m_J;
    double phi0 = bt.m_phi0;
    double ce = 0;
    
    // add contribution from charged solid-bound molecules
    for (int isbm=0; isbm<(int)m_pSBM.size(); ++isbm)
        ce += SBMChargeNumber(isbm)*spt.m_sbmr[isbm]/SBMMolarMass(isbm);
    
    // multiply material cFr with element material point cFr to account
    // for loadcurve associated with material cFr.
    double cF = (m_cFr*fpt.m_cFr*(1-phi0)+ce)/(J-phi0);
    
    return cF;
}

