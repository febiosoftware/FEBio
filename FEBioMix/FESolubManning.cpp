//
//  FESolubManning.cpp
//  FEBioMix
//
//  Created by Gerard Ateshian on 1/7/16.
//  Copyright © 2016 febio.org. All rights reserved.
//

#include "FESolubManning.h"
#include "FECore/FEModel.h"

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_PARAMETER_LIST(FESolubManning, FESoluteSolubility)
    ADD_PARAMETER2(m_ksi, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "ksi"  );
    ADD_PARAMETER (m_sol, FE_PARAM_INT, "co_ion");
    ADD_PARAMETER2(m_solub, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "solub");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! Constructor.
FESolubManning::FESolubManning(FEModel* pfem) : FESoluteSolubility(pfem)
{
    m_ksi = 1;
    m_sol = -1;
    m_lsol = -1;
    m_solub = 1;
    m_bcoi = false;
    m_lc = -1;
    m_pMP = nullptr;
}

//-----------------------------------------------------------------------------
bool FESolubManning::Init()
{
    if (FESoluteSolubility::Init() == false) return false;
    
    // get the parent which must be a solute material
    FESolute* m_pSol = dynamic_cast<FESolute*>(GetParent());
    
    // set m_bcoion flag
    if (m_pSol->GetSoluteID() == m_sol)
        m_bcoi = true;
    else
        m_bcoi = false;
    
    // get the ancestor material which must be a multiphasic material
    m_pMP = dynamic_cast<FEMultiphasic*> (GetAncestor());
    if (m_pMP == nullptr) return MaterialError("Ancestor material must be multiphasic");
    
    // extract the local id of the solute from the global id
    m_lsol = m_pMP->FindLocalSoluteID(m_sol-1); // m_sol must be zero-based
    if (m_lsol == -1) return MaterialError("Invalid value for sol");
    
    // extract optional load curve for Wells analysis
    if (m_lc >= 0 && m_pLC == nullptr)
    {
        m_pLC = GetFEModel()->GetLoadCurve(m_lc);
        if (m_pLC == 0) return false;
    }
    
    return true;
}

//-----------------------------------------------------------------------------
//! Solubility
double FESolubManning::Solubility(FEMaterialPoint& mp)
{
    double kPM = Solubility_Manning(mp);
    double kMM = Solubility_Wells(mp);
    
    return kPM*kMM;
}

//-----------------------------------------------------------------------------
//! Tangent of solubility with respect to strain
double FESolubManning::Tangent_Solubility_Strain(FEMaterialPoint &mp)
{
    double kPM = Solubility_Manning(mp);
    double kMM = Solubility_Wells(mp);
    
    double dkPMdJ = Tangent_Solubility_Strain_Manning(mp);
    double dkMMdJ = Tangent_Solubility_Strain_Wells(mp);
    
    return dkPMdJ*kMM + kPM*dkMMdJ;
}

//-----------------------------------------------------------------------------
//! Tangent of solubility with respect to strain
double FESolubManning::Tangent_Solubility_Concentration(FEMaterialPoint& mp, const int isol)
{
    double kPM = Solubility_Manning(mp);
    double kMM = Solubility_Wells(mp);
    
    double dkPMdc = Tangent_Solubility_Concentration_Manning(mp,isol);
    double dkMMdc = Tangent_Solubility_Concentration_Wells(mp,isol);
    
    return dkPMdc*kMM + kPM*dkMMdc;
}

//-----------------------------------------------------------------------------
//! Cross derivative of solubility with respect to strain and concentration
double FESolubManning::Tangent_Solubility_Strain_Concentration(FEMaterialPoint &mp, const int isol)
{
    // assume 0
    return 0;
}

//-----------------------------------------------------------------------------
//! Second derivative of solubility with respect to strain
double FESolubManning::Tangent_Solubility_Strain_Strain(FEMaterialPoint &mp)
{
    // assume 0
    return 0;
}

//-----------------------------------------------------------------------------
//! Second derivative of solubility with respect to concentration
double FESolubManning::Tangent_Solubility_Concentration_Concentration(FEMaterialPoint &mp, const int isol, const int jsol)
{
    // assume 0
    return 0;
}

//-----------------------------------------------------------------------------
//! Solubility
double FESolubManning::Solubility_Manning(FEMaterialPoint& mp)
{
    FESolutesMaterialPoint& spt = *mp.ExtractData<FESolutesMaterialPoint>();
    
    // evaluate X = FCD/co-ion actual concentration
    double ca = spt.m_ca[m_lsol];
    double cF = fabs(spt.m_cF);
    double X = 0;
    if (ca > 0) X = cF/ca;
    
    // --- Manning activity coefficient ---
    double kh;
    if (m_ksi <= 1)
        kh = exp(0.5*m_ksi*X/(X+2));
    else
    {
        double Y = X/m_ksi;
        if (m_bcoi)
            kh = exp(0.5*Y/(Y+2));
        else
            kh = (X+1)/(Y+1)*exp(0.5*Y/(Y+2));
    }
    
    assert(kh>0);
    
    return kh;
}

//-----------------------------------------------------------------------------
//! Tangent of solubility with respect to strain
double FESolubManning::Tangent_Solubility_Strain_Manning(FEMaterialPoint &mp)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    FEBiphasicMaterialPoint& bpt = *mp.ExtractData<FEBiphasicMaterialPoint>();
    FESolutesMaterialPoint& spt = *mp.ExtractData<FESolutesMaterialPoint>();
    
    // evaluate X = FCD/co-ion actual concentration
    double ca = spt.m_ca[m_lsol];
    double cF = fabs(spt.m_cF);
    double X = 0;
    if (ca > 0) X = cF/ca;
    
    // evaluate dX/dJ
    double J = pt.m_J;
    double phisr = bpt.m_phi0;
    double kt = spt.m_k[m_lsol];
    double dktdJ = spt.m_dkdJ[m_lsol];
    double dXdJ = -(1./(J-phisr)+dktdJ/kt)*X;
    
    // evaluate dkhdX
    double dkhdX = 0;
    if (m_ksi <= 1)
        dkhdX = m_ksi*exp(0.5*m_ksi*X/(X+2))/pow(X+2,2);
    else
    {
        double Y = X/m_ksi;
        double isk = 1./m_ksi;
        if (m_bcoi)
            dkhdX = isk*exp(0.5*Y/(Y+2))/pow(Y+2, 2);
        else
            dkhdX = exp(0.5*Y/(Y+2))*(Y*Y*(2-isk)+Y*(5-3*isk)+4-3*isk)/pow((Y+1)*(Y+2), 2);
    }
    
    // evaluate dkhdJ
    double dkhdJ = dkhdX*dXdJ;
    
    return dkhdJ;
}

//-----------------------------------------------------------------------------
//! Tangent of solubility with respect to concentration
double FESolubManning::Tangent_Solubility_Concentration_Manning(FEMaterialPoint &mp, const int isol)
{
    FESolutesMaterialPoint& spt = *mp.ExtractData<FESolutesMaterialPoint>();
    
    // evaluate X = FCD/co-ion actual concentration
    double ca = spt.m_ca[m_lsol];
    double cF = fabs(spt.m_cF);
    double X = 0;
    if (ca > 0) X = cF/ca;
    
    // evaluate dX/dc
    double kta = spt.m_k[m_lsol];
    double kt = spt.m_k[isol];
    int zt = m_pMP->GetSolute(isol)->ChargeNumber();
    double dXdc = -zt*kt/ca;
    if (isol == m_lsol) dXdc -= kta*X/ca;
    
    // evaluate dkhdX
    double dkhdX = 0;
    if (m_ksi <= 1)
        dkhdX = m_ksi*exp(0.5*m_ksi*X/(X+2))/pow(X+2,2);
    else
    {
        double Y = X/m_ksi;
        double isk = 1./m_ksi;
        if (m_bcoi)
            dkhdX = isk*exp(0.5*Y/(Y+2))/pow(Y+2, 2);
        else
            dkhdX = exp(0.5*Y/(Y+2))*(Y*Y*(2-isk)+Y*(5-3*isk)+4-3*isk)/pow((Y+1)*(Y+2), 2);
    }
    
    // evaluate dkhdc
    double dkhdc = dkhdX*dXdc;
    
    return dkhdc;
}

//-----------------------------------------------------------------------------
double FESolubManning::Solubility_Wells(FEMaterialPoint& mp)
{
    FESolutesMaterialPoint& spt = *mp.ExtractData<FESolutesMaterialPoint>();
    
    double ca = spt.m_ca[m_lsol];
    double solub = m_solub;
    
    if (m_pLC != nullptr)
        solub *= m_pLC->Value(ca);
    
    assert(solub>0);
    
    return solub;
}

//-----------------------------------------------------------------------------
double FESolubManning::Tangent_Solubility_Strain_Wells(FEMaterialPoint& mp)
{
    FESolutesMaterialPoint& spt = *mp.ExtractData<FESolutesMaterialPoint>();
    
    double ca = spt.m_ca[m_lsol];
    
    double dsolub = 0;
    
    if (m_pLC != nullptr)
        dsolub = m_solub*m_pLC->Deriv(ca);
    
    double f = spt.m_dkdJ[m_lsol]*spt.m_c[m_lsol];
    dsolub *= f;
    
    return dsolub;
}

//-----------------------------------------------------------------------------
double FESolubManning::Tangent_Solubility_Concentration_Wells(FEMaterialPoint& mp, const int isol)
{
    FESolutesMaterialPoint& spt = *mp.ExtractData<FESolutesMaterialPoint>();
    
    double ca = spt.m_ca[m_lsol];
    
    double dsolub = 0;
    
    if (m_pLC != nullptr)
        dsolub = m_solub*m_pLC->Deriv(ca);
    
    double f = spt.m_dkdc[m_lsol][isol]*spt.m_c[m_lsol];
    if (isol == m_lsol) f += spt.m_k[m_lsol];
    
    dsolub *= f;
    
    return dsolub;
}

//-----------------------------------------------------------------------------
void FESolubManning::Serialize(DumpStream& ar)
{
    if (ar.IsSaving())
    {
        ar << m_solub << m_lc;
    }
    else
    {
        ar >> m_solub >> m_lc;
        m_pLC = ar.GetFEModel().GetLoadCurve(m_lc);
    }
}

//-----------------------------------------------------------------------------
bool FESolubManning::SetParameterAttribute(FEParam& p, const char* szatt, const char* szval)
{
    if (strcmp(p.name(), "solub") == 0)
    {
        if (strcmp(szatt, "lc") == 0)
        {
            m_lc = atoi(szval) - 1;
            return true;
        }
    }
    return false;
}
