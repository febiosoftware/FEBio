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
#include "FESolubManning.h"
#include "FEMultiphasic.h"
#include <FECore/log.h>
#include "FEBioFluid/FEFluidSolutes.h"
#include "FEBioFluid/FESolutesMaterial.h"
#include "FEBioFluid/FEMultiphasicFSI.h"

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FESolubManning, FESoluteSolubility)
    ADD_PARAMETER(m_ksi  , FE_RANGE_GREATER_OR_EQUAL(0.0), "ksi"  );
    ADD_PARAMETER(m_sol  , "co_ion");
    ADD_PROPERTY(m_solub, "solub" );
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor.
FESolubManning::FESolubManning(FEModel* pfem) : FESoluteSolubility(pfem)
{
    m_ksi = 1;
    m_sol = -1;
    m_lsol = -1;
    m_bcoi = false;
    m_pMP = nullptr;
    m_pFS = nullptr;
    m_pSM = nullptr;
    m_pMF = nullptr;
	m_solub = nullptr;
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
    m_pFS = dynamic_cast<FEFluidSolutes*> (GetAncestor());
    m_pSM = dynamic_cast<FESolutesMaterial*> (GetAncestor());
    m_pMF = dynamic_cast<FEMultiphasicFSI*> (GetAncestor());
    if ((m_pMP == nullptr) && (m_pFS == nullptr) && (m_pSM == nullptr) && (m_pMF == nullptr)) {
        feLogError("Ancestor material must have solutes");
        return false;
    }
    
    // extract the local id of the solute from the global id
    // m_sol must be zero-based
    if (m_pMP)
        m_lsol = m_pMP->FindLocalSoluteID(m_sol);
    if (m_pFS)
        m_lsol = m_pFS->FindLocalSoluteID(m_sol);
    if (m_pSM)
        m_lsol = m_pSM->FindLocalSoluteID(m_sol);
    if (m_pMF)
        m_lsol = m_pMF->FindLocalSoluteID(m_sol);
    if (m_lsol == -1) {
        feLogError("Invalid value for sol");
        return false;
    }

	if (m_solub == nullptr) {
		feLogError("Function for solub not assigned");
		return false;
	}
    m_solub->Init();
    
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
    FEFluidSolutesMaterialPoint& fspt = *mp.ExtractData<FEFluidSolutesMaterialPoint>();
    FESolutesMaterial::Point& smpt = *mp.ExtractData<FESolutesMaterial::Point>();
    FEMultiphasicFSIMaterialPoint& mfpt = *mp.ExtractData<FEMultiphasicFSIMaterialPoint>();
    
    // evaluate X = FCD/co-ion actual concentration
    double ca = 0;
    double cF = 0;
    if (m_pMP)
    {
        ca = spt.m_ca[m_lsol];
        cF = fabs(spt.m_cF);
    }
    else if (m_pFS)
    {
        ca = fspt.m_ca[m_lsol];
    }
    else if (m_pSM)
    {
        ca = smpt.m_ca[m_lsol];
    }
    else if (m_pMF)
    {
        ca = mfpt.m_ca[m_lsol];
        cF = fabs(mfpt.m_cF);
    }
    
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
    FEBiphasicFSIMaterialPoint& bfpt = *mp.ExtractData<FEBiphasicFSIMaterialPoint>();
    FESolutesMaterialPoint& spt = *mp.ExtractData<FESolutesMaterialPoint>();
    FEFluidSolutesMaterialPoint& fspt = *mp.ExtractData<FEFluidSolutesMaterialPoint>();
    FESolutesMaterial::Point& smpt = *mp.ExtractData<FESolutesMaterial::Point>();
    FEMultiphasicFSIMaterialPoint& mfpt = *mp.ExtractData<FEMultiphasicFSIMaterialPoint>();
    
    // evaluate X = FCD/co-ion actual concentration
    double ca = 0.0;
    double cF = 0.0;
    double kt = 0.0;
    double dktdJ = 0.0;
    if (m_pMP)
    {
        ca = spt.m_ca[m_lsol];
        cF = fabs(spt.m_cF);
        kt = spt.m_k[m_lsol];
        dktdJ = spt.m_dkdJ[m_lsol];
    }
    else if (m_pFS)
    {
        ca = fspt.m_ca[m_lsol];
        kt = fspt.m_k[m_lsol];
        dktdJ = fspt.m_dkdJ[m_lsol];
    }
    else if (m_pSM)
    {
        ca = smpt.m_ca[m_lsol];
        kt = smpt.m_k[m_lsol];
    }
    else if (m_pMF)
    {
        ca = mfpt.m_ca[m_lsol];
        cF = fabs(mfpt.m_cF);
        kt = mfpt.m_k[m_lsol];
        dktdJ = mfpt.m_dkdJ[m_lsol];
    }
    double X = 0;
    if (ca > 0) X = cF/ca;
    
    // evaluate dX/dJ
    double J = pt.m_J;
    double phisr = 0.0;
    if (m_pMP)
        phisr = bpt.m_phi0;
    if (m_pMF)
        phisr = bfpt.m_phi0;
    
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
    FEFluidSolutesMaterialPoint& fspt = *mp.ExtractData<FEFluidSolutesMaterialPoint>();
    FESolutesMaterial::Point& smpt = *mp.ExtractData<FESolutesMaterial::Point>();
    FEMultiphasicFSIMaterialPoint& mfpt = *mp.ExtractData<FEMultiphasicFSIMaterialPoint>();
    
    // evaluate X = FCD/co-ion actual concentration
    double ca = 0;
    double cF = 0;
    double kta = 0;
    double kt = 0;
    int zt = 0;
    if (m_pMP)
    {
        ca = spt.m_ca[m_lsol];
        cF = fabs(spt.m_cF);
        kta = spt.m_k[m_lsol];
        kt = spt.m_k[isol];
        zt = m_pMP->GetSolute(isol)->ChargeNumber();
    }
    else if (m_pFS)
    {
        ca = fspt.m_ca[m_lsol];
        kta = fspt.m_k[m_lsol];
        kt = fspt.m_k[isol];
        zt = m_pFS->GetSolute(isol)->ChargeNumber();
    }
    else if (m_pSM)
    {
        ca = smpt.m_ca[m_lsol];
        kta = smpt.m_k[m_lsol];
        kt = smpt.m_k[isol];
        zt = m_pSM->GetSolute(isol)->ChargeNumber();
    }
    else if (m_pMF)
    {
        ca = mfpt.m_ca[m_lsol];
        cF = fabs(mfpt.m_cF);
        kta = mfpt.m_k[m_lsol];
        kt = mfpt.m_k[isol];
        zt = m_pMF->GetSolute(isol)->ChargeNumber();
    }
    double X = 0;
    if (ca > 0) X = cF/ca;
    
    // evaluate dX/dc
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
    FEFluidSolutesMaterialPoint& fspt = *mp.ExtractData<FEFluidSolutesMaterialPoint>();
    FESolutesMaterial::Point& smpt = *mp.ExtractData<FESolutesMaterial::Point>();
    FEMultiphasicFSIMaterialPoint& mfpt = *mp.ExtractData<FEMultiphasicFSIMaterialPoint>();
    
    double ca = 0;
    if (m_pMP)
    {
        ca = spt.m_ca[m_lsol];
    }
    else if (m_pFS)
    {
        ca = fspt.m_ca[m_lsol];
    }
    else if (m_pSM)
    {
        ca = smpt.m_ca[m_lsol];
    }
    else if (m_pMF)
    {
        ca = mfpt.m_ca[m_lsol];
    }
    
    double solub = m_solub->value(ca);
    
    assert(solub>0);
    
    return solub;
}

//-----------------------------------------------------------------------------
double FESolubManning::Tangent_Solubility_Strain_Wells(FEMaterialPoint& mp)
{
    FESolutesMaterialPoint& spt = *mp.ExtractData<FESolutesMaterialPoint>();
    FEFluidSolutesMaterialPoint& fspt = *mp.ExtractData<FEFluidSolutesMaterialPoint>();
    FESolutesMaterial::Point& smpt = *mp.ExtractData<FESolutesMaterial::Point>();
    FEMultiphasicFSIMaterialPoint& mfpt = *mp.ExtractData<FEMultiphasicFSIMaterialPoint>();
    
    double f = 0;
    double ca = 0;
    if (m_pMP)
    {
        ca = spt.m_ca[m_lsol];
        f = spt.m_dkdJ[m_lsol]*spt.m_c[m_lsol];
    }
    else if (m_pFS)
    {
        ca = fspt.m_ca[m_lsol];
        f = fspt.m_dkdJ[m_lsol]*fspt.m_c[m_lsol];
    }
    else if (m_pSM)
    {
        ca = smpt.m_ca[m_lsol];
        f = smpt.m_dkdJ[m_lsol]*smpt.m_c[m_lsol];
    }
    else if (m_pMF)
    {
        ca = mfpt.m_ca[m_lsol];
        f = mfpt.m_dkdJ[m_lsol]*mfpt.m_c[m_lsol];
    }
    
    double dsolub = m_solub->derive(ca);
    
    dsolub *= f;
    
    return dsolub;
}

//-----------------------------------------------------------------------------
double FESolubManning::Tangent_Solubility_Concentration_Wells(FEMaterialPoint& mp, const int isol)
{
    FESolutesMaterialPoint& spt = *mp.ExtractData<FESolutesMaterialPoint>();
    FEFluidSolutesMaterialPoint& fspt = *mp.ExtractData<FEFluidSolutesMaterialPoint>();
    FESolutesMaterial::Point& smpt = *mp.ExtractData<FESolutesMaterial::Point>();
    FEMultiphasicFSIMaterialPoint& mfpt = *mp.ExtractData<FEMultiphasicFSIMaterialPoint>();
    
    double ca = 0;
    double f = 0;
    if (m_pMP)
    {
        ca = spt.m_ca[m_lsol];
        f = spt.m_dkdc[m_lsol][isol]*spt.m_c[m_lsol];
        if (isol == m_lsol) f += spt.m_k[m_lsol];
    }
    else if (m_pFS)
    {
        ca = fspt.m_ca[m_lsol];
        f = fspt.m_dkdc[m_lsol][isol]*fspt.m_c[m_lsol];
        if (isol == m_lsol) f += fspt.m_k[m_lsol];
    }
    else if (m_pSM)
    {
        ca = smpt.m_ca[m_lsol];
        f = smpt.m_dkdc[m_lsol][isol]*smpt.m_c[m_lsol];
        if (isol == m_lsol) f += smpt.m_k[m_lsol];
    }
    else if (m_pMF)
    {
        ca = mfpt.m_ca[m_lsol];
        f = mfpt.m_dkdc[m_lsol][isol]*mfpt.m_c[m_lsol];
        if (isol == m_lsol) f += mfpt.m_k[m_lsol];
    }
    
    double dsolub = m_solub->derive(ca);
    
    
    dsolub *= f;
    
    return dsolub;
}
