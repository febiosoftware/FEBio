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
#include "FEOsmCoefManning.h"
#include "FEMultiphasic.h"
#include "FEBioFluid/FEFluidSolutes.h"
#include "FEBioFluid/FESolutesMaterial.h"
#include "FEBioFluid/FEMultiphasicFSI.h"
#include <FECore/log.h>

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FEOsmCoefManning, FEOsmoticCoefficient)
    ADD_PARAMETER(m_ksi , FE_RANGE_GREATER_OR_EQUAL(0.0), "ksi"  );
    ADD_PARAMETER (m_sol , "co_ion");
    ADD_PROPERTY(m_osmc, "osmc"  );
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor.
FEOsmCoefManning::FEOsmCoefManning(FEModel* pfem) : FEOsmoticCoefficient(pfem)
{
    m_ksi = 1;
    m_sol = -1;
    m_lsol = -1;
    m_pMP = nullptr;
    m_pFS = nullptr;
    m_pSM = nullptr;
    m_pMF = nullptr;
    // by default, set the osmotic coefficient to a constant = 1
    m_osmc = new FELinearFunction(pfem,0,1);
}

//-----------------------------------------------------------------------------
bool FEOsmCoefManning::Init()
{
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

	if (m_osmc == nullptr) {
		feLogError("function for osmc not specified");
		return false;
	}
    
    return FEOsmoticCoefficient::Init();
}

//-----------------------------------------------------------------------------
//! Osmotic coefficient
double FEOsmCoefManning::OsmoticCoefficient(FEMaterialPoint& mp)
{
    double phiPM = OsmoticCoefficient_Manning(mp);
    double phiMM = OsmoticCoefficient_Wells(mp);
    
    return phiPM + phiMM - 1;
}

//-----------------------------------------------------------------------------
//! Tangent of osmotic coefficient with respect to strain
double FEOsmCoefManning::Tangent_OsmoticCoefficient_Strain(FEMaterialPoint &mp)
{
    double dphiPMdJ = Tangent_OsmoticCoefficient_Strain_Manning(mp);
    double dphiMMdJ = Tangent_OsmoticCoefficient_Strain_Wells(mp);
    
    return dphiPMdJ + dphiMMdJ;
}

//-----------------------------------------------------------------------------
//! Tangent of osmotic coefficient with respect to concentration
double FEOsmCoefManning::Tangent_OsmoticCoefficient_Concentration(FEMaterialPoint &mp, const int isol)
{
    double dphiPMdc = Tangent_OsmoticCoefficient_Concentration_Manning(mp,isol);
    double dphiMMdc = Tangent_OsmoticCoefficient_Concentration_Wells(mp,isol);
    
    return dphiPMdc + dphiMMdc;
}

//-----------------------------------------------------------------------------
//! Osmotic coefficient
double FEOsmCoefManning::OsmoticCoefficient_Manning(FEMaterialPoint& mp)
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
    
    // --- Manning osmotic coefficient ---
    double osmcoef;
    if (m_ksi <= 1)
        osmcoef = 1 - 0.5*m_ksi*X/(X+2);
    else
        osmcoef = (0.5*X/m_ksi+2)/(X+2);

    assert(osmcoef>0);
    
    return osmcoef;
}

//-----------------------------------------------------------------------------
//! Tangent of osmotic coefficient with respect to strain
double FEOsmCoefManning::Tangent_OsmoticCoefficient_Strain_Manning(FEMaterialPoint &mp)
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
    
    double dosmdX;
    if (m_ksi <= 1)
        dosmdX = -m_ksi/pow(X+2, 2);
    else
        dosmdX = (1-2*m_ksi)/m_ksi/pow(X+2, 2);
    
    return dosmdX*dXdJ;
}

//-----------------------------------------------------------------------------
//! Tangent of osmotic coefficient with respect to concentration
double FEOsmCoefManning::Tangent_OsmoticCoefficient_Concentration_Manning(FEMaterialPoint &mp, const int isol)
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
    
    double dosmdX;
    if (m_ksi <= 1)
        dosmdX = -m_ksi/pow(X+2, 2);
    else
        dosmdX = (1./m_ksi-2)/pow(X+2, 2);
    
    return dosmdX*dXdc;
}

//-----------------------------------------------------------------------------
//! Osmotic coefficient
double FEOsmCoefManning::OsmoticCoefficient_Wells(FEMaterialPoint& mp)
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
    double osmc = m_osmc->value(ca);
    
    assert(osmc>0);
    
    return osmc;
}

//-----------------------------------------------------------------------------
//! Tangent of osmotic coefficient with respect to strain
double FEOsmCoefManning::Tangent_OsmoticCoefficient_Strain_Wells(FEMaterialPoint &mp)
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
    
    double dosmc = m_osmc->derive(ca);
    
    return dosmc*f;
}

//-----------------------------------------------------------------------------
//! Tangent of osmotic coefficient with respect to concentration
double FEOsmCoefManning::Tangent_OsmoticCoefficient_Concentration_Wells(FEMaterialPoint &mp, const int isol)
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
    
    double dosmc = m_osmc->derive(ca);
    
    return dosmc*f;
}
