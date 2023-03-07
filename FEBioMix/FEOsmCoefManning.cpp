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
#include "FEOsmCoefManning.h"
#include "FEMultiphasic.h"
#include <FECore/log.h>
#include <FECore/FECoreKernel.h>

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FEOsmCoefManning, FEOsmoticCoefficient)
    ADD_PARAMETER(m_ksi , FE_RANGE_GREATER_OR_EQUAL(0.0), "ksi"  );
    ADD_PARAMETER (m_sol , "co_ion")->setEnums("$(solutes)");
    ADD_PROPERTY(m_osmc, "osmc")->SetLongName("Wells osmotic coefficient");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor.
FEOsmCoefManning::FEOsmCoefManning(FEModel* pfem) : FEOsmoticCoefficient(pfem)
{
    m_ksi = 1;
    m_sol = -1;
    m_lsol = -1;
    m_osmc = nullptr;
}

//-----------------------------------------------------------------------------
bool FEOsmCoefManning::Init()
{
    // get the ancestor material which must be a multiphasic material
    FESoluteInterface* psm = dynamic_cast<FESoluteInterface*>(GetAncestor());
	if (psm == nullptr) {
		feLogError("Ancestor material must have solutes");
		return false;
	}
    
    // extract the local id of the solute from the global id
    // m_sol must be zero-based
    m_lsol = psm->FindLocalSoluteID(m_sol);
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
    FESoluteInterface* psm = dynamic_cast<FESoluteInterface*>(GetAncestor());

    // evaluate X = FCD/co-ion actual concentration
    double ca = psm->GetActualSoluteConcentration(mp, m_lsol);
    double cF = fabs(psm->GetFixedChargeDensity(mp));
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
    FESoluteInterface* psm = dynamic_cast<FESoluteInterface*>(GetAncestor());
    FEBiphasicInterface* pbm = dynamic_cast<FEBiphasicInterface*>(GetAncestor());

    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    
    // evaluate X = FCD/co-ion actual concentration
    double ca = psm->GetActualSoluteConcentration(mp, m_lsol);
    double cF = fabs(psm->GetFixedChargeDensity(mp));
    double kt = psm->GetPartitionCoefficient(mp, m_lsol);
    double dktdJ = psm->dkdJ(mp, m_lsol);

    double X = 0;
    if (ca > 0) X = cF/ca;
    
    // evaluate dX/dJ
    double J = pt.m_J;
    
    double phisr = pbm->GetReferentialSolidVolumeFraction(mp);
    
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
    FESoluteInterface* psm = dynamic_cast<FESoluteInterface*>(GetAncestor());

    // evaluate X = FCD/co-ion actual concentration
    double ca = psm->GetActualSoluteConcentration(mp, m_lsol);
    double cF = fabs(psm->GetFixedChargeDensity(mp));
    double kta = psm->GetPartitionCoefficient(mp, m_lsol);
    double kt = psm->GetPartitionCoefficient(mp, isol);
    int zt = psm->GetSolute(isol)->ChargeNumber();
    
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
    FESoluteInterface* psm = dynamic_cast<FESoluteInterface*>(GetAncestor());

    double ca = psm->GetActualSoluteConcentration(mp, m_lsol);
    double osmc = m_osmc->value(ca);
    
    assert(osmc>0);
    
    return osmc;
}

//-----------------------------------------------------------------------------
//! Tangent of osmotic coefficient with respect to strain
double FEOsmCoefManning::Tangent_OsmoticCoefficient_Strain_Wells(FEMaterialPoint &mp)
{
    FESoluteInterface* psm = dynamic_cast<FESoluteInterface*>(GetAncestor());

    double ca = psm->GetActualSoluteConcentration(mp, m_lsol);
    double c = psm->GetEffectiveSoluteConcentration(mp, m_lsol);
    double dkdJ = psm->dkdJ(mp, m_lsol);
    double f = dkdJ*c;
    
    double dosmc = m_osmc->derive(ca);
    
    return dosmc*f;
}

//-----------------------------------------------------------------------------
//! Tangent of osmotic coefficient with respect to concentration
double FEOsmCoefManning::Tangent_OsmoticCoefficient_Concentration_Wells(FEMaterialPoint &mp, const int isol)
{
    FESoluteInterface* psm = dynamic_cast<FESoluteInterface*>(GetAncestor());

    double ca = psm->GetActualSoluteConcentration(mp, m_lsol);
    double c = psm->GetEffectiveSoluteConcentration(mp, m_lsol);
    double k = psm->GetPartitionCoefficient(mp, m_lsol);
    double dkdc = psm->dkdc(mp, m_lsol, isol);
    double f = dkdc*c;
    if (isol == m_lsol) f += k;
    
    double dosmc = m_osmc->derive(ca);
    
    return dosmc*f;
}
