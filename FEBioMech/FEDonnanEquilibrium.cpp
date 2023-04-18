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
#include "FEDonnanEquilibrium.h"
#include <FECore/FEModel.h>
#include <FECore/log.h>

FEMaterialPointData* FEDonnanEquilibriumMaterialPoint::Copy()
{
    FEDonnanEquilibriumMaterialPoint* pt = new FEDonnanEquilibriumMaterialPoint(*this);
    if (m_pNext) pt->m_pNext = m_pNext->Copy();
    return pt;
}

void FEDonnanEquilibriumMaterialPoint::Init()
{
	FEMaterialPointData::Init();
    
    // intialize data to zero
    m_cF = 0;
    m_osm = 0;
    m_p = 0;
    m_bpi = 0;
}

void FEDonnanEquilibriumMaterialPoint::Serialize(DumpStream& ar)
{
    if (ar.IsSaving())
    {
        ar << m_cF << m_osm << m_p << m_bpi;
    }
    else
    {
        ar >> m_cF >> m_osm >> m_p >> m_bpi;
    }
	FEMaterialPointData::Serialize(ar);
}

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FEDonnanEquilibrium, FEElasticMaterial)
	ADD_PARAMETER(m_phiwr, FE_RANGE_LEFT_OPEN(0.0, 1.0), "phiw0");
    ADD_PARAMETER(m_phisr, "phis0");
	ADD_PARAMETER(m_cFr  , "cF0")->setUnits(UNIT_CONCENTRATION);
	ADD_PARAMETER(m_bosm , FE_RANGE_GREATER_OR_EQUAL(0.0), "bosm")->setUnits(UNIT_CONCENTRATION);
    ADD_PARAMETER(m_Phi  , FE_RANGE_GREATER_OR_EQUAL(0.0), "Phi");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEDonnanEquilibrium::FEDonnanEquilibrium(FEModel* pfem) : FEElasticMaterial(pfem) 
{
	m_Rgas = 0; m_Tabs = 0; m_cFr = 0; m_phiwr = -1; m_phisr = -1;
	m_bnew = false; m_binit = false; m_Phi = 1;
}

//-----------------------------------------------------------------------------
// FEDonnanEquilibrium
bool FEDonnanEquilibrium::Init()
{
    if (!m_binit) {
        if (m_phisr >= 0) {
            m_bnew = true;
            m_phiwr = 1 - m_phisr;  // use value at t=0 to initialize
        }
        m_binit = true;
    }

	m_Rgas = GetFEModel()->GetGlobalConstant("R");
	m_Tabs = GetFEModel()->GetGlobalConstant("T");
	
	if (m_Rgas <= 0) { feLogError("A positive universal gas constant R must be defined in Globals section"); return false; }
	if (m_Tabs <= 0) { feLogError("A positive absolute temperature T must be defined in Globals section");   return false; }
	
	return FEElasticMaterial::Init();
}

//-----------------------------------------------------------------------------
mat3ds FEDonnanEquilibrium::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// jacobian
	double J = pt.m_J;
	
	// calculate fixed charge density in current configuration
    double cF;
    if (m_bnew)
        cF = m_phiwr*m_cFr(mp)/(J-m_phisr);
    else
        cF = m_phiwr*m_cFr(mp)/(J-1+m_phiwr);
	
	// calculate osmotic pressure
	double p = m_Rgas*m_Tabs*m_Phi*(sqrt(cF*cF+m_bosm*m_bosm) - m_bosm);
	
	// calculate T = -p*I
	mat3dd I(1.0);	// identity tensor
	mat3ds s = -p*I;
	return s;
}

//-----------------------------------------------------------------------------
tens4ds FEDonnanEquilibrium::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// jacobian
	double J = pt.m_J;

	// calculate fixed charge density in current configuration
    double cF;
    if (m_bnew)
        cF = m_phiwr*m_cFr(mp)/(J-m_phisr);
    else
        cF = m_phiwr*m_cFr(mp)/(J-1+m_phiwr);
	
	// calculate osmotic pressure
	double tosm = sqrt(cF*cF+m_bosm*m_bosm);	// tissue osmolarity
	double p = m_Rgas*m_Tabs*m_Phi*(tosm - m_bosm);	// osmotic pressure
	
	// calculate derivative of osmotic pressure w.r.t. J
    double bpi;
    if (m_bnew)
        bpi = m_Rgas*m_Tabs*m_Phi*J*cF*cF/(J-m_phisr)/tosm;
    else
        bpi = m_Rgas*m_Tabs*m_Phi*J*cF*cF/(J-1+m_phiwr)/tosm;
	
	mat3dd I(1.0);	// Identity
	
	tens4ds IxI = dyad1s(I);
	tens4ds I4  = dyad4s(I);
	
	// calculate tangent osmotic modulus
	tens4ds c = bpi*IxI + p*(2.0*I4 - IxI);
	return c;
}

//-----------------------------------------------------------------------------
double FEDonnanEquilibrium::FixedChargeDensity(FEMaterialPoint& mp)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // jacobian
    double J = pt.m_J;
    
    // calculate fixed charge density in current configuration
    double cF;
    if (m_bnew)
        cF = m_phiwr*m_cFr(mp)/(J-m_phisr);
    else
        cF = m_phiwr*m_cFr(mp)/(J-1+m_phiwr);
    
    return cF;
}

//-----------------------------------------------------------------------------
double FEDonnanEquilibrium::OsmoticPressure(FEMaterialPoint& mp)
{
    double cF = FixedChargeDensity(mp);
    
    // calculate osmotic pressure
    double p = m_Rgas*m_Tabs*m_Phi*(sqrt(cF*cF+m_bosm*m_bosm) - m_bosm);
    
    return p;
}

//-----------------------------------------------------------------------------
double FEDonnanEquilibrium::OsmoticPressureTangent(FEMaterialPoint& mp)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // jacobian
    double J = pt.m_J;
    
    double cF = FixedChargeDensity(mp);
    double tosm = Osmolarity(mp);
    
    // calculate derivative of osmotic pressure w.r.t. J
    double bpi;
    if (m_bnew)
        bpi = m_Rgas*m_Tabs*m_Phi*J*cF*cF/(J-m_phisr)/tosm;
    else
        bpi = m_Rgas*m_Tabs*m_Phi*J*cF*cF/(J-1+m_phiwr)/tosm;
    
    return bpi;
}

//-----------------------------------------------------------------------------
double FEDonnanEquilibrium::Osmolarity(FEMaterialPoint& mp)
{
    double cF = FixedChargeDensity(mp);
    double tosm = sqrt(cF*cF+m_bosm*m_bosm);    // tissue osmolarity
    
    return tosm;
}

//-----------------------------------------------------------------------------
// update Donnan equilibrium material point at each iteration
void FEDonnanEquilibrium::UpdateSpecializedMaterialPoints(FEMaterialPoint& pt, const FETimeInfo& tp)
{
    // get the Donnan equilibrium material point data
    FEDonnanEquilibriumMaterialPoint& pd = *pt.ExtractData<FEDonnanEquilibriumMaterialPoint>();
    
    pd.m_cFr = m_cFr(pt);
    pd.m_cF = FixedChargeDensity(pt);
    pd.m_osm = Osmolarity(pt);
    pd.m_p = OsmoticPressure(pt);
    pd.m_bpi = OsmoticPressureTangent(pt);
}
