/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, The Trustees of Columbia University in 
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

