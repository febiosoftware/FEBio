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
#include "FEReactivePlasticDamageMaterialPoint.h"
#include "FEElasticMaterial.h"

#ifndef max
#define max(a, b) ((a)>(b)?(a):(b))
#endif

////////////////////// PLASTIC DAMAGE MATERIAL POINT //////////////////////////////
//-----------------------------------------------------------------------------
FEMaterialPointData* FEReactivePlasticDamageMaterialPoint::Copy()
{
    FEReactivePlasticDamageMaterialPoint* pt = new FEReactivePlasticDamageMaterialPoint(*this);
    if (m_pNext) pt->m_pNext = m_pNext->Copy();

    pt->m_Fusi = m_Fusi;
    pt->m_Fvsi = m_Fvsi;
    pt->m_Kv = m_Kv;
    pt->m_Ku = m_Ku;
    pt->m_gp = m_gp;
    pt->m_gpp = m_gpp;
    pt->m_gc = m_gc;
    pt->m_wy = m_wy;
    pt->m_Eyt = m_Eyt;
    pt->m_Eym = m_Eym;
    pt->m_di = m_di;
    pt->m_dy = m_dy;
    pt->m_d = m_d;
    pt->m_byld = m_byld;
    pt->m_byldt = m_byldt;

    return pt;
}

//-----------------------------------------------------------------------------
void FEReactivePlasticDamageMaterialPoint::Init()
{
    FEPlasticFlowCurveMaterialPoint& fp = *ExtractData<FEPlasticFlowCurveMaterialPoint>();
    
    if (fp.m_binit) {
        size_t n = fp.m_Ky.size();
        // intialize data
        m_Fusi.assign(n, mat3d(1,0,0,
                               0,1,0,
                               0,0,1));
        m_Fvsi.assign(n, mat3d(1,0,0,
                               0,1,0,
                               0,0,1));
        m_Fp = mat3dd(1);
        m_Ku.assign(n, 0);
        m_Kv.assign(n, 0);
        m_gp.assign(n, 0);
        m_gpp.assign(n, 0);
        m_gc.assign(n, 0);
        m_Rhat = 0;
        m_wy.assign(n,0);
        m_gp.assign(n, 0);
        m_Eyt.assign(n, 0);
        m_Eym.assign(n, 0);
        m_di.assign(n+1, 0);
        m_dy.assign(n, 0);
        m_d.assign(n+1, 0);
        m_byld.assign(n, false);
        m_byldt.assign(n, false);
        
        // don't forget to initialize the base class
        FEDamageMaterialPoint::Init();
    }
}

//-----------------------------------------------------------------------------
void FEReactivePlasticDamageMaterialPoint::Update(const FETimeInfo& timeInfo)
{
    FEElasticMaterialPoint& pt = *m_pNext->ExtractData<FEElasticMaterialPoint>();
    for (int i=0; i<m_Fusi.size(); ++i) {
        m_Fusi[i] = m_Fvsi[i];
        m_Ku[i] = m_Kv[i];
        m_gc[i] += fabs(m_gp[i] - m_gpp[i]);
        m_gpp[i] = m_gp[i];
        if (m_wy[i] > 0) m_byld[i] = true;
        m_Eym[i] = max(m_Eym[i], m_Eyt[i]);
    }
    m_Emax = max(m_Emax, m_Etrial);
    m_Fp = pt.m_F;
    
    // don't forget to update the base class
    FEDamageMaterialPoint::Update(timeInfo);
}

//-----------------------------------------------------------------------------
void FEReactivePlasticDamageMaterialPoint::Serialize(DumpStream& ar)
{
    FEMaterialPointData::Serialize(ar);
    
    if (ar.IsSaving())
    {
        ar << m_Fp << m_D << m_Etrial << m_Emax;
        ar << m_Fusi << m_Fvsi << m_Ku << m_Kv << m_gp << m_gpp << m_gc;
        ar << m_wy << m_Eyt << m_Eym;
        ar << m_di << m_dy << m_d << m_byld << m_byldt;
    }
    else
    {
        ar >> m_Fp >> m_D >> m_Etrial >> m_Emax;
        ar >> m_Fusi >> m_Fvsi >> m_Ku >> m_Kv >> m_gp >> m_gpp >> m_gc;
        ar >> m_wy >> m_Eyt >> m_Eym;
        ar >> m_di >> m_dy >> m_d >> m_byld >> m_byldt;
    }
}

//! Evaluate net mass fraction of yielded bonds
double FEReactivePlasticDamageMaterialPoint::YieldedBonds() const
{
    double w = 0;
    for (int i=0; i<m_wy.size(); ++i) w += m_wy[i];
    return w;
}

//! Evaluate net mass fraction of intact bonds
double FEReactivePlasticDamageMaterialPoint::IntactBonds() const
{
    double w = 0;
    int n = (int) m_wy.size();
    if (n == 0) return 1.0;
    for (int i=0; i<n; ++i) w += m_wy[i] + m_d[i];
    w += m_d[n];
    return 1.0-w;
}

