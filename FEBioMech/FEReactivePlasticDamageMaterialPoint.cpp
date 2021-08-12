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
FEMaterialPoint* FEReactivePlasticDamageMaterialPoint::Copy()
{
    FEReactivePlasticDamageMaterialPoint* pt = new FEReactivePlasticDamageMaterialPoint(*this);
    if (m_pNext) pt->m_pNext = m_pNext->Copy();
    return pt;
}

//-----------------------------------------------------------------------------
void FEReactivePlasticDamageMaterialPoint::Init()
{
    FEReactivePlasticDamage* prp = dynamic_cast<FEReactivePlasticDamage*>(m_pMat);
    
    // get size of vectors
    int n = prp ? prp->m_n : 1;
    
    // intialize data
    m_Fusi.assign(n, mat3d(1,0,0,
                           0,1,0,
                           0,0,1));
    m_Fvsi.assign(n, mat3d(1,0,0,
                           0,1,0,
                           0,0,1));
    m_Fp = mat3dd(1);
    m_Ku.resize(n);
    m_Kv.resize(n);
    m_gp.assign(n, 0);
    m_gpp.assign(n, 0);
    m_gc.assign(n, 0);
    m_Rhat = 0;
    m_wi.assign(n, 0);
    m_wy.assign(n,0);
    m_gp.assign(n, 0);
    m_Eyt.assign(n, 0);
    m_Eym.assign(n, 0);
    m_di.assign(n, 0);
    m_dy.assign(n, 0);
    m_d.assign(n, 0);
    m_yld.assign(n, 0);
    
    // Need to check this now that wmax/bias exist
    
    // initialize initial intact bond fraction values
    if (n == 1) {
        m_wi[0] = m_pMat->m_wmin;
    }
    else {
        m_wi[0] = m_pMat->m_wmin;
        double sw = m_wi[0];
        for (int i=1; i<n; ++i) {
            m_wi[i] = (1 - m_pMat->m_wmin)/(n-1);
            sw += m_wi[i];
        }
    }
    
    if (n == 1) {
        m_wi[0] = m_pMat->m_wmin;
    }
    else {
        // use bias r to reduce intervals in Ky and w as they increase proportionally
        double r = m_pMat->m_bias;
        // r= 1 uses uniform intervals
        if (r == 1) {
            m_wi[0] = m_pMat->m_wmin;
            double sw = m_wi[0];
            for (int i=1; i<n; ++i) {
                m_wi[i] = (m_pMat->m_wmax - m_pMat->m_wmin)/(n-1);
                sw += m_wi[i];
            }
        }
        else {
            double c = (1-r)/(1-pow(r, n-1));
            m_wi[0] = m_pMat->m_wmin;
            m_wi[1] = c*(m_pMat->m_wmax-m_pMat->m_wmin);
            double sw = m_wi[0];
            sw += m_wi[1];
            for (int i=2; i<n; ++i) {
                m_wi[i] = m_wi[i-1]*r;
                sw += m_wi[i];
            }
        }
    }
    
    
    
    
    // don't forget to initialize the base class
    FEDamageMaterialPoint::Init();
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
    FEMaterialPoint::Serialize(ar);
    FEReactivePlasticDamage* prp = dynamic_cast<FEReactivePlasticDamage*>(m_pMat);
    int n = prp ? prp->m_n : 1;
    
    if (ar.IsSaving())
    {
        ar << m_Fp << m_D << m_Etrial << m_Emax;
        ar << m_Fusi << m_Fvsi << m_Ku << m_Kv << m_gp << m_gpp << m_gc;
        ar << m_wi << m_wy << m_Eyt << m_Eym;
        ar << m_di << m_dy << m_d << m_yld;
    }
    else
    {
        ar >> m_Fp >> m_D >> m_Etrial >> m_Emax;
        ar >> m_Fusi >> m_Fvsi >> m_Ku >> m_Kv >> m_gp >> m_gpp >> m_gc;
        ar >> m_wi >> m_wy >> m_Eyt >> m_Eym;
        ar >> m_di >> m_dy >> m_d >> m_yld;
    }
}

//! Evaluate net mass fraction of yielded bonds
double FEReactivePlasticDamageMaterialPoint::YieldedBonds()
{
    double w = 0;
    for (int i=0; i<m_wy.size(); ++i) w += m_wy[i];
    return w;
}

//! Evaluate net mass fraction of intact bonds
double FEReactivePlasticDamageMaterialPoint::IntactBonds()
{
    double w = 0;
    for (int i=0; i<m_wi.size(); ++i) w += m_wy[i] + m_d[i];
    return 1.0-w;
}

