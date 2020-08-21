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
    m_Fusi.resize(n, mat3d(1,0,0,
                           0,1,0,
                           0,0,1));
    m_Fvsi.resize(n, mat3d(1,0,0,
                           0,1,0,
                           0,0,1));
    m_Fp = mat3dd(1);
    m_Ku.resize(n);
    m_Kv.resize(n);
    m_gp.resize(n, 0);
    m_gpp.resize(n, 0);
    m_gc.resize(n, 0);
    m_Rhat = 0;
    m_wi.resize(n, 0);
    m_wy.resize(n,0);
    m_gp.resize(n, 0);
    m_Eyt.resize(n, 0);
    m_Eym.resize(n, 0);
    m_di.resize(n, 0);
    m_dy.resize(n, 0);
    m_d.resize(n, 0);
    m_yld.resize(n, 0);
    m_D = 0.0;
    m_Eit = 0.0;
    m_Eim = 0.0;
    
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
    FEMaterialPoint::Init();
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
    m_Eim = max(m_Eim, m_Eit);
    m_Fp = pt.m_F;
    
    // don't forget to update the base class
    FEMaterialPoint::Update(timeInfo);
}

//-----------------------------------------------------------------------------
void FEReactivePlasticDamageMaterialPoint::Serialize(DumpStream& ar)
{
    FEMaterialPoint::Serialize(ar);
    FEReactivePlasticDamage* prp = dynamic_cast<FEReactivePlasticDamage*>(m_pMat);
    int n = prp ? prp->m_n : 1;
    
    if (ar.IsSaving())
    {
        ar << m_Fp << m_D << m_Eit << m_Eim;
        for (int i=0; i<n; ++i) {
            ar << m_Fusi[i] << m_Fvsi[i] << m_Ku[i] << m_Kv[i] << m_gp[i] << m_gpp[i];
            ar << m_wi[i] << m_wy[i] << m_Eyt[i] << m_Eym[i];
            ar << m_di[i] << m_dy[i] << m_d[i] << m_yld[i];
        }
    }
    else
    {
        ar >> m_Fp >> m_D >> m_Eit >> m_Eim;
        for (int i=0; i<n; ++i) {
            ar >> m_Fusi[i] >> m_Fvsi[i] >> m_Ku[i] >> m_Kv[i] >> m_gp[i] >> m_gpp[i];
            ar >> m_wi[i] >> m_wy[i] >> m_Eyt[i] >> m_Eym[i];
            ar >> m_di[i] >> m_dy[i] >> m_d[i] >> m_yld[i];
        }
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

