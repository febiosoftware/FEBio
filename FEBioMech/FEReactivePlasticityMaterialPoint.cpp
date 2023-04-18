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
#include "FEReactivePlasticityMaterialPoint.h"
#include "FEElasticMaterial.h"

////////////////////// PLASTICITY MATERIAL POINT //////////////////////////////
//-----------------------------------------------------------------------------
FEMaterialPointData* FEReactivePlasticityMaterialPoint::Copy()
{
    FEReactivePlasticityMaterialPoint* pt = new FEReactivePlasticityMaterialPoint(*this);
    if (m_pNext) pt->m_pNext = m_pNext->Copy();

	pt->m_Fusi = m_Fusi;
	pt->m_Fvsi = m_Fvsi;
	pt->m_Ku = m_Ku;
	pt->m_Kv = m_Kv;
	pt->m_w = m_w;
	pt->m_gp = m_gp;
	pt->m_gpp = m_gpp;
	pt->m_gc = m_gc;
    pt->m_byld = m_byld;

    return pt;
}

//-----------------------------------------------------------------------------
void FEReactivePlasticityMaterialPoint::Init()
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
        m_w.assign(n,0);
        m_gp.assign(n, 0);
        m_gpp.assign(n, 0);
        m_gc.assign(n, 0);
        m_byld.assign(n,false);
        m_Rhat = 0;

        // don't forget to initialize the base class
        FEMaterialPointData::Init();
    }
}

//-----------------------------------------------------------------------------
void FEReactivePlasticityMaterialPoint::Update(const FETimeInfo& timeInfo)
{
    FEElasticMaterialPoint& pt = *m_pNext->ExtractData<FEElasticMaterialPoint>();
    for (int i=0; i<m_Fusi.size(); ++i) {
        m_Fusi[i] = m_Fvsi[i];
        m_Ku[i] = m_Kv[i];
        m_gc[i] += fabs(m_gp[i] - m_gpp[i]);
        m_gpp[i] = m_gp[i];
        if (m_w[i] > 0) m_byld[i] = true;
    }
    m_Fp = pt.m_F;
    
    // don't forget to update the base class
	FEMaterialPointData::Update(timeInfo);
}

//-----------------------------------------------------------------------------
void FEReactivePlasticityMaterialPoint::Serialize(DumpStream& ar)
{
	FEMaterialPointData::Serialize(ar);

    if (ar.IsSaving())
    {
        ar << m_Rhat << m_Fp;
        ar << m_w << m_Fusi << m_Fvsi << m_Ku << m_Kv << m_gp << m_gpp << m_gc;
        ar << m_byld;
    }
    else
    {
        ar >> m_Rhat >> m_Fp;
        ar >> m_w >> m_Fusi >> m_Fvsi >> m_Ku >> m_Kv >> m_gp >> m_gpp >> m_gc;
        ar >> m_byld;
    }
}

//-----------------------------------------------------------------------------
//! Evaluate net mass fraction of yielded bonds
double FEReactivePlasticityMaterialPoint::YieldedBonds() const
{
    double w = 0;
    for (int i=0; i<m_w.size(); ++i) w += m_w[i];
    return w;
}

