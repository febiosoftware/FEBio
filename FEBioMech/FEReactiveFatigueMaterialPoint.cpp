/*This file is part of the FEBio source code and is licensed under the MIT license
 listed below.
 
 See Copyright-FEBio.txt for details.
 
 Copyright (c) 2022 University of Utah, The Trustees of Columbia University in
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



#include "FEReactiveFatigueMaterialPoint.h"
#include <FECore/DumpStream.h>

////////////////////// FATIGUE BOND /////////////////////////////////
FatigueBond::FatigueBond()
{
    m_wft = m_wfp = m_Xfmax = m_Xftrl = m_Fft = m_Ffp = m_time = 0;
    m_erase = false;
}

FatigueBond::FatigueBond(const FatigueBond& fb)
{
    m_wft = fb.m_wft;
    m_wfp = fb.m_wfp;
    m_Xfmax = fb.m_Xfmax;
    m_Xftrl = fb.m_Xftrl;
    m_Fft = fb.m_Fft;
    m_Ffp = fb.m_Ffp;
    m_time = fb.m_time;
    m_erase = fb.m_erase;
}

FatigueBond::FatigueBond(FatigueBond& fb)
{
    m_wft = fb.m_wft;
    m_wfp = fb.m_wfp;
    m_Xfmax = fb.m_Xfmax;
    m_Xftrl = fb.m_Xftrl;
    m_Fft = fb.m_Fft;
    m_Ffp = fb.m_Ffp;
    m_time = fb.m_time;
    m_erase = fb.m_erase;
}

void FatigueBond::Update()
{
    m_wfp = m_wft;
    m_Ffp = m_Fft;
    if (m_Xftrl > m_Xfmax) m_Xfmax = m_Xftrl;
}

////////////////////// FATIGUE MATERIAL POINT /////////////////////////////////
//-----------------------------------------------------------------------------
// default constructor
FEReactiveFatigueMaterialPoint::FEReactiveFatigueMaterialPoint(FEMaterialPointData*pt) : FEDamageMaterialPoint(pt)
{
	// TODO: Is this working properly? I don't see how this ExtractData could work?
	 
    // get the reactive fatigue material point data
    FEReactiveFatigueMaterialPoint& rfmp = *pt->ExtractData<FEReactiveFatigueMaterialPoint>();
    
    m_D = rfmp.m_D;
    m_wit = rfmp.m_wit;
    m_wip = rfmp.m_wip;
    m_Ximax = rfmp.m_Ximax;
    m_Xitrl = rfmp.m_Xitrl;
    m_Xip = rfmp.m_Xip;
    m_aXit = rfmp.m_aXit;
    m_aXip = rfmp.m_aXip;
    m_Fit = rfmp.m_Fit;
    m_Fip = rfmp.m_Fip;
    m_wbt = rfmp.m_wbt;
    m_wbp = rfmp.m_wbp;
    m_wft = rfmp.m_wft;
    
    m_fb.clear();
    for (int ig=0; ig<m_fb.size(); ++ig)
        m_fb.push_back(rfmp.m_fb[ig]);
}

//-----------------------------------------------------------------------------
// copy constructor
FEReactiveFatigueMaterialPoint::FEReactiveFatigueMaterialPoint(const FEReactiveFatigueMaterialPoint& rfmp) : FEDamageMaterialPoint(rfmp)
{
    m_D = rfmp.m_D;
    m_wit = rfmp.m_wit;
    m_wip = rfmp.m_wip;
    m_Ximax = rfmp.m_Ximax;
    m_Xitrl = rfmp.m_Xitrl;
    m_Xip = rfmp.m_Xip;
    m_aXit = rfmp.m_aXit;
    m_aXip = rfmp.m_aXip;
    m_Fit = rfmp.m_Fit;
    m_Fip = rfmp.m_Fip;
    m_wbt = rfmp.m_wbt;
    m_wbp = rfmp.m_wbp;
    m_wft = rfmp.m_wft;
    
    m_fb.clear();
    for (int ig=0; ig<m_fb.size(); ++ig)
        m_fb.push_back(rfmp.m_fb[ig]);
}

FEReactiveFatigueMaterialPoint::FEReactiveFatigueMaterialPoint(FEReactiveFatigueMaterialPoint& rfmp) : FEDamageMaterialPoint(rfmp)
{
    m_D = rfmp.m_D;
    m_wit = rfmp.m_wit;
    m_wip = rfmp.m_wip;
    m_Ximax = rfmp.m_Ximax;
    m_Xitrl = rfmp.m_Xitrl;
    m_Xip = rfmp.m_Xip;
    m_aXit = rfmp.m_aXit;
    m_aXip = rfmp.m_aXip;
    m_Fit = rfmp.m_Fit;
    m_Fip = rfmp.m_Fip;
    m_wbt = rfmp.m_wbt;
    m_wbp = rfmp.m_wbp;
    m_wft = rfmp.m_wft;
    
    m_fb.clear();
    for (int ig=0; ig<m_fb.size(); ++ig)
        m_fb.push_back(rfmp.m_fb[ig]);
}

//-----------------------------------------------------------------------------
FEMaterialPointData* FEReactiveFatigueMaterialPoint::Copy()
{
    FEReactiveFatigueMaterialPoint* pt = new FEReactiveFatigueMaterialPoint(*this);
    if (m_pNext) pt->m_pNext = m_pNext->Copy();
    return pt;
}

//-----------------------------------------------------------------------------
void FEReactiveFatigueMaterialPoint::Init()
{
	FEMaterialPointData::Init();
    
    // intialize total damate
    m_D = 0;
    
    // initialize intact bond fraction to 1
    m_wip = m_wit = 1.0;
    
    // initialize intact damage criterion
    m_Ximax = m_Xitrl = m_Xip = 0;
    m_Fip = m_Fit = 0;
    
    // initialize broken and fatigue bond fraction to 0
    m_wbp = m_wbt = m_wft = 0;
    
    // clear the fatigue bond structure
    m_fb.clear();
}

//-----------------------------------------------------------------------------
void FEReactiveFatigueMaterialPoint::Update(const FETimeInfo& timeInfo)
{
	FEMaterialPointData::Update(timeInfo);
    
    // let's check overlapping generations of fatigued bonds
    if (m_fb.size() > 1) {
        for (int ig=0; ig < m_fb.size() - 1; ++ig) {
            double Xfmax = max(m_fb[ig].m_Xftrl,m_fb[ig].m_Xfmax);
            if (Xfmax <= m_fb.back().m_Xftrl) {
                m_fb.back().m_wft += m_fb[ig].m_wft;
                m_fb[ig].m_erase = true;
            }
        }
    }
    // cull generations that have been marked for erasure
    std::deque<FatigueBond>::iterator it = m_fb.begin();
    while (it != m_fb.end()) {
        if (it->m_erase) it = m_fb.erase(it);
        else ++it;
    }
    
    // update intact and damage bonds
    m_wip = m_wit;
    m_wbp = m_wbt;
    
    // update damage response for intact bonds
    if (m_Xitrl > m_Ximax) m_Ximax = m_Xitrl;
    m_Xip = m_Xitrl;
    m_aXip = m_aXit;
    m_Fip = m_Fit;
    
    // update damage response for fatigues bonds
    for (int ig=0; ig<m_fb.size(); ++ig) m_fb[ig].Update();

    // evaluate total fatigue bond fraction
    m_wft = 0;
    for (int ig=0; ig<m_fb.size(); ++ig) m_wft += m_fb[ig].m_wft;
    
    // evaluate total damage at current time
    m_D = 1 - m_wit - m_wft;
    if (m_D < 0) m_D = 0;
    if (m_D > 1) m_D = 1;
}

//-----------------------------------------------------------------------------
void FEReactiveFatigueMaterialPoint::Serialize(DumpStream& ar)
{
	FEMaterialPointData::Serialize(ar);
    ar & m_D;
    ar & m_wit & m_wip & m_Ximax & m_Xitrl & m_Xip & m_aXit & m_aXip & m_Fit & m_Fip;
    ar & m_wbt & m_wbp & m_wft;
    
    // handle deques and boolean
    if (ar.IsSaving()) {
        int n = (int)m_fb.size();
        ar << n;
        for (int i=0; i<n; ++i) {
            ar << m_fb[i].m_wft << m_fb[i].m_wfp;
            ar << m_fb[i].m_Xfmax << m_fb[i].m_Xftrl;
            ar << m_fb[i].m_Fft << m_fb[i].m_Ffp;
            ar << m_fb[i].m_time << m_fb[i].m_erase;
        }
    } else {
        int n;
        ar >> n;
        m_fb.clear();
        for (int i=0; i<n; ++i) {
            FatigueBond fb;
            ar >> fb.m_wft >> fb.m_wfp;
            ar >> fb.m_Xfmax >> fb.m_Xftrl;
            ar >> fb.m_Fft >> fb.m_Ffp;
            ar >> fb.m_time >> fb.m_erase;
            m_fb.push_back(fb);
        }
    }
}
