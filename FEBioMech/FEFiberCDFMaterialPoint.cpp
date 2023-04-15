/*This file is part of the FEBio source code and is licensed under the MIT license
 listed below.
 
 See Copyright-FEBio.txt for details.
 
 Copyright (c) 2023 University of Utah, The Trustees of Columbia University in
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



#include "FEFiberCDFMaterialPoint.h"
#include <FECore/DumpStream.h>

//-----------------------------------------------------------------------------
// FEFiberCDFMaterialPoint
//-----------------------------------------------------------------------------

FEFiberCDFMaterialPoint::FEFiberCDFMaterialPoint(FEMaterialPointData* pt) : FEMaterialPointData(pt)
{
    m_In_1_p = m_In_1_t = 0;
    m_sed_t = m_sed_p = 0;
    m_ds_t = m_ds_p = 0;
    m_d2s_p = m_d2s_t = 0;
}

//-----------------------------------------------------------------------------
FEMaterialPointData* FEFiberCDFMaterialPoint::Copy()
{
    FEFiberCDFMaterialPoint* pt = new FEFiberCDFMaterialPoint(*this);
    if (m_pNext) pt->m_pNext = m_pNext->Copy();
    return pt;
}

//-----------------------------------------------------------------------------
void FEFiberCDFMaterialPoint::Init()
{
    // initialize data
    m_In_1_p = m_In_1_t = 0;
    m_sed_t = m_sed_p = 0;
    m_ds_t = m_ds_p = 0;
    m_d2s_p = m_d2s_t = 0;
    
    // don't forget to intialize the nested data
    FEMaterialPointData::Init();
}

//-----------------------------------------------------------------------------
void FEFiberCDFMaterialPoint::Update(const FETimeInfo& timeInfo)
{
    m_In_1_p = m_In_1_t;
    m_sed_p = m_sed_t;
    m_ds_p = m_ds_t;
    m_d2s_p = m_d2s_t;
}

//-----------------------------------------------------------------------------
void FEFiberCDFMaterialPoint::Serialize(DumpStream& ar)
{
    FEMaterialPointData::Serialize(ar);
    ar & m_In_1_t & m_In_1_p;
    ar & m_sed_t & m_sed_p;
    ar & m_ds_t & m_ds_p;
    ar & m_d2s_t & m_d2s_p;
}

//-----------------------------------------------------------------------------
// Perform integration in one step
void FEFiberCDFMaterialPoint::Integrate(const double cdf)
{
    m_d2s_t = (cdf > 0) ? cdf : 0;

    double dIn = m_In_1_t - m_In_1_p;
    m_sed_t = m_sed_p + m_ds_p*dIn + 0.25*dIn*dIn*(m_d2s_t + m_d2s_p);
    m_ds_t = m_ds_p + 0.5*dIn*(m_d2s_t + m_d2s_p);

}

//-----------------------------------------------------------------------------
void FEFiberCDFMaterialPoint::SetFiberStrain(const double In_1) {
    m_In_1_t = In_1;
    if (m_In_1_t <= 0) {
        m_d2s_t = 0;
    }
}
