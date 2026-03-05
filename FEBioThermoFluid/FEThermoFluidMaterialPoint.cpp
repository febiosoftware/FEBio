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

#include "FEThermoFluidMaterialPoint.h"

//============================================================================
FEThermoFluidMaterialPoint::FEThermoFluidMaterialPoint(FEMaterialPointData* pt) : FEMaterialPointData(pt) {}

//-----------------------------------------------------------------------------
FEMaterialPointData* FEThermoFluidMaterialPoint::Copy()
{
    FEThermoFluidMaterialPoint* pt = new FEThermoFluidMaterialPoint(*this);
    if (m_pNext) pt->m_pNext = m_pNext->Copy();
    return pt;
}

//-----------------------------------------------------------------------------
void FEThermoFluidMaterialPoint::Serialize(DumpStream& ar)
{
	FEMaterialPointData::Serialize(ar);
    ar & m_T & m_Tdot & m_gradT & m_k & m_K & m_dKJ & m_dKT & m_cv & m_dcvJ & m_dcvT & m_cp & m_q;
}

//-----------------------------------------------------------------------------
void FEThermoFluidMaterialPoint::Init()
{
    m_T = m_Tdot = 0;
    m_gradT = vec3d(0);
    m_k = 0;
    m_K = m_dKJ = m_dKT = 0;
    m_cv = m_dcvJ = m_dcvT = 0;
    m_cp = 0;
    m_q = vec3d(0);
    
    // don't forget to initialize the base class
	FEMaterialPointData::Init();
}
