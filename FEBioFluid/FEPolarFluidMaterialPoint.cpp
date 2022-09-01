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

#include "FEPolarFluidMaterialPoint.h"

//============================================================================
// FEPolarFluidMaterialPoint
//============================================================================
//-----------------------------------------------------------------------------
FEPolarFluidMaterialPoint::FEPolarFluidMaterialPoint(FEMaterialPointData* pt) : FEMaterialPointData(pt)
{
    m_gf = m_gfdot = m_hf = vec3d(0,0,0);
    m_sfa = mat3da(0,0,0);
    m_Psi = mat3d(0,0,0,0,0,0,0,0,0);
}

//-----------------------------------------------------------------------------
FEMaterialPointData* FEPolarFluidMaterialPoint::Copy()
{
    FEPolarFluidMaterialPoint* pt = new FEPolarFluidMaterialPoint(*this);
    if (m_pNext) pt->m_pNext = m_pNext->Copy();
    return pt;
}

//-----------------------------------------------------------------------------
void FEPolarFluidMaterialPoint::Serialize(DumpStream& ar)
{
	FEMaterialPointData::Serialize(ar);
    ar & m_gf & m_gfdot & m_hf & m_sfa & m_Psi;
}

//-----------------------------------------------------------------------------
void FEPolarFluidMaterialPoint::Init()
{
    m_gf = m_gfdot = m_hf = vec3d(0,0,0);
    m_sfa = mat3da(0,0,0);
    m_Psi = mat3d(0,0,0,0,0,0,0,0,0);

    // don't forget to initialize the base class
	FEMaterialPointData::Init();
}
