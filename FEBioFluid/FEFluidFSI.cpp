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
#include "FEFluidFSI.h"
#include <FECore/FECoreKernel.h>
#include <FECore/DumpStream.h>

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEFluidFSI, FEMaterial)
	// material properties
	ADD_PROPERTY(m_pSolid, "solid", FEProperty::Required | FEProperty::TopLevel);
	ADD_PROPERTY(m_pFluid, "fluid");
END_FECORE_CLASS();

//============================================================================
// FEFSIMaterialPoint
//============================================================================
FEFSIMaterialPoint::FEFSIMaterialPoint(FEMaterialPointData* pt) : FEMaterialPointData(pt) {}

//-----------------------------------------------------------------------------
FEMaterialPointData* FEFSIMaterialPoint::Copy()
{
    FEFSIMaterialPoint* pt = new FEFSIMaterialPoint(*this);
    if (m_pNext) pt->m_pNext = m_pNext->Copy();
    return pt;
}

//-----------------------------------------------------------------------------
void FEFSIMaterialPoint::Serialize(DumpStream& ar)
{
	FEMaterialPointData::Serialize(ar);
	ar & m_w & m_aw & m_Jdot & m_ss;
}

//-----------------------------------------------------------------------------
void FEFSIMaterialPoint::Init()
{
    m_w = m_aw = vec3d(0,0,0);
    m_Jdot = 0;
    m_ss.zero();
    
	FEMaterialPointData::Init();
}

//============================================================================
// FEFluidFSI
//============================================================================

//-----------------------------------------------------------------------------
//! FEFluidFSI constructor

FEFluidFSI::FEFluidFSI(FEModel* pfem) : FEMaterial(pfem)
{
	m_pSolid = 0;
	m_pFluid = 0;
}

//-----------------------------------------------------------------------------
// returns a pointer to a new material point object
FEMaterialPointData* FEFluidFSI::CreateMaterialPointData()
{
    FEFluidMaterialPoint* fpt = new FEFluidMaterialPoint(m_pSolid->CreateMaterialPointData());
    return new FEFSIMaterialPoint(fpt);
}

//-----------------------------------------------------------------------------
// initialize
bool FEFluidFSI::Init()
{
    m_pFluid->Init();
    m_pSolid->Init();
    // set the solid density to zero (required for the solid of a FSI domain)
    m_pSolid->SetDensity(0.0);
    
    return FEMaterial::Init();
}
