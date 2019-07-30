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
#include "FEFluidSolutes.h"
#include <FECore/FECoreKernel.h>
#include <FECore/DumpStream.h>

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEFluidSolutes, FEMaterial)
// material properties
ADD_PROPERTY(m_pFluid, "fluid");
ADD_PROPERTY(m_pSolute, "solute"             , FEProperty::Optional);
END_FECORE_CLASS();

//============================================================================
// FEFluidSolutesMaterialPoint
//============================================================================
FEFluidSolutesMaterialPoint::FEFluidSolutesMaterialPoint(FEMaterialPoint* pt) : FEMaterialPoint(pt) {}

//-----------------------------------------------------------------------------
FEMaterialPoint* FEFluidSolutesMaterialPoint::Copy()
{
    FEFluidSolutesMaterialPoint* pt = new FEFluidSolutesMaterialPoint(*this);
    if (m_pNext) pt->m_pNext = m_pNext->Copy();
    return pt;
}

//-----------------------------------------------------------------------------
void FEFluidSolutesMaterialPoint::Serialize(DumpStream& ar)
{
    FEMaterialPoint::Serialize(ar);
    ar & m_nsol;
    ar & m_c & m_gradc & m_j & m_cdot;
}

//-----------------------------------------------------------------------------
void FEFluidSolutesMaterialPoint::Init()
{
    m_nsol = 0;
    m_c.clear();
    m_gradc.clear();
    m_j.clear();
    m_cdot.clear();

    FEMaterialPoint::Init();
}

//============================================================================
// FEFluidSolutes
//============================================================================

//-----------------------------------------------------------------------------
//! FEFluidSolutes constructor

FEFluidSolutes::FEFluidSolutes(FEModel* pfem) : FEMaterial(pfem)
{
    m_pFluid = 0;
}

//-----------------------------------------------------------------------------
// returns a pointer to a new material point object
FEMaterialPoint* FEFluidSolutes::CreateMaterialPointData()
{
    FEFluidMaterialPoint* fpt = new FEFluidMaterialPoint();
    return new FEFluidSolutesMaterialPoint(fpt);
}

//-----------------------------------------------------------------------------
// initialize
bool FEFluidSolutes::Init()
{
    // set the solute IDs first, since they are referenced in FESolute::Init()
    for (int i = 0; i<Solutes(); ++i) {
        m_pSolute[i]->SetSoluteLocalID(i);
    }
    
    // call the base class.
    // This also initializes all properties
    return FEMaterial::Init();
}

//-----------------------------------------------------------------------------
//! actual concentration
double FEFluidSolutes::Concentration(FEMaterialPoint& pt, const int sol)
{
    FEFluidSolutesMaterialPoint& spt = *pt.ExtractData<FEFluidSolutesMaterialPoint>();
    
    // effective concentration
    double ca = spt.m_c[sol];
    
    return ca;
}

//-----------------------------------------------------------------------------
//! Calculate solute molar flux

vec3d FEFluidSolutes::SoluteFlux(FEMaterialPoint& pt, const int sol)
{
    FEFluidSolutesMaterialPoint& spt = *pt.ExtractData<FEFluidSolutesMaterialPoint>();
    
    // concentration gradient
    vec3d gradc = spt.m_gradc[sol];
    
    // solute free diffusivity
    double D0 = m_pSolute[sol]->m_pDiff->Free_Diffusivity(pt);
    
    // solute flux j
    vec3d j = -gradc*D0;
    
    return j;
}
