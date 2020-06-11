/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2020 University of Utah, The Trustees of Columbia University in 
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
#include "FESolutesMaterial.h"
#include <FECore/FECoreKernel.h>
#include <FECore/DumpStream.h>

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FESolutesMaterial, FEMaterial)
	ADD_PROPERTY(m_pSolute, "solute", FEProperty::Optional);
END_FECORE_CLASS();

//============================================================================
// FEFluidSolutesMaterialPoint
//============================================================================
FESolutesMaterial::Point::Point(FEMaterialPoint* pt) : FEMaterialPoint(pt) 
{
	m_vft = vec3d(0.0, 0.0, 0.0);
	m_divf = 0.0;
}

//-----------------------------------------------------------------------------
FEMaterialPoint* FESolutesMaterial::Point::Copy()
{
	FESolutesMaterial::Point* pt = new FESolutesMaterial::Point(*this);
    if (m_pNext) pt->m_pNext = m_pNext->Copy();
    return pt;
}

//-----------------------------------------------------------------------------
void FESolutesMaterial::Point::Serialize(DumpStream& ar)
{
    FEMaterialPoint::Serialize(ar);
	ar & m_vft & m_divf;
    ar & m_nsol;
    ar & m_c & m_gradc & m_j & m_cdot;
}

//-----------------------------------------------------------------------------
void FESolutesMaterial::Point::Init()
{
    m_nsol = 0;
    m_c.clear();
    m_gradc.clear();
    m_j.clear();
    m_cdot.clear();

    FEMaterialPoint::Init();
}

//============================================================================
// FESolutesMaterial
//============================================================================

//-----------------------------------------------------------------------------
//! FEFluidSolutes constructor

FESolutesMaterial::FESolutesMaterial(FEModel* pfem) : FEMaterial(pfem)
{
}

//-----------------------------------------------------------------------------
// returns a pointer to a new material point object
FEMaterialPoint* FESolutesMaterial::CreateMaterialPointData()
{
    return new FESolutesMaterial::Point(nullptr);
}

//-----------------------------------------------------------------------------
// initialize
bool FESolutesMaterial::Init()
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
double FESolutesMaterial::Concentration(FEMaterialPoint& pt, const int sol)
{
	FESolutesMaterial::Point& spt = *pt.ExtractData<FESolutesMaterial::Point>();
    
    // effective concentration
    double ca = spt.m_c[sol];
    
    return ca;
}

//-----------------------------------------------------------------------------
//! Calculate solute molar flux

vec3d FESolutesMaterial::SoluteFlux(FEMaterialPoint& pt, const int sol)
{
	FESolutesMaterial::Point& spt = *pt.ExtractData<FESolutesMaterial::Point>();
    
    // concentration gradient
    vec3d gradc = spt.m_gradc[sol];
    
    // solute free diffusivity
    double D0 = m_pSolute[sol]->m_pDiff->Free_Diffusivity(pt);
    
    // solute flux j
    vec3d j = -gradc*D0;
    
    return j;
}
