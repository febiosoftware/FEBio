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
#include "FEThermoFluidPressureLoad.h"
#include "FEFluid.h"
#include "FEBioThermoFluid.h"

//=============================================================================
BEGIN_FECORE_CLASS(FEThermoFluidPressureLoad, FESurfaceLoad)
    ADD_PARAMETER(m_p0, "pressure");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! constructor
FEThermoFluidPressureLoad::FEThermoFluidPressureLoad(FEModel* pfem) : FESurfaceLoad(pfem)
{
    m_alpha = 1.0;
    m_p0 = 0;
    
    m_dofEF = pfem->GetDOFIndex(FEBioThermoFluid::GetVariableName(FEBioThermoFluid::FLUID_DILATATION), 0);
    m_dofT = pfem->GetDOFIndex(FEBioThermoFluid::GetVariableName(FEBioThermoFluid::TEMPERATURE), 0);
}

//-----------------------------------------------------------------------------
//! initialize
bool FEThermoFluidPressureLoad::Init()
{
    if (FESurfaceLoad::Init() == false) return false;

    // get fluid from first surface element
    // assuming the entire surface bounds the same fluid
    FESurfaceElement& el = m_psurf->Element(0);
    FEElement* pe = el.m_elem[0];
    if (pe == nullptr) return false;

    // get the material
    FEMaterial* pm = GetFEModel()->GetMaterial(pe->GetMatID());
    m_pfluid = pm->ExtractProperty<FEThermoFluid>();
    if (m_pfluid == nullptr) return false;
    
    return true;
}

//-----------------------------------------------------------------------------
//! Activate the degrees of freedom for this BC
void FEThermoFluidPressureLoad::Activate()
{
    FESurface* ps = &GetSurface();
    
    for (int i=0; i<ps->Nodes(); ++i)
    {
        FENode& node = ps->Node(i);
        // mark node as having prescribed DOF
        node.set_bc(m_dofEF, DOF_PRESCRIBED);
    }
}

//-----------------------------------------------------------------------------
//! Evaluate and prescribe the resistance pressure
void FEThermoFluidPressureLoad::Update()
{
    // prescribe this dilatation at the nodes
    FESurface* ps = &GetSurface();

    for (int i=0; i<ps->Nodes(); ++i)
    {
        if (ps->Node(i).m_ID[m_dofEF] < -1)
        {
            FENode& node = ps->Node(i);
            double T = node.get(m_dofT);
            // calculate the dilatation
            double e = m_pfluid->Dilatation(T,m_p0);
            
            // set node as having prescribed DOF
            node.set(m_dofEF, e);
        }
    }
}

//-----------------------------------------------------------------------------
//! calculate residual
void FEThermoFluidPressureLoad::LoadVector(FEGlobalVector& R, const FETimeInfo& tp)
{
    m_alpha = tp.alpha; m_alphaf = tp.alphaf;
}

//-----------------------------------------------------------------------------
//! serialization
void FEThermoFluidPressureLoad::Serialize(DumpStream& ar)
{
    FESurfaceLoad::Serialize(ar);
    ar & m_alpha & m_alphaf;
    ar & m_pfluid;
}
