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
#include "FEFluidPressureLoad.h"
#include "FECore/FECoreKernel.h"
#include "FEBioFluid.h"
#include <FECore/FEModel.h>

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEFluidPressureLoad, FESurfaceLoad)
ADD_PARAMETER(m_p0    , "pressure");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEFluidPressureLoad::FEFluidPressureLoad(FEModel* pfem) : FESurfaceLoad(pfem)
{
    m_pfluid = nullptr;
    m_p0 = 0;
    
    m_dofEF = pfem->GetDOFIndex(FEBioFluid::GetVariableName(FEBioFluid::FLUID_DILATATION), 0);
    
    m_dof.Clear();
    m_dof.AddDof(m_dofEF);
}

//-----------------------------------------------------------------------------
bool FEFluidPressureLoad::Init()
{
    if (FESurfaceLoad::Init() == false) return false;
    
    // get fluid from first surface element
    // assuming the entire surface bounds the same fluid
    FESurfaceElement& el = m_psurf->Element(0);
    FEElement* pe = el.m_elem[0];
    if (pe == nullptr) return false;
    
    // get the material
    FEMaterial* pm = GetFEModel()->GetMaterial(pe->GetMatID());
    m_pfluid = pm->ExtractProperty<FEFluidMaterial>();
    if (m_pfluid == nullptr) return false;
    
    return true;
}

//-----------------------------------------------------------------------------
//! Activate the degrees of freedom for this BC
void FEFluidPressureLoad::Activate()
{
    FESurface* ps = &GetSurface();
    
    for (int i=0; i<ps->Nodes(); ++i)
    {
        FENode& node = ps->Node(i);
        // mark node as having prescribed DOF
        node.set_bc(m_dofEF, DOF_PRESCRIBED);
    }
    
    FESurfaceLoad::Activate();
}

//-----------------------------------------------------------------------------
//! Evaluate and prescribe the resistance pressure
void FEFluidPressureLoad::Update()
{
    // prescribe this dilatation at the nodes
    FESurface* ps = &GetSurface();
    
    for (int i=0; i<ps->Nodes(); ++i)
    {
        if (ps->Node(i).m_ID[m_dofEF] < -1)
        {
            FENode& node = ps->Node(i);
            // calculate the dilatation
            double e = node.get(m_dofEF);
            bool good = m_pfluid->Dilatation(0,m_p0,e);
            assert(good);
            // set node as having prescribed DOF
            node.set(m_dofEF, e);
        }
    }
    
    GetFEModel()->SetMeshUpdateFlag(true);
}

//-----------------------------------------------------------------------------
//! serialization
void FEFluidPressureLoad::Serialize(DumpStream& ar)
{
    FESurfaceLoad::Serialize(ar);
    if (ar.IsShallow()) return;
    ar & m_dofEF;
    ar & m_pfluid;
}
