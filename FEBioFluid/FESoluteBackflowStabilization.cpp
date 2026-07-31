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
#include "FESoluteBackflowStabilization.h"
#include "FEFluidSolutes.h"
#include "FEBioFluidSolutes.h"
#include <FECore/FENodeNodeList.h>
#include <FECore/FEModel.h>

//=============================================================================
BEGIN_FECORE_CLASS(FESoluteBackflowStabilization, FESurfaceLoad)
    ADD_PARAMETER(m_isol   , "solute_id")->setEnums("$(solutes)");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! constructor
FESoluteBackflowStabilization::FESoluteBackflowStabilization(FEModel* pfem) : FESurfaceLoad(pfem), m_dofW(pfem)
{
    m_isol = -1;
    m_dofC = (pfem ? pfem->GetDOFIndex(FEBioFluidSolutes::GetVariableName(FEBioFluidSolutes::FLUID_CONCENTRATION), 0) : -1);
}

//-----------------------------------------------------------------------------
//! allocate storage
void FESoluteBackflowStabilization::SetSurface(FESurface* ps)
{
    FESurfaceLoad::SetSurface(ps);
}

//-----------------------------------------------------------------------------
//! initialize
bool FESoluteBackflowStabilization::Init()
{
    // determine the nr of concentration equations
    FEModel& fem = *GetFEModel();
    DOFS& fedofs = fem.GetDOFS();
    m_dofW.AddVariable(FEBioFluidSolutes::GetVariableName(FEBioFluidSolutes::RELATIVE_FLUID_VELOCITY));
    int MAX_CDOFS = fedofs.GetVariableSize(FEBioFluidSolutes::GetVariableName(FEBioFluidSolutes::FLUID_CONCENTRATION));
    if ((m_isol < 1) || (m_isol > MAX_CDOFS)) return false;
    
    m_dof.AddDofs(m_dofW);
	m_dof.AddVariable(FEBioFluidSolutes::GetVariableName(FEBioFluidSolutes::FLUID_CONCENTRATION));

    FESurface* ps = &GetSurface();
    ps->UpdateNodeNormals();
    
    return FESurfaceLoad::Init();
}

//-----------------------------------------------------------------------------
//! Activate the degrees of freedom for this BC
void FESoluteBackflowStabilization::Activate()
{
    FESurface* ps = &GetSurface();
    
    int dofc = m_dofC + m_isol - 1;
    
    m_ndof.assign(ps->Nodes(),DOF_OPEN);
    
    for (int i=0; i<ps->Nodes(); ++i)
    {
        FENode& node = ps->Node(i);
        m_ndof[i] = node.get_bc(dofc);
    }
    
    FESurfaceLoad::Activate();
}

//-----------------------------------------------------------------------------
//! Evaluate and prescribe the resistance pressure
void FESoluteBackflowStabilization::Update()
{
    // determine backflow conditions
    MarkBackFlow();
    
    // prescribe solute backflow constraint at the nodes
    FESurface* ps = &GetSurface();
    
    int dofc = m_dofC + m_isol - 1;
    
    for (int i=0; i<ps->Nodes(); ++i)
    {
        if (m_ndof[i] == DOF_OPEN) {
            FENode& node = ps->Node(i);
            // set node as having prescribed DOF (concentration at previous time)
            if (node.get_bc(dofc) == DOF_PRESCRIBED) node.set(dofc,node.get_prev(dofc));
        }
    }
}

//-----------------------------------------------------------------------------
//! evaluate the flow rate across this surface
void FESoluteBackflowStabilization::MarkBackFlow()
{
    const FETimeInfo& tp = GetTimeInfo();
    
    // Mark all nodes on this surface to have open concentration DOF
    FESurface* ps = &GetSurface();
    int dofc = m_dofC + m_isol - 1;
    for (int i=0; i<ps->Nodes(); ++i)
    {
        FENode& node = ps->Node(i);
        if (m_ndof[i] == DOF_OPEN) {
            node.set_bc(dofc, DOF_OPEN);
            if (node.m_ID[dofc] < -1) node.m_ID[dofc] = -node.m_ID[dofc] - 2;
            // check velocity normal to surface at this node
            vec3d v = node.get_vec3d(m_dofW[0], m_dofW[1], m_dofW[2]);
            double vn = v*ps->NodeNormal(i);
            if (vn < 0) {
                node.set_bc(dofc, DOF_PRESCRIBED);
                node.m_ID[dofc] = -node.m_ID[dofc] - 2;
            }
        }
    }
}

//-----------------------------------------------------------------------------
//! calculate residual
void FESoluteBackflowStabilization::LoadVector(FEGlobalVector& R)
{
}

//-----------------------------------------------------------------------------
//! serialization
void FESoluteBackflowStabilization::Serialize(DumpStream& ar)
{
    FESurfaceLoad::Serialize(ar);
    if (ar.IsShallow()) return;
	ar & m_dofW & m_dofC;
}
