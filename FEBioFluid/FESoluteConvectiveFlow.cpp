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
#include "FESoluteConvectiveFlow.h"
#include "FECore/FEElemElemList.h"
#include "FECore/FEGlobalMatrix.h"
#include "FECore/FEGlobalVector.h"
#include "FECore/LinearSolver.h"
#include <FECore/FEModel.h>
#include "FEBioFluidSolutes.h"

//=============================================================================
BEGIN_FECORE_CLASS(FESoluteConvectiveFlow, FESurfaceLoad)
    ADD_PARAMETER(m_sol   , "solute_id")->setEnums("$(solutes)");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! constructor
FESoluteConvectiveFlow::FESoluteConvectiveFlow(FEModel* pfem) : FESurfaceLoad(pfem), m_dofW(pfem)
{
    m_sol = -1;
    
    m_dofW.AddVariable(FEBioFluidSolutes::GetVariableName(FEBioFluidSolutes::RELATIVE_FLUID_VELOCITY));
    m_dofEF = GetDOFIndex(FEBioFluidSolutes::GetVariableName(FEBioFluidSolutes::FLUID_DILATATION), 0);
    m_dofC = GetDOFIndex(FEBioFluidSolutes::GetVariableName(FEBioFluidSolutes::FLUID_CONCENTRATION), 0);
}

//-----------------------------------------------------------------------------
//! initialize
bool FESoluteConvectiveFlow::Init()
{
    // determine the nr of concentration equations
    FEModel& fem = *GetFEModel();
    DOFS& fedofs = fem.GetDOFS();
    int MAX_CDOFS = fedofs.GetVariableSize(FEBioFluidSolutes::GetVariableName(FEBioFluidSolutes::FLUID_CONCENTRATION));
    if ((m_sol < 1) || (m_sol > MAX_CDOFS)) return false;
    
    FEMesh& mesh = GetMesh();
    m_octree = new FEOctreeSearch(&mesh);
    m_octree->Init();
    
    FESurface* ps = &GetSurface();
    m_np = new FENormalProjection(*ps);
    m_np->SetTolerance(0.01);
    m_np->SetSearchRadius(1.0);
    m_np->Init();
    
    int NN = mesh.Nodes();
    m_bexclude.assign(NN, false);
    
    return FESurfaceLoad::Init();
}

//-----------------------------------------------------------------------------
//! Activate the degrees of freedom for this BC
void FESoluteConvectiveFlow::Activate()
{
    int dofc = m_dofC + m_sol - 1;
    
    FEModel& fem = *GetFEModel();
    FEMesh& mesh = GetMesh();
    
    for (int i=0; i<mesh.Nodes(); ++i)
    {
        FENode& node = mesh.Node(i);
        if (node.get_bc(dofc) == DOF_PRESCRIBED) {
            m_bexclude[i] = true;
        }
        else if (node.get_bc(dofc) == DOF_OPEN) {
            // mark node as having prescribed DOF
            node.set_bc(dofc, DOF_PRESCRIBED);
        }
    }
    
    FESurfaceLoad::Activate();
}

//-----------------------------------------------------------------------------
// return nodal value
void FESoluteConvectiveFlow::Update()
{
    FEMesh& mesh = GetMesh();

    const FETimeInfo& tp = GetTimeInfo();
    double gamma = tp.gamma;
    double dt = tp.timeIncrement;
    
    for (int i=0; i<mesh.Nodes(); ++i)
    {
        if (!m_bexclude[i]) {
            FENode& node = mesh.Node(i);
            vec3d x = node.m_rt;
            vec3d vt = node.get_vec3d(m_dofW[0], m_dofW[1], m_dofW[2]);
            vec3d vp = node.get_vec3d_prev(m_dofW[0], m_dofW[1], m_dofW[2]);
            
            vec3d X = x - (vt*gamma + vp*(1-gamma))*dt;
            
            int dofc = m_dofC + m_sol - 1;
            double r[3] = { 0 };
            double c = 0;
            
            // search for the solid element in which X lies
            FESolidElement* el = (FESolidElement*)m_octree->FindElement(X, r);
            if (el) {
                const int NELN = FESolidElement::MAX_NODES;
                double ep[NELN], cp[NELN];
                int neln = el->Nodes();
                for (int j=0; j<neln; ++j) {
                    FENode& node = mesh.Node(el->m_node[j]);
                    ep[j] = node.get_prev(m_dofEF);
                    cp[j] = node.get_prev(dofc);
                }
                double Jt = 1 + node.get(m_dofEF);
                double Jp = 1 + el->evaluate(ep, r[0], r[1], r[2]);
                c = Jp*el->evaluate(cp, r[0], r[1], r[2])/Jt;
            }
            // if solid element is not found, project x onto the solute inlet surface
            else {
                FESurfaceElement* pme;
                vec3d n = x - X;
                n.unit();
                pme = m_np->Project(x, n, r);
                if (pme) {
                    const int NELN = FEShellElement::MAX_NODES;
                    double ep[NELN], cp[NELN];
                    int neln = pme->Nodes();
                    for (int j=0; j<neln; ++j) {
                        FENode& mode = mesh.Node(pme->m_node[j]);
                        ep[j] = mode.get_prev(m_dofEF);
                        cp[j] = mode.get_prev(dofc);
                    }
                    double Jt = 1 + node.get(m_dofEF);
                    double Jp = 1 + pme->eval(ep, r[0], r[1]);
                    c = Jp*pme->eval(cp, r[0], r[1])/Jt;
                }
                else
                    c = node.get_prev(dofc);
            }
            
            if (node.m_ID[dofc] < -1)
                node.set(dofc, c);
        }
    }
}

//-----------------------------------------------------------------------------
//! serialization
void FESoluteConvectiveFlow::Serialize(DumpStream& ar)
{
    FESurfaceLoad::Serialize(ar);
    if (ar.IsShallow()) return;
    ar & m_dofW & m_dofEF & m_dofC;
}
