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
#include "FEFluidRCRLoad.h"
#include "FEFluid.h"
#include "FEBioFluid.h"
#include <FECore/FEAnalysis.h>
#include <FECore/log.h>
#include <FECore/FEModel.h>

//=============================================================================
BEGIN_FECORE_CLASS(FEFluidRCRLoad, FESurfaceLoad)
    ADD_PARAMETER(m_R , "R")->setUnits("F.t/L^5");
    ADD_PARAMETER(m_Rd, "Rd")->setUnits("F.t/L^5");
    ADD_PARAMETER(m_p0, "initial_pressure")->setUnits(UNIT_PRESSURE);
    ADD_PARAMETER(m_pd, "pressure_offset")->setUnits(UNIT_PRESSURE);
    ADD_PARAMETER(m_C , "capacitance")->setUnits("L^5/F");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! constructor
FEFluidRCRLoad::FEFluidRCRLoad(FEModel* pfem) : FESurfaceLoad(pfem), m_dofW(pfem)
{
    m_R = 0.0;
    m_pfluid = nullptr;
    m_p0 = 0;
    m_Rd = 0.0;
    m_pd = 0.0;
    m_C = 0.0;

    m_dofW.AddVariable(FEBioFluid::GetVariableName(FEBioFluid::RELATIVE_FLUID_VELOCITY));
    m_dofEF = pfem->GetDOFIndex(FEBioFluid::GetVariableName(FEBioFluid::FLUID_DILATATION), 0);
}

//-----------------------------------------------------------------------------
//! initialize
//! TODO: Generalize to include the initial conditions
bool FEFluidRCRLoad::Init()
{
    if (FESurfaceLoad::Init() == false) return false;

    m_dof.Clear();
    m_dof.AddDofs(m_dofW);
    m_dof.AddDof(m_dofEF);

    // get fluid from first surface element
    // assuming the entire surface bounds the same fluid
    FESurfaceElement& el = m_psurf->Element(0);
    FEElement* pe = el.m_elem[0];
    if (pe == nullptr) return false;

    // get the material
    FEMaterial* pm = GetFEModel()->GetMaterial(pe->GetMatID());
    m_pfluid = pm->ExtractProperty<FEFluidMaterial>();
    if (m_pfluid == nullptr) return false;

    m_pn = m_pp = m_p0;
    m_pdn = m_pdp = m_pd;
    m_qn = m_qp = 0;
    m_tp = 0;

    return true;
}

//-----------------------------------------------------------------------------
//! Activate the degrees of freedom for this BC
void FEFluidRCRLoad::Activate()
{
    FESurface* ps = &GetSurface();

    for (int i = 0; i < ps->Nodes(); ++i)
    {
        FENode& node = ps->Node(i);
        // mark node as having prescribed DOF
        node.set_bc(m_dofEF, DOF_PRESCRIBED);
    }
    
    FESurfaceLoad::Activate();
}

//-----------------------------------------------------------------------------
//! Evaluate and prescribe the resistance pressure
void FEFluidRCRLoad::Update()
{
    // Check if we started a new time, if so, update variables
    FETimeInfo& timeInfo = GetFEModel()->GetTime();
    double time = timeInfo.currentTime;
    int iter = timeInfo.currentIteration;
    double dt = timeInfo.timeIncrement;
    if ((time > m_tp) && (iter == 0)) {
        m_pp = m_pn;
        m_qp = m_qn;
        m_pdp = m_pdn;
        m_tp = time;
    }

    // evaluate the flow rate at the current time
    m_qn = FlowRate();
    m_pdn = m_pd;

    double tau = m_Rd * m_C;

    // calculate the RCR pressure
    m_pn = m_pdn + (m_Rd / (1 + tau / dt) + m_R) * m_qn + tau / (dt + tau) * (m_pp - m_pdp - m_R * m_qp);

    // calculate the dilatation
    double e = 0;
    bool good = m_pfluid->Dilatation(0, m_pn, e);
    assert(good);

    // prescribe this dilatation at the nodes
    FESurface* ps = &GetSurface();

    for (int i = 0; i < ps->Nodes(); ++i)
    {
        if (ps->Node(i).m_ID[m_dofEF] < -1)
        {
            FENode& node = ps->Node(i);
            // set node as having prescribed DOF
            node.set(m_dofEF, e);
        }
    }

    // Force a mesh update after loads have been updated
    ForceMeshUpdate();
}

//-----------------------------------------------------------------------------
//! evaluate the flow rate across this surface at current time
double FEFluidRCRLoad::FlowRate()
{
    double Q = 0;

    const FETimeInfo& tp = GetTimeInfo();

    vec3d rt[FEElement::MAX_NODES];
    vec3d vt[FEElement::MAX_NODES];

    for (int iel = 0; iel < m_psurf->Elements(); ++iel)
    {
        FESurfaceElement& el = m_psurf->Element(iel);

        // nr integration points
        int nint = el.GaussPoints();

        // nr of element nodes
        int neln = el.Nodes();

        // nodal coordinates
        for (int i = 0; i < neln; ++i) {
            FENode& node = m_psurf->GetMesh()->Node(el.m_node[i]);
            rt[i] = node.m_rt;
            vt[i] = node.get_vec3d(m_dofW[0], m_dofW[1], m_dofW[2]);
        }

        double* Nr, * Ns;
        double* N;
        double* w = el.GaussWeights();

        vec3d dxr, dxs, v;

        // repeat over integration points
        for (int n = 0; n < nint; ++n)
        {
            N = el.H(n);
            Nr = el.Gr(n);
            Ns = el.Gs(n);

            // calculate the velocity and tangent vectors at integration point
            dxr = dxs = v = vec3d(0, 0, 0);
            for (int i = 0; i < neln; ++i)
            {
                v += vt[i] * N[i];
                dxr += rt[i] * Nr[i];
                dxs += rt[i] * Ns[i];
            }

            vec3d normal = dxr ^ dxs;
            double q = normal * v;
            Q += q * w[n];
        }
    }

    return Q;
}

//-----------------------------------------------------------------------------
//! calculate residual
void FEFluidRCRLoad::LoadVector(FEGlobalVector& R)
{
}

//-----------------------------------------------------------------------------
//! serialization
void FEFluidRCRLoad::Serialize(DumpStream& ar)
{
    FESurfaceLoad::Serialize(ar);
    ar & m_pn& m_pp& m_qn& m_qp& m_pdn& m_pdp& m_tp;
    if (ar.IsShallow()) return;
    ar & m_pfluid;
    ar & m_dofW & m_dofEF;
}
