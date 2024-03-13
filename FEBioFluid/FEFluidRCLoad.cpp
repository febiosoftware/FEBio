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
#include "FEFluidRCLoad.h"
#include "FEFluid.h"
#include "FEBioFluid.h"
#include <FECore/FEAnalysis.h>
#include <FECore/FEModel.h>

//=============================================================================
BEGIN_FECORE_CLASS(FEFluidRCLoad, FESurfaceLoad)
    ADD_PARAMETER(m_R, "R");
    ADD_PARAMETER(m_p0, "initial_pressure");
    ADD_PARAMETER(m_C, "capacitance");
    ADD_PARAMETER(m_Bern, "Bernoulli");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! constructor
FEFluidRCLoad::FEFluidRCLoad(FEModel* pfem) : FESurfaceLoad(pfem), m_dofW(pfem)
{
    m_R = 0.0;
    m_pfluid = nullptr;
    m_p0 = 0;
    m_C = 0.0;
    m_Bern = false;

    m_pt = m_dpt = 0;
    m_qt = m_dqt = 0;
    m_pp = m_dpp = 0;
    m_qp = m_dqp = 0;

    // TODO: Can this be done in Init, since  there is no error checking
    if (pfem)
    {
        m_dofW.AddVariable(FEBioFluid::GetVariableName(FEBioFluid::RELATIVE_FLUID_VELOCITY));
        m_dofEF = pfem->GetDOFIndex(FEBioFluid::GetVariableName(FEBioFluid::FLUID_DILATATION), 0);
    }
}

//-----------------------------------------------------------------------------
//! initialize
//! TODO: Generalize to include the initial conditions
bool FEFluidRCLoad::Init()
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

    return true;
}

//-----------------------------------------------------------------------------
//! Activate the degrees of freedom for this BC
void FEFluidRCLoad::Activate()
{
    FESurface* ps = &GetSurface();

    for (int i = 0; i < ps->Nodes(); ++i)
    {
        FENode& node = ps->Node(i);
        // mark node as having prescribed DOF
        node.set_bc(m_dofEF, DOF_PRESCRIBED);
    }
}

//-----------------------------------------------------------------------------
//! Evaluate and prescribe the resistance pressure
void FEFluidRCLoad::Update()
{
    const FETimeInfo& tp = GetTimeInfo();

    double dt = tp.timeIncrement;
    double gamma = tp.gamma;
    double omoog = 1 - 1. / gamma;
    int niter = GetFEModel()->GetCurrentStep()->GetFESolver()->m_niter;

    // if this is the start of a new time step, update the flow/pressure parameters
    if (niter == 0) {
        m_qp = m_qt;
        m_dqp = m_dqt * omoog;
        m_pp = m_pt;
        m_dpp = m_dpt * omoog;
    }

    // evaluate the flow rate
    m_qt = FlowRate();

    // evaluate the outflow pressure
    double qt = m_Bern ? m_qt * fabs(m_qt) : m_qt;
    double qp = m_Bern ? m_qp * fabs(m_qp) : m_qp;
    m_dqt = m_dqp * omoog + (qt - qp) / gamma / dt;
    m_dpt = m_dpp * omoog + (m_pt - m_pp) / gamma / dt;
    m_pt = m_pp + (m_R * m_dqt + m_qt / m_C - m_dpp * omoog) * gamma * dt;

    // calculate the dilatation
    double e = 0;
    bool good = m_pfluid->Dilatation(0, m_pt + m_p0, e);
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

    GetFEModel()->SetMeshUpdateFlag(true);
}

//-----------------------------------------------------------------------------
//! evaluate the flow rate across this surface
double FEFluidRCLoad::FlowRate()
{
    double Q = 0;

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
void FEFluidRCLoad::LoadVector(FEGlobalVector& R)
{
}

//-----------------------------------------------------------------------------
//! serialization
void FEFluidRCLoad::Serialize(DumpStream& ar)
{
	FESurfaceLoad::Serialize(ar);
    ar & m_pt & m_dpt & m_qt & m_dqt;
    if (ar.IsShallow()) return;
    ar & m_pfluid;
    ar & m_dofW & m_dofEF;
}
