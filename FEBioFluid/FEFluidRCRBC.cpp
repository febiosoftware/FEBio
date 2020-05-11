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
#include "FEFluidRCRBC.h"
#include "FEFluid.h"
#include "FEBioFluid.h"
#include <FECore/FEAnalysis.h>

//=============================================================================
BEGIN_FECORE_CLASS(FEFluidRCRBC, FESurfaceLoad)
ADD_PARAMETER(m_R , "R");
ADD_PARAMETER(m_Rd , "Rd");
ADD_PARAMETER(m_p0, "initial_pressure");
ADD_PARAMETER(m_pd, "pressure_offset");
ADD_PARAMETER(m_C, "capacitance");
ADD_PARAMETER(m_Bern, "Bernoulli");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! constructor
FEFluidRCRBC::FEFluidRCRBC(FEModel* pfem) : FESurfaceLoad(pfem), m_dofW(pfem)
{
    m_R = 0.0;
    m_pfluid = nullptr;
    m_alpha = 1.0;
    m_p0 = 0;
    m_Rd = 0.0;
    m_pd = 0.0;
    m_C = 0.0;
    m_Bern = 0;
    
    m_stepHist.clear();
    m_timeHist.clear();
    m_flowHist.clear();
    
    m_dofW.AddVariable(FEBioFluid::GetVariableName(FEBioFluid::RELATIVE_FLUID_VELOCITY));
    m_dofEF = pfem->GetDOFIndex(FEBioFluid::GetVariableName(FEBioFluid::FLUID_DILATATION), 0);
}

//-----------------------------------------------------------------------------
//! initialize
//! TODO: Generalize to include the initial conditions
bool FEFluidRCRBC::Init()
{
    if (FESurfaceLoad::Init() == false) return false;
    
    // get fluid from first surface element
    // assuming the entire surface bounds the same fluid
    FESurfaceElement& el = m_psurf->Element(0);
    FEMesh* mesh = m_psurf->GetMesh();
    FEElement* pe = el.m_elem[0];
    if (pe == nullptr) return false;
    
    // get the material
    FEMaterial* pm = GetFEModel()->GetMaterial(pe->GetMatID());
    m_pfluid = pm->ExtractProperty<FEFluidMaterial>();
    if (m_pfluid == nullptr) return false;
    
    m_stepHist.resize(1,0.0);
    m_timeHist.resize(1,0.0);
    m_flowHist.resize(1,0.0);
    
    return true;
}

//-----------------------------------------------------------------------------
//! Activate the degrees of freedom for this BC
void FEFluidRCRBC::Activate()
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
void FEFluidRCRBC::Update()
{
    // evaluate the flow rate
    double Q = FlowRate();
    
    int numsteps = GetFEModel()->GetCurrentStep()->m_ntimesteps;
    FETimeInfo& tp = GetFEModel()->GetTime();
    
    m_flowHist.resize(numsteps + 1, 0);
    m_timeHist.resize(numsteps + 1);
    m_stepHist.resize(numsteps + 1);
    
    m_timeHist[numsteps] = tp.currentTime;
    m_stepHist[numsteps] = tp.timeIncrement;
    
    m_flowHist[numsteps] = Q;
    
    double tau = m_Rd*m_C;
    
    // calculate the resistance pressure
    double pR = 0;
    if (m_Bern == 1 || m_Bern == 2)
        pR = m_R*Q*abs(Q);
    else
        pR = m_R*Q;
    
    //calculate initial pressure contribution
    double pi = 0.0;
    if (tau > 0)
        pi = m_p0*exp(-m_timeHist[numsteps]/tau);
    
    //calculate pressure from capacitor contribution
    double pC = 0.0;
    if (m_C > 0 && m_Rd > 0)
    {
        for (int i = 0; i<=numsteps; ++i)
        {
            double p1 = 0;
            double p2 = 0;
            if (m_Bern == 1)
            {
                if (i != 0)
                    p1 = exp(-(m_timeHist[numsteps]-m_timeHist[i-1])/tau)/m_C*m_flowHist[i-1]*abs(m_flowHist[i-1]);
                p2 = exp(-(m_timeHist[numsteps]-m_timeHist[i])/tau)/m_C*m_flowHist[i]*abs(m_flowHist[i]);
            }
            else
            {
                if (i != 0)
                    p1 = exp(-(m_timeHist[numsteps]-m_timeHist[i-1])/tau)/m_C*m_flowHist[i-1];
                p2 = exp(-(m_timeHist[numsteps]-m_timeHist[i])/tau)/m_C*m_flowHist[i];
            }
            pC += (p1+p2)/2.0*m_stepHist[i];
        }
    }
    
    double p = pR + pi + m_pd + pC;
    
    // calculate the dilatation
    double e = m_pfluid->Dilatation(0,p);
    
    // prescribe this dilatation at the nodes
    FESurface* ps = &GetSurface();
    
    for (int i=0; i<ps->Nodes(); ++i)
    {
        if (ps->Node(i).m_ID[m_dofEF] < -1)
        {
            FENode& node = ps->Node(i);
            // set node as having prescribed DOF
            node.set(m_dofEF, e);
        }
    }
}

//-----------------------------------------------------------------------------
//! evaluate the flow rate across this surface
double FEFluidRCRBC::FlowRate()
{
    double Q = 0;
    
    vec3d rt[FEElement::MAX_NODES];
    vec3d vt[FEElement::MAX_NODES];
    
    for (int iel=0; iel<m_psurf->Elements(); ++iel)
    {
        FESurfaceElement& el = m_psurf->Element(iel);
        
        // nr integration points
        int nint = el.GaussPoints();
        
        // nr of element nodes
        int neln = el.Nodes();
        
        // nodal coordinates
        for (int i=0; i<neln; ++i) {
            FENode& node = m_psurf->GetMesh()->Node(el.m_node[i]);
            rt[i] = node.m_rt*m_alpha + node.m_rp*(1-m_alpha);
            vt[i] = node.get_vec3d(m_dofW[0], m_dofW[1], m_dofW[2])*m_alphaf + node.get_vec3d_prev(m_dofW[0], m_dofW[1], m_dofW[2])*(1-m_alphaf);
        }
        
        double* Nr, *Ns;
        double* N;
        double* w  = el.GaussWeights();
        
        vec3d dxr, dxs, v;
        
        // repeat over integration points
        for (int n=0; n<nint; ++n)
        {
            N  = el.H(n);
            Nr = el.Gr(n);
            Ns = el.Gs(n);
            
            // calculate the velocity and tangent vectors at integration point
            dxr = dxs = v = vec3d(0,0,0);
            for (int i=0; i<neln; ++i)
            {
                v += vt[i]*N[i];
                dxr += rt[i]*Nr[i];
                dxs += rt[i]*Ns[i];
            }
            
            vec3d normal = dxr ^ dxs;
            double q = normal*v;
            Q += q*w[n];
        }
    }
    
    return Q;
}

//-----------------------------------------------------------------------------
//! calculate residual
void FEFluidRCRBC::LoadVector(FEGlobalVector& R, const FETimeInfo& tp)
{
    m_alpha = tp.alpha; m_alphaf = tp.alphaf;
    
    /*
    int numsteps = GetFEModel()->GetCurrentStep()->m_ntimesteps;
    m_flowHist.resize(numsteps + 1, 0);
    m_timeHist.resize(numsteps + 1);
    m_stepHist.resize(numsteps + 1);
    
    m_timeHist[numsteps] = tp.currentTime;
    m_stepHist[numsteps] = tp.timeIncrement;
    */
}

//-----------------------------------------------------------------------------
//! serialization
void FEFluidRCRBC::Serialize(DumpStream& ar)
{
    FESurfaceLoad::Serialize(ar);
    ar & m_alpha & m_alphaf;
    ar & m_pfluid;
    ar & m_timeHist;
    ar & m_flowHist;
    ar & m_stepHist;
}
