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
#include "FEFluidSolutesRCRBC.h"
#include "FEBioFluidSolutes.h"
#include <FECore/FEAnalysis.h>
#include <FECore/log.h>
#include <FECore/FEModel.h>

//=============================================================================
BEGIN_FECORE_CLASS(FEFluidSolutesRCRBC, FEPrescribedSurface)
    ADD_PARAMETER(m_R , "R")->setUnits("F.t/L^5");
    ADD_PARAMETER(m_Rd, "Rd")->setUnits("F.t/L^5");
    ADD_PARAMETER(m_p0, "initial_pressure")->setUnits(UNIT_PRESSURE);
    ADD_PARAMETER(m_pd, "pressure_offset")->setUnits(UNIT_PRESSURE);
    ADD_PARAMETER(m_C , "capacitance")->setUnits("L^5/F");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! constructor
FEFluidSolutesRCRBC::FEFluidSolutesRCRBC(FEModel* pfem) : FEPrescribedSurface(pfem), m_dofW(pfem)
{
    m_R = 0.0;
    m_p0 = 0;
    m_Rd = 0.0;
    m_pd = 0.0;
    m_C = 0.0;
    m_dofEF = -1;
    m_dofC = -1;
    m_Rgas = 0;
    m_Tabs = 0;
}

//-----------------------------------------------------------------------------
//! initialize
//! TODO: Generalize to include the initial conditions
bool FEFluidSolutesRCRBC::Init()
{
    m_dofW.AddVariable(FEBioFluidSolutes::GetVariableName(FEBioFluidSolutes::RELATIVE_FLUID_VELOCITY));
    m_dofEF = GetDOFIndex(FEBioFluidSolutes::GetVariableName(FEBioFluidSolutes::FLUID_DILATATION), 0);
    m_dofC = GetDOFIndex(FEBioFluidSolutes::GetVariableName(FEBioFluidSolutes::FLUID_CONCENTRATION), 0);
    SetDOFList(m_dofEF);

    if (FEPrescribedSurface::Init() == false) return false;
    
    m_Rgas = GetFEModel()->GetGlobalConstant("R");
    m_Tabs = GetFEModel()->GetGlobalConstant("T");

    m_e.assign(GetSurface()->Nodes(), 0.0);
    
    m_pn = m_pp = m_p0;
    m_pdn = m_pdp = m_pd;
    m_qn = m_qp = 0;
    m_tp = 0;
    
    return true;
}

//-----------------------------------------------------------------------------
void FEFluidSolutesRCRBC::Update()
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
    
    // prescribe this dilatation at the nodes
    FESurface* ps = GetSurface();
    
    int N = ps->Nodes();
    std::vector<vector<double>> efNodes(N, vector<double>());
    std::vector<vector<double>> caNodes(N, vector<double>());
    
    // Project mean of osmotic coefficient and partition coefficient from int points to nodes on surface
    // Use these to evaluate osmolarity at each node on surface
    for (int i=0; i<ps->Elements(); ++i)
    {
        FESurfaceElement& el = ps->Element(i);
        FEElement* e = el.m_elem[0];
        FEMaterial* pm = GetFEModel()->GetMaterial(e->GetMatID());
        FEFluid* pfl = pm->ExtractProperty<FEFluid>();
        FESoluteInterface* psi = pm->ExtractProperty<FESoluteInterface>();
        FESolidElement* se = dynamic_cast<FESolidElement*>(e);
        if (se) {
            double efo[FEElement::MAX_NODES] = {0};
            if (psi) {
                const int nsol = psi->Solutes();
                std::vector<double> kappa(nsol,0);
                double osc = 0;
                const int nint = se->GaussPoints();
                // get the average osmotic coefficient and partition coefficients in the solid element
                for (int j=0; j<nint; ++j) {
                    FEMaterialPoint* pt = se->GetMaterialPoint(j);
                    osc += psi->GetOsmoticCoefficient()->OsmoticCoefficient(*pt);
                    for (int k=0; k<nsol; ++k)
                        kappa[k] += psi->GetPartitionCoefficient(*pt, k);
                }
                osc /= nint;
                for (int k=0; k<nsol; ++k) kappa[k] /= nint;
                // loop over face nodes
                for (int j=0; j<el.Nodes(); ++j) {
                    double osm = 0;
                    FENode& node = ps->Node(el.m_lnode[j]);
                    // calculate osmolarity at this node, using nodal effective solute concentrations
                    for (int k=0; k<nsol; ++k)
                        osm += node.get(m_dofC+psi->GetSolute(k)->GetSoluteID()-1)*kappa[k];
                    // evaluate dilatation at this node
                    double c = m_Rgas*osc*osm;
                    bool good = pfl->Dilatation(0, m_pn - m_Rgas*m_Tabs*osc*osm, efo[j]);
                    assert(good);
                }
            }
            else {
                // loop over face nodes
                for (int j=0; j<el.Nodes(); ++j) {
                    FENode& node = ps->Node(el.m_lnode[j]);
                    // evaluate dilatation at this node
                    bool good = pfl->Dilatation(0, m_pn, efo[j]);
                    assert(good);
                }
            }
            // only keep the dilatations at the nodes of the surface face
            for (int j=0; j<el.Nodes(); ++j)
                efNodes[el.m_lnode[j]].push_back(efo[j]);
        }
        //If no solid element, insert all 0s
        else {
            for (int j=0; j<el.Nodes(); ++j)
                efNodes[el.m_lnode[j]].push_back(0);
        }
    }
    
    //For each node, average the nodal ef
    for (int i=0; i<ps->Nodes(); ++i)
    {
        double ef = 0;
        for (int j = 0; j < efNodes[i].size(); ++j)
            ef += efNodes[i][j];
        ef /= efNodes[i].size();
        
        // store value for now
        m_e[i] = ef;
    }
    
    FEPrescribedSurface::Update();
}

//-----------------------------------------------------------------------------
//! evaluate the flow rate across this surface at current time
double FEFluidSolutesRCRBC::FlowRate()
{
    double Q = 0;
    
    vec3d rt[FEElement::MAX_NODES];
    vec3d vt[FEElement::MAX_NODES];
    
    const FETimeInfo& tp = GetTimeInfo();
    double alpha = tp.alpha;
    double alphaf = tp.alphaf;
    
    // prescribe this dilatation at the nodes
    FESurface* ps = GetSurface();
    
    for (int iel=0; iel<ps->Elements(); ++iel)
    {
        FESurfaceElement& el = ps->Element(iel);
        
        // nr integration points
        int nint = el.GaussPoints();
        
        // nr of element nodes
        int neln = el.Nodes();
        
        // nodal coordinates
        for (int i=0; i<neln; ++i) {
            FENode& node = ps->Node(el.m_lnode[i]);
            rt[i] = node.m_rt;
            vt[i] = node.get_vec3d(m_dofW[0], m_dofW[1], m_dofW[2]);
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
void FEFluidSolutesRCRBC::GetNodalValues(int nodelid, std::vector<double>& val)
{
    val[0] = m_e[nodelid];
}

//-----------------------------------------------------------------------------
// copy data from another class
void FEFluidSolutesRCRBC::CopyFrom(FEBoundaryCondition* pbc)
{
    // TODO: implement this
    assert(false);
}

//-----------------------------------------------------------------------------
//! serialization
void FEFluidSolutesRCRBC::Serialize(DumpStream& ar)
{
    FEPrescribedSurface::Serialize(ar);
    ar & m_e;
    ar & m_pn & m_pp & m_qn & m_qp & m_pdn & m_pdp & m_tp & m_e;
    if (ar.IsShallow()) return;
    ar & m_dofW & m_dofEF & m_dofC;
    ar & m_Rgas & m_Tabs;
}
