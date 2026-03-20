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
#include "FEFluidCOBC.h"
#include "FEFluid.h"
#include "FEBioFluid.h"
#include <FECore/FEAnalysis.h>
#include <FECore/log.h>
#include <FECore/FEModel.h>

//=============================================================================
BEGIN_FECORE_CLASS(FEFluidCOBC, FEPrescribedSurface)
    ADD_PARAMETER(m_Ra, "Ra")->setLongName("arterial resistance")->setUnits("F.t/L^5");
    ADD_PARAMETER(m_Ca, "Ca")->setLongName("arterial compliance")->setUnits("L^5/F");
    ADD_PARAMETER(m_Ram,"Ra-micro")->setLongName("microvascular arterial resistance")->setUnits("F.t/L^5");
    ADD_PARAMETER(m_Cmi,"Cmi")->setLongName("myocardial compliance")->setUnits("L^5/F");
    ADD_PARAMETER(m_Rv, "Rv")->setLongName("ventricular resistance")->setUnits("F.t/L^5");
    ADD_PARAMETER(m_Rvm,"Rv-micro")->setLongName("microvascular ventricular resistance")->setUnits("F.t/L^5");
    ADD_PARAMETER(m_p0, "initial_pressure")->setUnits(UNIT_PRESSURE);
//    ADD_PARAMETER(m_pd, "pressure_offset")->setUnits(UNIT_PRESSURE);

    ADD_PROPERTY(m_Pmi,"Pmi",FEProperty::Optional)->SetLongName("mycoardial pressure");
    ADD_PROPERTY(m_PRA,"PRA",FEProperty::Optional)->SetLongName("right atrium pressure");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! constructor
FEFluidCOBC::FEFluidCOBC(FEModel* pfem) : FEPrescribedSurface(pfem), m_dofW(pfem)
{
    m_Ra = m_Ram = m_Rv = m_Rvm = 0.0;
    m_Ca = m_Cmi = 0.0;
    m_pfluid = nullptr;
    m_psurf = nullptr;
    m_p0 = 0;
    m_pd = 0.0;
    m_e = 0.0;
    
 }

//-----------------------------------------------------------------------------
//! initialize
//! TODO: Generalize to include the initial conditions
bool FEFluidCOBC::Init()
{
    m_dofW.AddVariable(FEBioFluid::GetVariableName(FEBioFluid::RELATIVE_FLUID_VELOCITY));
    m_dofEF = GetDOFIndex(FEBioFluid::GetVariableName(FEBioFluid::FLUID_DILATATION), 0);
    SetDOFList(m_dofEF);

    if (FEPrescribedSurface::Init() == false) return false;
    
    m_psurf = GetSurface();
    
    // get fluid from first surface element
    // assuming the entire surface bounds the same fluid
    FESurfaceElement& el = m_psurf->Element(0);
    FEElement* pe = el.m_elem[0].pe;
    if (pe == nullptr) return false;
    
    // get the material
    FEMaterial* pm = GetFEModel()->GetMaterial(pe->GetMatID());
    m_pfluid = pm->ExtractProperty<FEFluidMaterial>();
    if (m_pfluid == nullptr) return false;
    
    m_pn = m_pp = m_ppp = m_p0;
//    m_pdn = m_pdp = m_pd;
    m_qn = m_qp = m_qpp = 0;
    m_tp = m_tpp = 0;

    if (m_Pmi) { if (!m_Pmi->Init()) return false; }
    if (m_PRA) { if (!m_PRA->Init()) return false; }
    
    return true;
}

//-----------------------------------------------------------------------------
//! Evaluate and prescribe the resistance pressure
void FEFluidCOBC::UpdateDilatation()
{
	// Check if we started a new time, if so, update variables
	FETimeInfo& timeInfo = GetFEModel()->GetTime();
	double time = timeInfo.currentTime;
	int iter = timeInfo.currentIteration;
	double dt = timeInfo.timeIncrement;
	if ((time > m_tp) && (iter == 0)) {
        m_ppp = m_pp;
		m_pp = m_pn;
        m_qpp = m_qp;
		m_qp = m_qn;
		m_pdp = m_pdn;
        m_tpp = m_tp;
		m_tp = time;
	}

	// evaluate the flow rate at the current time
	m_qn = FlowRate();
	m_pdn = m_pd;
    
    double Rve   = m_Rv + m_Rvm;
    double Rae   = m_Ram + Rve;
    double Rt    = m_Ra + Rae;
    double tauvi = Rve*m_Cmi;
    double taua  = m_Ram*m_Ca;
    double tauae = Rae*m_Ca;
    
    double dtp = m_tp - m_tpp;
    
    double c0 = tauvi*taua/(dt*dt);
    double c1 = m_Ra*c0;
    double c2 = (m_Ra*(tauae + tauvi)+m_Ram*tauvi)/dt;
    double a = c1 + c2 + Rt;
    double b = (c1+c2)*m_qp + m_Ra*tauvi*taua/dt/dtp*(m_qp-m_qpp);
    double c = c0 + (tauvi+tauae)/dt + 1;
    double d = c0*m_pp + tauvi*taua/dt/dtp*(m_pp - m_ppp) + (tauvi+tauae)/dt*m_pp;
    if (m_PRA) d += m_PRA->value(time);
    if (m_Pmi) d += tauvi*m_Pmi->derive(time);

	// calculate the RCR pressure
    m_pn = (a*m_qn + d - b)/c;

	// calculate the dilatation
	m_e = 0.0;
	bool good = m_pfluid->Dilatation(0, m_pn, m_e);
	assert(good);
}

void FEFluidCOBC::UpdateModel() { Update(); }

void FEFluidCOBC::Update()
{
	UpdateDilatation();

    // the base class handles mapping the values to the nodal dofs
    FEPrescribedSurface::Update();

	// TODO: Is this necessary?
	GetFEModel()->SetMeshUpdateFlag(true);
}

//-----------------------------------------------------------------------------
//! evaluate the flow rate across this surface at current time
double FEFluidCOBC::FlowRate()
{
    double Q = 0;
    
    const FETimeInfo& tp = GetTimeInfo();

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
            FENode& node = m_psurf->Node(el.m_lnode[i]);
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
void FEFluidCOBC::PrepStep(std::vector<double>& ui, bool brel)
{
	UpdateDilatation();
	FEPrescribedSurface::PrepStep(ui, brel);
}

//-----------------------------------------------------------------------------
void FEFluidCOBC::GetNodalValues(int nodelid, std::vector<double>& val)
{
    val[0] = m_e;
	FENode& node = GetMesh().Node(m_nodeList[nodelid]);
	node.set(m_dofEF, m_e);
}

//-----------------------------------------------------------------------------
// copy data from another class
void FEFluidCOBC::CopyFrom(FEBoundaryCondition* pbc)
{
    // TODO: implement this
    assert(false);
}

//-----------------------------------------------------------------------------
//! serialization
void FEFluidCOBC::Serialize(DumpStream& ar)
{
    FEPrescribedSurface::Serialize(ar);
    ar & m_pn & m_pp & m_pp & m_qn & m_qp & m_qpp & m_tp & m_tpp & m_e;
    if (ar.IsShallow()) return;
    ar & m_pfluid;
    ar & m_dofW & m_dofEF;
    ar & m_psurf;
    ar & m_Pmi & m_PRA;
}
