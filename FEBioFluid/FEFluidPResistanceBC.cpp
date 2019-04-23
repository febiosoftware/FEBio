/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, Columbia University, and others.

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
#include "FEFluidPResistanceBC.h"
#include "FEFluid.h"
#include "FEFluidP.h"
#include "FEFluidFSI.h"
#include "stdafx.h"
#include "FECore/FEModel.h"
#include "FECore/DOFS.h"
#include "FEBioFluid.h"

//=============================================================================
BEGIN_FECORE_CLASS(FEFluidPResistanceBC, FESurfaceLoad)
	ADD_PARAMETER(m_R , "R");
	ADD_PARAMETER(m_p0, "pressure_offset");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! constructor
FEFluidPResistanceBC::FEFluidPResistanceBC(FEModel* pfem) : FESurfaceLoad(pfem), m_dofW(pfem)
{
    m_R = 0.0;
    m_alpha = 1.0;
    m_p0 = 0;
    
	m_dofW.AddVariable(FEBioFluid::GetVariableName(FEBioFluid::RELATIVE_FLUID_VELOCITY));
    m_dofEF  = pfem->GetDOFIndex(FEBioFluid::GetVariableName(FEBioFluid::FLUID_DILATATION), 0);
}

//-----------------------------------------------------------------------------
//! initialize
bool FEFluidPResistanceBC::Init()
{
    FEModelComponent::Init();
    
    FESurface* ps = &GetSurface();
    ps->Init();
    
    return true;
}

//-----------------------------------------------------------------------------
//! Activate the degrees of freedom for this BC
void FEFluidPResistanceBC::Activate()
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
void FEFluidPResistanceBC::Update()
{
    // evaluate the flow rate
    double Q = FlowRate();
    
    // calculate the resistance pressure
    double p = m_R*Q + m_p0;
    
    // prescribe this pressure at the nodes
    FESurface* ps = &GetSurface();
    
    for (int i=0; i<ps->Nodes(); ++i)
    {
        if (ps->Node(i).m_ID[m_dofEF] < -1)
        {
            FENode& node = ps->Node(i);
            // set node as having prescribed DOF
            node.set(m_dofEF, p);
        }
    }
}

//-----------------------------------------------------------------------------
//! evaluate the flow rate across this surface
double FEFluidPResistanceBC::FlowRate()
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

//! calculate residual
void FEFluidPResistanceBC::Residual(FEGlobalVector& R, const FETimeInfo& tp) 
{ 
	m_alpha = tp.alpha; m_alphaf = tp.alphaf; 
}
