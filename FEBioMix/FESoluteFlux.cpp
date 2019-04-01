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
#include "FESoluteFlux.h"
#include "FECore/FEModel.h"
#include "FECore/FEAnalysis.h"

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FESoluteFlux, FESurfaceLoad)
	ADD_PARAMETER(m_flux   , "flux");
	ADD_PARAMETER(m_blinear, "linear");
    ADD_PARAMETER(m_bshellb, "shell_bottom");
	ADD_PARAMETER(m_isol   , "solute_id");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! constructor
FESoluteFlux::FESoluteFlux(FEModel* pfem) : FESurfaceLoad(pfem)
{ 
	m_flux = 1.0;
	m_blinear = false; 
    m_bshellb = false;
	m_isol = 0;

	m_dofX = pfem->GetDOFIndex("x");
	m_dofY = pfem->GetDOFIndex("y");
	m_dofZ = pfem->GetDOFIndex("z");
	m_dofC = pfem->GetDOFIndex("concentration", 0);
    m_dofSX = pfem->GetDOFIndex("sx");
    m_dofSY = pfem->GetDOFIndex("sy");
    m_dofSZ = pfem->GetDOFIndex("sz");
    m_dofD = pfem->GetDOFIndex("shell concentration", 0);
}
	
//-----------------------------------------------------------------------------
//! allocate storage
void FESoluteFlux::SetSurface(FESurface* ps)
{ 
	FESurfaceLoad::SetSurface(ps);
	m_flux.SetItemList(ps->GetFacetSet());
}

//-----------------------------------------------------------------------------
void FESoluteFlux::UnpackLM(FEElement& el, vector<int>& lm)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

    // get nodal DOFS
    DOFS& fedofs = fem.GetDOFS();
    int MAX_CDOFS = fedofs.GetVariableSize("concentration");
    
	int N = el.Nodes();
	lm.resize(N*(3+MAX_CDOFS));

    if (!m_bshellb) {
        for (int i=0; i<N; ++i)
        {
            int n = el.m_node[i];
            
            FENode& node = mesh.Node(n);
            vector<int>& id = node.m_ID;
            
            // first the displacement dofs
            lm[3*i  ] = id[m_dofX];
            lm[3*i+1] = id[m_dofY];
            lm[3*i+2] = id[m_dofZ];
            
            // concentration dofs
            for (int k=0; k<MAX_CDOFS; ++k)
                lm[(3+k)*N + i] = id[m_dofC+k];
        }
    }
    else {
        for (int i=0; i<N; ++i)
        {
            int n = el.m_node[i];
            
            FENode& node = mesh.Node(n);
            vector<int>& id = node.m_ID;
            
            // first the displacement dofs
            lm[3*i  ] = id[m_dofSX];
            lm[3*i+1] = id[m_dofSY];
            lm[3*i+2] = id[m_dofSZ];
            
            // concentration dofs
            for (int k=0; k<MAX_CDOFS; ++k)
                lm[(3+k)*N + i] = id[m_dofD+k];
        }
    }
}

//-----------------------------------------------------------------------------
//! calculates the stiffness contribution due to solute flux
//!
void FESoluteFlux::FluxStiffness(FESurfaceElement& el, matrix& ke, double dt)
{
	int i, j, n;
	
	int nint = el.GaussPoints();
	int neln = el.Nodes();
	
	vec3d dxr, dxs;
	
	// gauss weights
	double* w = el.GaussWeights();
	
	// get the element's nodal coordinates
	vec3d rt[FEElement::MAX_NODES];
	for (j=0; j<neln; ++j) rt[j] = m_psurf->GetMesh()->Node(el.m_node[j]).m_rt;
    if (m_bshellb) {
        for (int j=0; j<neln; ++j) {
            FENode& nd = m_psurf->GetMesh()->Node(el.m_node[j]);
            rt[j] -= nd.m_d0 + nd.get_vec3d(m_dofX, m_dofY, m_dofZ) - nd.get_vec3d(m_dofSX, m_dofSY, m_dofSZ);
            
        }
    }
	
	vec3d kab, t1, t2;
	
	ke.zero();
	
	double* N, *Gr, *Gs;
	
	// repeat over integration points
	for (n=0; n<nint; ++n)
	{
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);

		N = el.H(n);
		Gr = el.Gr(n);
		Gs = el.Gs(n);
		
		// calculate velocities and covariant basis vectors at integration point
		dxr = dxs = vec3d(0,0,0);
		for (i=0; i<neln; ++i)
		{
			dxr += rt[i]*Gr[i];
			dxs += rt[i]*Gs[i];
		}

		double wr = m_flux(mp);
        if (m_bshellb) wr = -wr;
		
		// calculate surface normal
		vec3d dxt = dxr ^ dxs;
		
		// calculate stiffness component
		for (i=0; i<neln; ++i)
			for (j=0; j<neln; ++j)
			{
				t1 = dxt/dxt.norm()*wr;
				t2 = dxs*Gr[j] - dxr*Gs[j];
				kab = (t1 ^ t2)*(N[i]*w[n])*dt;
				
				ke[4*i+3][4*j  ] += kab.x;
				ke[4*i+3][4*j+1] += kab.y;
				ke[4*i+3][4*j+2] += kab.z;
			}
	}
}

//-----------------------------------------------------------------------------
//! calculates the equivalent nodal volumetric flow rates due to solute flux
//!
bool FESoluteFlux::FlowRate(FESurfaceElement& el, vector<double>& fe, double dt)
{
	int i, n;
	
	// nr integration points
	int nint = el.GaussPoints();
	
	// nr of element nodes
	int neln = el.Nodes();
	
	// get the element's nodal coordinates
	vec3d rt[FEElement::MAX_NODES];
	for (int j=0; j<neln; ++j) rt[j] = m_psurf->GetMesh()->Node(el.m_node[j]).m_rt;
    if (m_bshellb) {
        for (int j=0; j<neln; ++j) {
            FENode& nd = m_psurf->GetMesh()->Node(el.m_node[j]);
            rt[j] -= nd.m_d0 + nd.get_vec3d(m_dofX, m_dofY, m_dofZ) - nd.get_vec3d(m_dofSX, m_dofSY, m_dofSZ);
            
        }
    }
	
	double* Gr, *Gs;
	double* N;
	double* w  = el.GaussWeights();
	
	vec3d dxr, dxs, dxt;
	
	// volumetric flow rate
	double f;
	
	// repeat over integration points
	zero(fe);
	for (n=0; n<nint; ++n)
	{
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);

		N  = el.H(n);
		Gr = el.Gr(n);
		Gs = el.Gs(n);
		
		dxr = dxs = vec3d(0,0,0);
		for (i=0; i<neln; ++i) 
		{
			dxr += rt[i]*Gr[i];
			dxs += rt[i]*Gs[i];
		}
		double wr = m_flux(mp);
        if (m_bshellb) wr = -wr;
		dxt = dxr ^ dxs;
		
		f = dxt.norm()*wr*w[n]*dt;
		
		for (i=0; i<neln; ++i)
		{
			fe[i] += N[i]*f;
		}
	}
	
	return true;
}

//-----------------------------------------------------------------------------
//! calculates the equivalent nodal volumetric flow rates due to solute flux
//!
bool FESoluteFlux::LinearFlowRate(FESurfaceElement& el, vector<double>& fe, double dt)
{
	int i, n;
	
	// nr integration points
	int nint = el.GaussPoints();
	
	// nr of element nodes
	int neln = el.Nodes();
	
	// nodal coordinates
	vec3d r0[FEElement::MAX_NODES];
	for (i=0; i<neln; ++i) r0[i] = m_psurf->GetMesh()->Node(el.m_node[i]).m_r0;
    if (m_bshellb) {
        for (int j=0; j<neln; ++j) {
            FENode& nd = m_psurf->GetMesh()->Node(el.m_node[j]);
            r0[j] -= nd.m_d0;
        }
    }
	
	double* Gr, *Gs;
	double* N;
	double* w  = el.GaussWeights();
	
	vec3d dxr, dxs, dxt, vs;
	
	// volumetric flow rate
	double f;
	
	// repeat over integration points
	zero(fe);
	for (n=0; n<nint; ++n)
	{
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);

		N  = el.H(n);
		Gr = el.Gr(n);
		Gs = el.Gs(n);
		
		dxr = dxs = vec3d(0,0,0);
		for (i=0; i<neln; ++i) 
		{
			dxr += r0[i]*Gr[i];
			dxs += r0[i]*Gs[i];
		}
		double wr = m_flux(mp);

        if (m_bshellb) wr = -wr;
		dxt = dxr ^ dxs;
		
		f = dxt.norm()*wr*w[n]*dt;
		
		for (i=0; i<neln; ++i)
		{
			fe[i] += N[i]*f;
		}
	}
	
	return true;
}

//-----------------------------------------------------------------------------
void FESoluteFlux::StiffnessMatrix(const FETimeInfo& tp, FESolver* psolver)
{
	double dt = tp.timeIncrement;
	
	matrix ke;
	vector<int> elm;
	
	int N = m_psurf->Elements();
	for (int m=0; m<N; ++m)
	{
		// get the surface element
		FESurfaceElement& el = m_psurf->Element(m);
			
		// calculate nodal normal solute flux
		int neln = el.Nodes();
		if (m_blinear == false)
		{
					
			// get the element stiffness matrix
			int ndof = neln*4;
			ke.resize(ndof, ndof);
					
			// calculate pressure stiffness
			FluxStiffness(el, ke, dt);

			// get the element's LM vector
			UnpackLM(el, elm);
					
			// TODO: the problem here is that the LM array that is returned by the UnpackElement
			// function does not give the equation numbers in the right order. For this reason we
			// have to create a new lm array and place the equation numbers in the right order.
			// What we really ought to do is fix the UnpackElement function so that it returns
			// the LM vector in the right order for solute-solid elements.
			vector<int> lm(ndof);
			for (int i=0; i<neln; ++i)
			{
				lm[4*i  ] = elm[3*i];
				lm[4*i+1] = elm[3*i+1];
				lm[4*i+2] = elm[3*i+2];
				lm[4*i+3] = elm[(3+m_isol-1)*neln+i];  // m_isol is 1-based
			}
					
			// assemble element matrix in global stiffness matrix
			psolver->AssembleStiffness(el.m_node, lm, ke);
		}
	}
}

//-----------------------------------------------------------------------------
void FESoluteFlux::Residual(const FETimeInfo& tp, FEGlobalVector& R)
{
	double dt = tp.timeIncrement;
	
	vector<double> fe;
	vector<int> elm;

	int N = m_psurf->Elements();
	for (int i=0; i<N; ++i)
	{
		FESurfaceElement& el = m_psurf->Element(i);
			
		// calculate nodal normal solute flux
		int neln = el.Nodes();
		
		int ndof = neln;
		fe.resize(ndof);
			
		if (m_blinear) LinearFlowRate(el, fe, dt); else FlowRate(el, fe, dt);

		// get the element's LM vector
		UnpackLM(el, elm);

		// We only need the solute concentration dofs, so just extract these.
		vector<int> lm(ndof);
		for (int i=0; i<neln; ++i)
			lm[i] = elm[(3+m_isol-1)*neln+i];  // m_isol is 1-based
			
		// add element force vector to global force vector
		R.Assemble(el.m_node, lm, fe);
	}
}
