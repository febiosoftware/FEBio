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
#include "FEFluidFlux.h"
#include "FECore/FESolver.h"
#include "FECore/FEModel.h"
#include "FECore/FEAnalysis.h"

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEFluidFlux, FESurfaceLoad)
	ADD_PARAMETER(m_flux    , "flux"   );
	ADD_PARAMETER(m_blinear , "linear" );
    ADD_PARAMETER(m_bshellb , "shell_bottom");
	ADD_PARAMETER(m_bmixture, "mixture");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEFluidFlux::FEFluidFlux(FEModel* pfem) : FESurfaceLoad(pfem)
{ 
	m_blinear = false; 
	m_bmixture = false;
    m_bshellb = false;
	m_flux = 1.0;

	// get the degrees of freedom
	m_dofX = pfem->GetDOFIndex("x");
	m_dofY = pfem->GetDOFIndex("y");
	m_dofZ = pfem->GetDOFIndex("z");
	m_dofP = pfem->GetDOFIndex("p");
	m_dofVX = pfem->GetDOFIndex("vx");
	m_dofVY = pfem->GetDOFIndex("vy");
	m_dofVZ = pfem->GetDOFIndex("vz");
    m_dofSX = pfem->GetDOFIndex("sx");
    m_dofSY = pfem->GetDOFIndex("sy");
    m_dofSZ = pfem->GetDOFIndex("sz");
    m_dofQ = pfem->GetDOFIndex("q");
    m_dofSVX = pfem->GetDOFIndex("svx");
    m_dofSVY = pfem->GetDOFIndex("svy");
    m_dofSVZ = pfem->GetDOFIndex("svz");
}

//-----------------------------------------------------------------------------
void FEFluidFlux::SetSurface(FESurface* ps) 
{ 
	FESurfaceLoad::SetSurface(ps);
	m_flux.SetItemList(ps->GetFacetSet());
}

//-----------------------------------------------------------------------------
void FEFluidFlux::UnpackLM(FEElement& el, vector<int>& lm)
{
	FEMesh& mesh = GetFEModel()->GetMesh();
	int N = el.Nodes();
	lm.resize(N*4);
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
            
            // now the pressure dofs
            lm[3*N+i] = id[m_dofP];
        }
    }
    else {
        for (int i=0; i<N; ++i)
        {
            int n = el.m_node[i];
            FENode& node = mesh.Node(n);
            vector<int>& id = node.m_ID;
            
            // first the shell displacement dofs
            lm[3*i  ] = id[m_dofSX];
            lm[3*i+1] = id[m_dofSY];
            lm[3*i+2] = id[m_dofSZ];
            
            // now the shell pressure dofs
            lm[3*N+i] = id[m_dofQ];
        }
    }
}

//-----------------------------------------------------------------------------
//! calculates the stiffness contribution due to fluid flux

void FEFluidFlux::FluxStiffness(FESurfaceElement& el, matrix& ke, double dt, bool mixture)
{
	int i, j, n;

	int nint = el.GaussPoints();
	int neln = el.Nodes();

	// solid velocity at integration point
	vec3d vr;

	vec3d dxr, dxs, dxt;

	// gauss weights
	double* w = el.GaussWeights();

	// nodal coordinates and velocities
	vec3d rt[FEElement::MAX_NODES], vt[FEElement::MAX_NODES];
	for (i=0; i<neln; ++i)
	{
		rt[i] = m_psurf->GetMesh()->Node(el.m_node[i]).m_rt;
		vt[i] = m_psurf->GetMesh()->Node(el.m_node[i]).get_vec3d(m_dofVX, m_dofVY, m_dofVZ);
	}
    if (m_bshellb) {
        for (int j=0; j<neln; ++j) {
            FENode& nd = m_psurf->GetMesh()->Node(el.m_node[j]);
            rt[j] -= nd.m_d0 + nd.get_vec3d(m_dofX, m_dofY, m_dofZ) - nd.get_vec3d(m_dofSX, m_dofSY, m_dofSZ);
            vt[i] = nd.get_vec3d(m_dofSVX, m_dofSVY, m_dofSVZ);
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
		vr = dxr = dxs = vec3d(0,0,0);
		for (i=0; i<neln; ++i)
		{
			vr += vt[i]*N[i];
			dxr += rt[i]*Gr[i];
			dxs += rt[i]*Gs[i];
		}

		double wr = m_flux(mp);
        if (m_bshellb) wr = -wr;
		
		// calculate surface normal
		dxt = dxr ^ dxs;

		// calculate stiffness component
		for (i=0; i<neln; ++i)
			for (j=0; j<neln; ++j)
			{
				t1 = (dxt/dxt.norm())*wr - vr*mixture;
				t2 = dxs*Gr[j] - dxr*Gs[j];
				kab = ((t1^t2)*(!mixture) + dxt*mixture*N[j]/dt)*N[i]*w[n]*dt;
//				kab = (t1^t2 + dxt*mixture*N[j]/dt)*N[i]*w[n]*dt;

				ke[3*neln+i][3*j  ] += kab.x;
				ke[3*neln+i][3*j+1] += kab.y;
				ke[3*neln+i][3*j+2] += kab.z;
			}
	}
}

//-----------------------------------------------------------------------------
//! calculates the stiffness contribution due to fluid flux
//! for a steady-state analysis

void FEFluidFlux::FluxStiffnessSS(FESurfaceElement& el, matrix& ke, double dt, bool mixture)
{
	int i, j, n;
	
	int nint = el.GaussPoints();
	int neln = el.Nodes();
	
	vec3d dxr, dxs, dxt;
	
	// gauss weights
	double* w = el.GaussWeights();
	
	// nodal coordinates and velocities
	vec3d rt[FEElement::MAX_NODES];
	for (i=0; i<neln; ++i)
		rt[i] = m_psurf->GetMesh()->Node(el.m_node[i]).m_rt;
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
		dxt = dxr ^ dxs;
		
		// calculate stiffness component
		for (i=0; i<neln; ++i)
			for (j=0; j<neln; ++j)
			{
				t1 = (dxt/dxt.norm())*wr;
				t2 = dxs*Gr[j] - dxr*Gs[j];
				kab = (t1^t2)*((!mixture)*N[i]*w[n]*dt);
				
				ke[3*neln+i][3*j  ] += kab.x;
				ke[3*neln+i][3*j+1] += kab.y;
				ke[3*neln+i][3*j+2] += kab.z;
			}
	}
}

//-----------------------------------------------------------------------------
//! calculates the equivalent nodal volumetric flow rates due to fluid flux

bool FEFluidFlux::FlowRate(FESurfaceElement& el, vector<double>& fe, double dt, bool mixture)
{
	int i, n;

	// nr integration points
	int nint = el.GaussPoints();

	// nr of element nodes
	int neln = el.Nodes();

	// nodal coordinates and velocities
	vec3d rt[FEElement::MAX_NODES], vt[FEElement::MAX_NODES];
	for (i=0; i<neln; ++i)
	{
		rt[i] = m_psurf->GetMesh()->Node(el.m_node[i]).m_rt;
		vt[i] = m_psurf->GetMesh()->Node(el.m_node[i]).get_vec3d(m_dofVX, m_dofVY, m_dofVZ);
	}
    if (m_bshellb) {
        for (int j=0; j<neln; ++j) {
            FENode& nd = m_psurf->GetMesh()->Node(el.m_node[j]);
            rt[j] -= nd.m_d0 + nd.get_vec3d(m_dofX, m_dofY, m_dofZ) - nd.get_vec3d(m_dofSX, m_dofSY, m_dofSZ);
            vt[i] = nd.get_vec3d(m_dofSVX, m_dofSVY, m_dofSVZ);
        }
    }
	
	double* Gr, *Gs;
	double* N;
	double* w  = el.GaussWeights();

	// solid velocity at integration point
	vec3d vr;
	
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

		vr = dxr = dxs = vec3d(0,0,0);
		for (i=0; i<neln; ++i) 
		{
			vr += vt[i]*N[i];
			dxr += rt[i]*Gr[i];
			dxs += rt[i]*Gs[i];
		}

		double wr = m_flux(mp);
        if (m_bshellb) wr = -wr;
		dxt = dxr ^ dxs;

		f = (dxt.norm()*wr - (vr*dxt)*mixture)*w[n]*dt;

		for (i=0; i<neln; ++i)
		{
			fe[3*neln+i] += N[i]*f;
		}
	}

	return true;
}

//-----------------------------------------------------------------------------
//! calculates the equivalent nodal volumetric flow rates due to fluid flux
//! for a steady-state analysis

bool FEFluidFlux::FlowRateSS(FESurfaceElement& el, vector<double>& fe, double dt, bool mixture)
{
	int i, n;
	
	// nr integration points
	int nint = el.GaussPoints();
	
	// nr of element nodes
	int neln = el.Nodes();
	
	// nodal coordinates and velocities
	vec3d rt[FEElement::MAX_NODES];
	for (i=0; i<neln; ++i)
		rt[i] = m_psurf->GetMesh()->Node(el.m_node[i]).m_rt;
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
		
		f = (dxt.norm()*wr)*w[n]*dt;
		
		for (i=0; i<neln; ++i)
		{
			fe[3*neln+i] += N[i]*f;
		}
	}
	
	return true;
}

//-----------------------------------------------------------------------------
//! calculates the equivalent nodal volumetric flow rates due to fluid flux

bool FEFluidFlux::LinearFlowRate(FESurfaceElement& el, vector<double>& fe, double dt, bool mixture)
{
	int i, n;

	// nr integration points
	int nint = el.GaussPoints();

	// nr of element nodes
	int neln = el.Nodes();

	FEMesh& mesh = *m_psurf->GetMesh();

	// nodal coordinates and velocity
	vec3d r0[FEElement::MAX_NODES], rt[FEElement::MAX_NODES], vt[FEElement::MAX_NODES];
	for (i=0; i<neln; ++i)
	{
		r0[i] = mesh.Node(el.m_node[i]).m_r0;
		rt[i] = mesh.Node(el.m_node[i]).m_rt;
		vt[i] = mesh.Node(el.m_node[i]).get_vec3d(m_dofVX, m_dofVY, m_dofVZ);
	}
    if (m_bshellb) {
        for (int j=0; j<neln; ++j) {
            FENode& nd = mesh.Node(el.m_node[j]);
            r0[j] -= nd.m_d0;
            rt[j] -= nd.m_d0 + nd.get_vec3d(m_dofX, m_dofY, m_dofZ) - nd.get_vec3d(m_dofSX, m_dofSY, m_dofSZ);
            vt[i] = nd.get_vec3d(m_dofSVX, m_dofSVY, m_dofSVZ);
        }
    }

	double* Gr, *Gs;
	double* N;
	double* w  = el.GaussWeights();

	// solid velocity at integration points
	vec3d vr;
	
	vec3d dXr, dXs, dXt;
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

		dXr = dXs = vec3d(0,0,0);
		vr = dxr = dxs = vec3d(0,0,0);
		for (i=0; i<neln; ++i) 
		{
			vr += vt[i]*N[i];
			dXr += r0[i]*Gr[i];
			dXs += r0[i]*Gs[i];
			dxr += rt[i]*Gr[i];
			dxs += rt[i]*Gs[i];
		}

		double wr = m_flux(mp);
        if (m_bshellb) wr = -wr;
		dXt = dXr ^ dXs;
		dxt = dxr ^ dxs;
		
		f = (dXt.norm()*wr - (vr*dxt)*mixture)*w[n]*dt;

		for (i=0; i<neln; ++i)
		{
			fe[3*neln+i] += N[i]*f;
		}
	}

	return true;
}

//-----------------------------------------------------------------------------
//! calculates the equivalent nodal volumetric flow rates due to fluid flux
//! for a steady-state analysis

bool FEFluidFlux::LinearFlowRateSS(FESurfaceElement& el, vector<double>& fe, double dt, bool mixture)
{
	int i, n;
	
	// nr integration points
	int nint = el.GaussPoints();
	
	// nr of element nodes
	int neln = el.Nodes();
	
	FEMesh& mesh = *m_psurf->GetMesh();
	
	// nodal coordinates and velocity
	vec3d r0[FEElement::MAX_NODES], rt[FEElement::MAX_NODES];
	for (i=0; i<neln; ++i)
	{
		r0[i] = mesh.Node(el.m_node[i]).m_r0;
		rt[i] = mesh.Node(el.m_node[i]).m_rt;
	}
    if (m_bshellb) {
        for (int j=0; j<neln; ++j) {
            FENode& nd = mesh.Node(el.m_node[j]);
            r0[j] -= nd.m_d0;
            rt[j] -= nd.m_d0 + nd.get_vec3d(m_dofX, m_dofY, m_dofZ) - nd.get_vec3d(m_dofSX, m_dofSY, m_dofSZ);
        }
    }
	
	double* Gr, *Gs;
	double* N;
	double* w  = el.GaussWeights();
	
	vec3d dXr, dXs, dXt;
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
		
		dXr = dXs = vec3d(0,0,0);
		dxr = dxs = vec3d(0,0,0);
		for (i=0; i<neln; ++i) 
		{
			dXr += r0[i]*Gr[i];
			dXs += r0[i]*Gs[i];
			dxr += rt[i]*Gr[i];
			dxs += rt[i]*Gs[i];
		}

		double wr = m_flux(mp);
        if (m_bshellb) wr = -wr;
		dXt = dXr ^ dXs;
		dxt = dxr ^ dxs;
		
		f = (dXt.norm()*wr)*w[n]*dt;
		
		for (i=0; i<neln; ++i)
		{
			fe[3*neln+i] += N[i]*f;
		}
	}
	
	return true;
}

//-----------------------------------------------------------------------------
void FEFluidFlux::StiffnessMatrix(const FETimeInfo& tp, FESolver* psolver)
{
	FEModel& fem = *psolver->GetFEModel();
	double dt = tp.timeIncrement;

	matrix ke;

	vector<int> elm;

	int N = m_psurf->Elements();

	if (fem.GetCurrentStep()->m_nanalysis == FE_STEADY_STATE)
	{
		for (int m=0; m<N; ++m)
		{
			// get the surface element
			FESurfaceElement& el = m_psurf->Element(m);
			UnpackLM(el, elm);
				
			// calculate nodal normal fluid flux
			int neln = el.Nodes();
					
			if (!m_blinear || m_bmixture)
			{
					
				// get the element stiffness matrix
				int ndof = neln*4;
				ke.resize(ndof, ndof);
						
				// calculate pressure stiffness
				FluxStiffnessSS(el, ke, dt, m_bmixture);
						
				// assemble element matrix in global stiffness matrix
				psolver->AssembleStiffness(el.m_node, elm, ke);
			}
		}
	}
	else 
	{
		for (int m=0; m<N; ++m)
		{
			// get the surface element
			FESurfaceElement& el = m_psurf->Element(m);
			UnpackLM(el, elm);
					
			// calculate nodal normal fluid flux
			int neln = el.Nodes();
					
			if (!m_blinear || m_bmixture)
			{
				// get the element stiffness matrix
				int ndof = neln*4;
				ke.resize(ndof, ndof);
						
				// calculate pressure stiffness
				FluxStiffness(el, ke, dt, m_bmixture);
						
				// assemble element matrix in global stiffness matrix
				psolver->AssembleStiffness(el.m_node, elm, ke);
			}
		}
	}
}

//-----------------------------------------------------------------------------
void FEFluidFlux::Residual(const FETimeInfo& tp, FEGlobalVector& R)
{
	FEModel& fem = R.GetFEModel();
	double dt = tp.timeIncrement;

	vector<double> fe;

	vector<int> elm;

	int N = m_psurf->Elements();

	if (fem.GetCurrentStep()->m_nanalysis == FE_STEADY_STATE)
	{
		for (int i=0; i<N; ++i)
		{
			FESurfaceElement& el = m_psurf->Element(i);
			UnpackLM(el, elm);
				
			// calculate nodal normal fluid flux
			int neln = el.Nodes();
				
			int ndof = 4*neln;
			fe.resize(ndof);
				
			if (m_blinear == true) 
				LinearFlowRateSS(el, fe, dt, m_bmixture);
			else
				FlowRateSS(el, fe, dt, m_bmixture);
				
			// add element force vector to global force vector
			R.Assemble(el.m_node, elm, fe);
		}
	}
	else {
		for (int i=0; i<N; ++i)
		{
			FESurfaceElement& el = m_psurf->Element(i);
			UnpackLM(el, elm);
				
			// calculate nodal normal fluid flux
			int neln = el.Nodes();
				
			int ndof = 4*neln;
			fe.resize(ndof);
				
			if (m_blinear == true) 
				LinearFlowRate(el, fe, dt, m_bmixture);
			else
				FlowRate(el, fe, dt, m_bmixture);
				
			// add element force vector to global force vector
			R.Assemble(el.m_node, elm, fe);
		}
	}
}
