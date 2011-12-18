#include "stdafx.h"
#include "FEFluidFlux.h"
#include "FESolidSolver.h"

//-----------------------------------------------------------------------------
//! calculates the stiffness contribution due to fluid flux

void FEFluidFlux::FluxStiffness(FESurfaceElement& el, matrix& ke, vector<double>& wn, double dt, bool mixture)
{
	int i, j, n;

	int nint = el.GaussPoints();
	int neln = el.Nodes();

	// normal fluid flux at integration point
	double wr;
	
	// solid velocity at integration point
	vec3d vr;

	vec3d dxr, dxs, dxt;

	// gauss weights
	double* w = el.GaussWeights();

	// nodal coordinates and velocities
	vec3d rt[4], vt[4];
	for (i=0; i<neln; ++i)
	{
		rt[i] = m_psurf->GetMesh()->Node(el.m_node[i]).m_rt;
		vt[i] = m_psurf->GetMesh()->Node(el.m_node[i]).m_vt;
	}
	
	vec3d kab, t1, t2;

	ke.zero();

	double* N, *Gr, *Gs;
	
	// repeat over integration points
	for (n=0; n<nint; ++n)
	{
		N = el.H(n);
		Gr = el.Gr(n);
		Gs = el.Gs(n);

		// calculate velocities and covariant basis vectors at integration point
		wr = 0;
		vr = dxr = dxs = vec3d(0,0,0);
		for (i=0; i<neln; ++i)
		{
			wr += N[i]*wn[i];
			vr += vt[i]*N[i];
			dxr += rt[i]*Gr[i];
			dxs += rt[i]*Gs[i];
		}
		
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

void FEFluidFlux::FluxStiffnessSS(FESurfaceElement& el, matrix& ke, vector<double>& wn, double dt, bool mixture)
{
	int i, j, n;
	
	int nint = el.GaussPoints();
	int neln = el.Nodes();
	
	// normal fluid flux at integration point
	double wr;
	
	vec3d dxr, dxs, dxt;
	
	// gauss weights
	double* w = el.GaussWeights();
	
	// nodal coordinates and velocities
	vec3d rt[4];
	for (i=0; i<neln; ++i)
	{
		rt[i] = m_psurf->GetMesh()->Node(el.m_node[i]).m_rt;
	}
	
	vec3d kab, t1, t2;
	
	ke.zero();
	
	double* N, *Gr, *Gs;
	
	// repeat over integration points
	for (n=0; n<nint; ++n)
	{
		N = el.H(n);
		Gr = el.Gr(n);
		Gs = el.Gs(n);
		
		// calculate velocities and covariant basis vectors at integration point
		wr = 0;
		dxr = dxs = vec3d(0,0,0);
		for (i=0; i<neln; ++i)
		{
			wr += N[i]*wn[i];
			dxr += rt[i]*Gr[i];
			dxs += rt[i]*Gs[i];
		}
		
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

bool FEFluidFlux::FlowRate(FESurfaceElement& el, vector<double>& fe, vector<double>& wn, double dt, bool mixture)
{
	int i, n;

	// nr integration points
	int nint = el.GaussPoints();

	// nr of element nodes
	int neln = el.Nodes();

	// nodal coordinates and velocities
	vec3d rt[4], vt[4];
	for (i=0; i<neln; ++i)
	{
		rt[i] = m_psurf->GetMesh()->Node(el.m_node[i]).m_rt;
		vt[i] = m_psurf->GetMesh()->Node(el.m_node[i]).m_vt;
	}
	
	double* Gr, *Gs;
	double* N;
	double* w  = el.GaussWeights();

	// normal fluid flux at integration points
	double wr;

	// solid velocity at integration point
	vec3d vr;
	
	vec3d dxr, dxs, dxt;

	// volumetric flow rate
	double f;

	// repeat over integration points
	zero(fe);
	for (n=0; n<nint; ++n)
	{
		N  = el.H(n);
		Gr = el.Gr(n);
		Gs = el.Gs(n);

		wr = 0;
		vr = dxr = dxs = vec3d(0,0,0);
		for (i=0; i<neln; ++i) 
		{
			wr += N[i]*wn[i];
			vr += vt[i]*N[i];
			dxr += rt[i]*Gr[i];
			dxs += rt[i]*Gs[i];
		}
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

bool FEFluidFlux::FlowRateSS(FESurfaceElement& el, vector<double>& fe, vector<double>& wn, double dt, bool mixture)
{
	int i, n;
	
	// nr integration points
	int nint = el.GaussPoints();
	
	// nr of element nodes
	int neln = el.Nodes();
	
	// nodal coordinates and velocities
	vec3d rt[4];
	for (i=0; i<neln; ++i)
	{
		rt[i] = m_psurf->GetMesh()->Node(el.m_node[i]).m_rt;
	}
	
	double* Gr, *Gs;
	double* N;
	double* w  = el.GaussWeights();
	
	// normal fluid flux at integration points
	double wr;
	
	vec3d dxr, dxs, dxt;
	
	// volumetric flow rate
	double f;
	
	// repeat over integration points
	zero(fe);
	for (n=0; n<nint; ++n)
	{
		N  = el.H(n);
		Gr = el.Gr(n);
		Gs = el.Gs(n);
		
		wr = 0;
		dxr = dxs = vec3d(0,0,0);
		for (i=0; i<neln; ++i) 
		{
			wr += N[i]*wn[i];
			dxr += rt[i]*Gr[i];
			dxs += rt[i]*Gs[i];
		}
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

bool FEFluidFlux::LinearFlowRate(FESurfaceElement& el, vector<double>& fe, vector<double>& wn, double dt, bool mixture)
{
	int i, n;

	// nr integration points
	int nint = el.GaussPoints();

	// nr of element nodes
	int neln = el.Nodes();
	assert(neln <= 4);

	FEMesh& mesh = *m_psurf->GetMesh();

	// nodal coordinates and velocity
	vec3d r0[4], rt[4], vt[4];
	for (i=0; i<neln; ++i)
	{
		r0[i] = mesh.Node(el.m_node[i]).m_r0;
		rt[i] = mesh.Node(el.m_node[i]).m_rt;
		vt[i] = mesh.Node(el.m_node[i]).m_vt;
	}

	double* Gr, *Gs;
	double* N;
	double* w  = el.GaussWeights();

	// normal fluid flux at integration points
	double Wr;

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
		N  = el.H(n);
		Gr = el.Gr(n);
		Gs = el.Gs(n);

		Wr = 0;
		dXr = dXs = vec3d(0,0,0);
		vr = dxr = dxs = vec3d(0,0,0);
		for (i=0; i<neln; ++i) 
		{
			Wr += N[i]*wn[i];
			vr += vt[i]*N[i];
			dXr += r0[i]*Gr[i];
			dXs += r0[i]*Gs[i];
			dxr += rt[i]*Gr[i];
			dxs += rt[i]*Gs[i];
		}
		dXt = dXr ^ dXs;
		dxt = dxr ^ dxs;
		
		f = (dXt.norm()*Wr - (vr*dxt)*mixture)*w[n]*dt;

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

bool FEFluidFlux::LinearFlowRateSS(FESurfaceElement& el, vector<double>& fe, vector<double>& wn, double dt, bool mixture)
{
	int i, n;
	
	// nr integration points
	int nint = el.GaussPoints();
	
	// nr of element nodes
	int neln = el.Nodes();
	assert(neln <= 4);
	
	FEMesh& mesh = *m_psurf->GetMesh();
	
	// nodal coordinates and velocity
	vec3d r0[4], rt[4];
	for (i=0; i<neln; ++i)
	{
		r0[i] = mesh.Node(el.m_node[i]).m_r0;
		rt[i] = mesh.Node(el.m_node[i]).m_rt;
	}
	
	double* Gr, *Gs;
	double* N;
	double* w  = el.GaussWeights();
	
	// normal fluid flux at integration points
	double Wr;
	
	vec3d dXr, dXs, dXt;
	vec3d dxr, dxs, dxt;
	
	// volumetric flow rate
	double f;
	
	// repeat over integration points
	zero(fe);
	for (n=0; n<nint; ++n)
	{
		N  = el.H(n);
		Gr = el.Gr(n);
		Gs = el.Gs(n);
		
		Wr = 0;
		dXr = dXs = vec3d(0,0,0);
		dxr = dxs = vec3d(0,0,0);
		for (i=0; i<neln; ++i) 
		{
			Wr += N[i]*wn[i];
			dXr += r0[i]*Gr[i];
			dXs += r0[i]*Gs[i];
			dxr += rt[i]*Gr[i];
			dxs += rt[i]*Gs[i];
		}
		dXt = dXr ^ dXs;
		dxt = dxr ^ dxs;
		
		f = (dXt.norm()*Wr)*w[n]*dt;
		
		for (i=0; i<neln; ++i)
		{
			fe[3*neln+i] += N[i]*f;
		}
	}
	
	return true;
}

//-----------------------------------------------------------------------------
void FEFluidFlux::Serialize(DumpFile& ar)
{
	if (ar.IsSaving())
	{
		ar << m_blinear << m_bmixture;
		ar << (int) m_PC.size();
		for (int i=0; i<(int) m_PC.size(); ++i)
		{
			LOAD& fc = m_PC[i];
			ar << fc.lc;
			ar << fc.s[0] << fc.s[1] << fc.s[2] << fc.s[3];
			ar << fc.bc;
		}
	}
	else
	{
		int n;
		ar >> m_blinear >> m_bmixture;
		ar >> n;
		m_PC.resize(n);
		for (int i=0; i<n; ++i)
		{
			LOAD& fc = m_PC[i];
			ar >> fc.lc;
			ar >> fc.s[0] >> fc.s[1] >> fc.s[2] >> fc.s[3];
			ar >> fc.bc;
		}
	}
}

//-----------------------------------------------------------------------------
void FEFluidFlux::StiffnessMatrix(FESolver* psolver)
{
	FESolidSolver& solver = dynamic_cast<FESolidSolver&>(*psolver);

	FEM& fem = solver.m_fem;
	double dt = fem.m_pStep->m_dt;

	matrix ke;

	vector<int> elm;

	int nfr = m_PC.size();

	if (fem.m_pStep->m_nanalysis == FE_STEADY_STATE)
	{
		for (int m=0; m<nfr; ++m)
		{
			LOAD& fc = m_PC[m];
			if (fc.bc == 0)
			{
				// get the surface element
				FESurfaceElement& el = m_psurf->Element(m);
				
				// skip rigid surface elements
				// TODO: do we really need to skip rigid elements?
				if (!el.IsRigid())
				{
					m_psurf->UnpackLM(el, elm);
					
					// calculate nodal normal fluid flux
					int neln = el.Nodes();
					vector<double> wn(neln);
					
					if (!m_blinear || m_bmixture)
					{
						double g = fem.GetLoadCurve(fc.lc)->Value();
						
						for (int j=0; j<neln; ++j) wn[j] = g*fc.s[j];
						
						// get the element stiffness matrix
						int ndof = neln*4;
						ke.resize(ndof, ndof);
						
						// calculate pressure stiffness
						FluxStiffnessSS(el, ke, wn, dt, m_bmixture);
						
						// assemble element matrix in global stiffness matrix
						solver.AssembleStiffness(el.m_node, elm, ke);
					}
				}
			}
		}
	}
	else 
	{
		for (int m=0; m<nfr; ++m)
		{
			LOAD& fc = m_PC[m];
			if (fc.bc == 0)
			{
				// get the surface element
				FESurfaceElement& el = m_psurf->Element(m);
				
				// skip rigid surface elements
				// TODO: do we really need to skip rigid elements?
				if (!el.IsRigid())
				{
					m_psurf->UnpackLM(el, elm);
					
					// calculate nodal normal fluid flux
					int neln = el.Nodes();
					vector<double> wn(neln);
					
					if (!m_blinear || m_bmixture)
					{
						double g = fem.GetLoadCurve(fc.lc)->Value();
						
						for (int j=0; j<neln; ++j) wn[j] = g*fc.s[j];
						
						// get the element stiffness matrix
						int ndof = neln*4;
						ke.resize(ndof, ndof);
						
						// calculate pressure stiffness
						FluxStiffness(el, ke, wn, dt, m_bmixture);
						
						// assemble element matrix in global stiffness matrix
						solver.AssembleStiffness(el.m_node, elm, ke);
					}
				}
			}
		}
	}

}

void FEFluidFlux::Residual(FESolver* psolver, vector<double>& R)
{
	FESolidSolver& solver = dynamic_cast<FESolidSolver&>(*psolver);

	FEM& fem = solver.m_fem;
	double dt = fem.m_pStep->m_dt;

	vector<double> fe;

	vector<int> elm;

	int nfr = m_PC.size();

	if (fem.m_pStep->m_nanalysis == FE_STEADY_STATE)
	{
		for (int i=0; i<nfr; ++i)
		{
			LOAD& fc = m_PC[i];
			if (fc.bc == 0)
			{
				FESurfaceElement& el = m_psurf->Element(i);
				m_psurf->UnpackLM(el, elm);
				
				// calculate nodal normal fluid flux
				int neln = el.Nodes();
				vector<double> wn(neln);
				
				double g = fem.GetLoadCurve(fc.lc)->Value();
				
				for (int j=0; j<neln; ++j) wn[j] = g*fc.s[j];
				
				int ndof = 4*neln;
				fe.resize(ndof);
				
				if (m_blinear == true) 
					LinearFlowRateSS(el, fe, wn, dt, m_bmixture);
				else
					FlowRateSS(el, fe, wn, dt, m_bmixture);
				
				// add element force vector to global force vector
				solver.AssembleResidual(el.m_node, elm, fe, R);
			}
		}
	}
	else {
		for (int i=0; i<nfr; ++i)
		{
			LOAD& fc = m_PC[i];
			if (fc.bc == 0)
			{
				FESurfaceElement& el = m_psurf->Element(i);
				m_psurf->UnpackLM(el, elm);
				
				// calculate nodal normal fluid flux
				int neln = el.Nodes();
				vector<double> wn(neln);
				
				double g = fem.GetLoadCurve(fc.lc)->Value();
				
				for (int j=0; j<neln; ++j) wn[j] = g*fc.s[j];
				
				int ndof = 4*neln;
				fe.resize(ndof);
				
				if (m_blinear == true) 
					LinearFlowRate(el, fe, wn, dt, m_bmixture);
				else
					FlowRate(el, fe, wn, dt, m_bmixture);
				
				// add element force vector to global force vector
				solver.AssembleResidual(el.m_node, elm, fe, R);
			}
		}
	}

}
