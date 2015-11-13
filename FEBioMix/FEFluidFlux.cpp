#include "FEFluidFlux.h"
#include "FECore/FESolver.h"
#include "FECore/FEModel.h"

//-----------------------------------------------------------------------------
FEFluidFlux::LOAD::LOAD()
{ 
	s[0] = s[1] = s[2] = s[3] = s[4] = s[5] = s[6] = s[7] = s[8] = 1.0; 
	lc = -1; 
}

//-----------------------------------------------------------------------------
BEGIN_PARAMETER_LIST(FEFluidFlux, FESurfaceLoad)
	ADD_PARAMETER(m_flux    , FE_PARAM_DOUBLE, "flux"   );
	ADD_PARAMETER(m_blinear , FE_PARAM_BOOL  , "linear" );
	ADD_PARAMETER(m_bmixture, FE_PARAM_BOOL  , "mixture");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FEFluidFlux::FEFluidFlux(FEModel* pfem) : FESurfaceLoad(pfem)
{ 
	m_blinear = false; 
	m_bmixture = false; 
	m_flux = 1.0;
}

//-----------------------------------------------------------------------------
void FEFluidFlux::Create(int n) 
{ 
	m_PC.resize(n); 
}

//-----------------------------------------------------------------------------
//! \deprecated This is only used in the 1.2 file format which is obsolete
bool FEFluidFlux::SetAttribute(const char* szatt, const char* szval)
{
	if (strcmp(szatt, "type") == 0)
	{
		if      (strcmp(szval, "linear"   ) == 0) SetLinear(true );
		else if (strcmp(szval, "nonlinear") == 0) SetLinear(false);
		else return false;
		return true;
	}
	else if (strcmp(szatt, "flux") == 0)
	{
		if      (strcmp(szval, "mixture") == 0) SetMixture(true);
		else if (strcmp(szval, "fluid"  ) == 0) SetMixture(false); 
		else return false;
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
bool FEFluidFlux::SetFacetAttribute(int nface, const char* szatt, const char* szval)
{
	LOAD& pc = FluidFlux(nface);
	if      (strcmp(szatt, "id") == 0) {}
	else if (strcmp(szatt, "lc") == 0) pc.lc = atoi(szval) - 1;
	else if (strcmp(szatt, "scale") == 0)
	{
		double s = atof(szval);
		pc.s[0] = pc.s[1] = pc.s[2] = pc.s[3] = s;
		pc.s[4] = pc.s[5] = pc.s[6] = pc.s[7] = s;
		pc.s[8] = s;
	}
	else return false;

	return true;
}

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
	vec3d rt[FEElement::MAX_NODES], vt[FEElement::MAX_NODES];
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
	vec3d rt[FEElement::MAX_NODES];
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
	vec3d rt[FEElement::MAX_NODES], vt[FEElement::MAX_NODES];
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
	vec3d rt[FEElement::MAX_NODES];
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
	vec3d r0[FEElement::MAX_NODES], rt[FEElement::MAX_NODES], vt[FEElement::MAX_NODES];
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
	vec3d r0[FEElement::MAX_NODES], rt[FEElement::MAX_NODES];
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
		ar << m_flux;
		ar << m_blinear << m_bmixture;
		ar << (int) m_PC.size();
		for (int i=0; i<(int) m_PC.size(); ++i)
		{
			LOAD& fc = m_PC[i];
			ar << fc.lc;
			ar << fc.s[0] << fc.s[1] << fc.s[2] << fc.s[3];
			ar << fc.s[4] << fc.s[5] << fc.s[6] << fc.s[7] << fc.s[8];
		}
	}
	else
	{
		int n;
		ar >> m_flux;
		ar >> m_blinear >> m_bmixture;
		ar >> n;
		m_PC.resize(n);
		for (int i=0; i<n; ++i)
		{
			LOAD& fc = m_PC[i];
			ar >> fc.lc;
			ar >> fc.s[0] >> fc.s[1] >> fc.s[2] >> fc.s[3];
			ar >> fc.s[4] >> fc.s[5] >> fc.s[6] >> fc.s[7] >> fc.s[8];
		}
	}
}

//-----------------------------------------------------------------------------
void FEFluidFlux::StiffnessMatrix(FESolver* psolver)
{
	FEModel& fem = psolver->GetFEModel();
	double dt = fem.GetCurrentStep()->m_dt;

	matrix ke;

	vector<int> elm;

	int nfr = m_PC.size();

	if (fem.GetCurrentStep()->m_nanalysis == FE_STEADY_STATE)
	{
		for (int m=0; m<nfr; ++m)
		{
			LOAD& fc = m_PC[m];

			// get the surface element
			FESurfaceElement& el = m_psurf->Element(m);
			m_psurf->UnpackLM(el, elm);
				
			// calculate nodal normal fluid flux
			int neln = el.Nodes();
			vector<double> wn(neln);
					
			if (!m_blinear || m_bmixture)
			{
				double g = m_flux;
				if (fc.lc >= 0) g *= fem.GetLoadCurve(fc.lc)->Value();
						
				for (int j=0; j<neln; ++j) wn[j] = g*fc.s[j];
						
				// get the element stiffness matrix
				int ndof = neln*4;
				ke.resize(ndof, ndof);
						
				// calculate pressure stiffness
				FluxStiffnessSS(el, ke, wn, dt, m_bmixture);
						
				// assemble element matrix in global stiffness matrix
				psolver->AssembleStiffness(el.m_node, elm, ke);
			}
		}
	}
	else 
	{
		for (int m=0; m<nfr; ++m)
		{
			LOAD& fc = m_PC[m];

			// get the surface element
			FESurfaceElement& el = m_psurf->Element(m);
			m_psurf->UnpackLM(el, elm);
					
			// calculate nodal normal fluid flux
			int neln = el.Nodes();
			vector<double> wn(neln);
					
			if (!m_blinear || m_bmixture)
			{
				double g = m_flux;
				if (fc.lc >= 0) g *= fem.GetLoadCurve(fc.lc)->Value();
						
				for (int j=0; j<neln; ++j) wn[j] = g*fc.s[j];
						
				// get the element stiffness matrix
				int ndof = neln*4;
				ke.resize(ndof, ndof);
						
				// calculate pressure stiffness
				FluxStiffness(el, ke, wn, dt, m_bmixture);
						
				// assemble element matrix in global stiffness matrix
				psolver->AssembleStiffness(el.m_node, elm, ke);
			}
		}
	}
}

//-----------------------------------------------------------------------------
void FEFluidFlux::Residual(FEGlobalVector& R)
{
	FEModel& fem = R.GetFEModel();
	double dt = fem.GetCurrentStep()->m_dt;

	vector<double> fe;

	vector<int> elm;

	int nfr = m_PC.size();

	if (fem.GetCurrentStep()->m_nanalysis == FE_STEADY_STATE)
	{
		for (int i=0; i<nfr; ++i)
		{
			LOAD& fc = m_PC[i];

			FESurfaceElement& el = m_psurf->Element(i);
			m_psurf->UnpackLM(el, elm);
				
			// calculate nodal normal fluid flux
			int neln = el.Nodes();
			vector<double> wn(neln);
				
			double g = m_flux;
			if (fc.lc >= 0) g *= fem.GetLoadCurve(fc.lc)->Value();
				
			for (int j=0; j<neln; ++j) wn[j] = g*fc.s[j];
				
			int ndof = 4*neln;
			fe.resize(ndof);
				
			if (m_blinear == true) 
				LinearFlowRateSS(el, fe, wn, dt, m_bmixture);
			else
				FlowRateSS(el, fe, wn, dt, m_bmixture);
				
			// add element force vector to global force vector
			R.Assemble(el.m_node, elm, fe);
		}
	}
	else {
		for (int i=0; i<nfr; ++i)
		{
			LOAD& fc = m_PC[i];

			FESurfaceElement& el = m_psurf->Element(i);
			m_psurf->UnpackLM(el, elm);
				
			// calculate nodal normal fluid flux
			int neln = el.Nodes();
			vector<double> wn(neln);
				
			double g = m_flux;
			if (fc.lc >= 0) g *= fem.GetLoadCurve(fc.lc)->Value();
				
			for (int j=0; j<neln; ++j) wn[j] = g*fc.s[j];
				
			int ndof = 4*neln;
			fe.resize(ndof);
				
			if (m_blinear == true) 
				LinearFlowRate(el, fe, wn, dt, m_bmixture);
			else
				FlowRate(el, fe, wn, dt, m_bmixture);
				
			// add element force vector to global force vector
			R.Assemble(el.m_node, elm, fe);
		}
	}
}
