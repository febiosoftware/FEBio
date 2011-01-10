#include "stdafx.h"
#include "FEFluxSurface.h"
#include "FESolidSolver.h"

//-----------------------------------------------------------------------------
//! calculates the stiffness contribution due to fluid flux

void FEFluxSurface::FluxStiffness(FESurfaceElement& el, matrix& ke, vector<double>& wn, double dt, bool mixture)
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
	vec3d* rt = el.rt();
	vec3d* vt = el.vt();
	
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
//! calculates the equivalent nodal volumetric flow rates due to fluid flux

bool FEFluxSurface::FlowRate(FESurfaceElement& el, vector<double>& fe, vector<double>& wn, double dt, bool mixture)
{
	int i, n;

	// nr integration points
	int nint = el.GaussPoints();

	// nr of element nodes
	int neln = el.Nodes();

	// nodal coordinates and velocities
	vec3d* rt = el.rt();
	vec3d* vt = el.vt();
	
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

bool FEFluxSurface::LinearFlowRate(FESurfaceElement& el, vector<double>& fe, vector<double>& wn, double dt, bool mixture)
{
	int i, n;

	// nr integration points
	int nint = el.GaussPoints();

	// nr of element nodes
	int neln = el.Nodes();

	// nodal coordinates and velocity
	vec3d *r0 = el.r0();
	vec3d *rt = el.rt();
	vec3d *vt = el.vt();

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

void FEFluxSurface::Serialize(FEM& fem, DumpFile& ar)
{
	if (ar.IsSaving())
	{
		size_t i;
		for (i=0; i<m_el.size(); ++i)
		{
			FESurfaceElement& el = m_el[i];
			ar << el.Type();
			ar << el.GetMatID();
			ar << el.m_nID;
			ar << el.m_nrigid;
			ar << el.m_node;
			ar << el.m_lnode;
		}

		// fluid fluxes
		for (i=0; i<m_PC.size(); ++i)
		{
			FEFluidFlux& fc = m_PC[i];
			ar << fc.blinear << fc.mixture << fc.face << fc.lc;
			ar << fc.s[0] << fc.s[1] << fc.s[2] << fc.s[3];
			ar << fc.bc;
		}
	}
	else
	{
		size_t i; 
		int m, mat;
		for (i=0; i<m_el.size(); ++i)
		{
			FESurfaceElement& el = m_el[i];
			ar >> m;
			el.SetType(m);

			ar >> mat; el.SetMatID(mat);
			ar >> el.m_nID;
			ar >> el.m_nrigid;
			ar >> el.m_node;
			ar >> el.m_lnode;
		}

		// fluid fluxes
		for (i=0; i<m_PC.size(); ++i)
		{
			FEFluidFlux& fc = m_PC[i];
			ar >> fc.blinear  >> fc.mixture >> fc.face >> fc.lc;
			ar >> fc.s[0] >> fc.s[1] >> fc.s[2] >> fc.s[3];
			ar >> fc.bc;
		}

		// initialize surface data
		Init();
	}
}

void FEFluxSurface::StiffnessMatrix(FESolidSolver* psolver)
{
	FEM& fem = psolver->m_fem;
	double dt = fem.m_pStep->m_dt;

	matrix ke;

	int nfr = m_PC.size();
	for (int m=0; m<nfr; ++m)
	{
		FEFluidFlux& fc = m_PC[m];
		if ((fc.bc == 0) && fc.IsActive())
		{
			// get the surface element
			FESurfaceElement& el = m_el[m];

			// skip rigid surface elements
			// TODO: do we really need to skip rigid elements?
			if (!el.IsRigid())
			{
				UnpackElement(el);

				// calculate nodal normal fluid flux
				int neln = el.Nodes();
				vector<double> wn(neln);

				if (!fc.blinear)
				{
					double g = fem.GetLoadCurve(fc.lc)->Value();

					for (int j=0; j<neln; ++j) wn[j] = g*fc.s[j];
					
					// get the element stiffness matrix
					int ndof = neln*4;
					ke.Create(ndof, ndof);

					// calculate pressure stiffness
					FluxStiffness(el, ke, wn, dt, fc.mixture);

					// assemble element matrix in global stiffness matrix
					psolver->AssembleStiffness(el.m_node, el.LM(), ke);
				}
			}
		}
	}
}

void FEFluxSurface::Residual(FESolidSolver* psolver, vector<double>& R)
{
	FEM& fem = psolver->m_fem;
	double dt = fem.m_pStep->m_dt;

	vector<double> fe;

	int nfr = m_PC.size();
	for (int i=0; i<nfr; ++i)
	{
		FEFluidFlux& fc = m_PC[i];
		if ((fc.bc == 0) && fc.IsActive())
		{
			FESurfaceElement& el = m_el[i];
			UnpackElement(el);

			// calculate nodal normal fluid flux
			int neln = el.Nodes();
			vector<double> wn(neln);

			double g = fem.GetLoadCurve(fc.lc)->Value();

			for (int j=0; j<neln; ++j) wn[j] = g*fc.s[j];

			int ndof = 4*neln;
			fe.resize(ndof);

			if (fc.blinear) 
				LinearFlowRate(el, fe, wn, dt, fc.mixture);
			else
				FlowRate(el, fe, wn, dt, fc.mixture);

			// add element force vector to global force vector
			psolver->AssembleResidual(el.m_node, el.LM(), fe, R);
		}
	}
}

