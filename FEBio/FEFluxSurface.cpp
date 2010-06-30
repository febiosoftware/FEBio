#include "stdafx.h"
#include "FEFluxSurface.h"
#include "FESolidSolver.h"

//-----------------------------------------------------------------------------
//! calculates the stiffness contribution due to fluid flux

void FEFluxSurface::FluxStiffness(FESurfaceElement& el, matrix& ke)
{
	int i, j, n;

	int nint = el.GaussPoints();
	int neln = el.Nodes();

	// fluid flux at nodes
	double *wn = el.pt();

	// fluid flux at integration point
	double wr;

	// gauss weights
	double* w = el.GaussWeights();

	// nodal coordinates
	vec3d* rt = el.rt();

	vec3d kab;

	ke.zero();

	double* N, *Gr, *Gs;
	
	// repeat over integration points
	for (n=0; n<nint; ++n)
	{
		N = el.H(n);
		Gr = el.Gr(n);
		Gs = el.Gs(n);

		// calculate fluid flux at integration point
		wr = 0;
		for (i=0; i<neln; ++i) wr += N[i]*wn[i];
		
		// calculate covariant basis vectors (surface tangents)
		vec3d gr, gs;
		for (i=0; i<neln; ++i)
		{
			gr += rt[i]*Gr[i];
			gs += rt[i]*Gs[i];
		}
		
		// calculate surface normal
		vec3d gn = gr ^ gs;
		double J = gn.unit();

		// calculate stiffness component
		for (i=0; i<neln; ++i)
			for (j=0; j<neln; ++j)
			{
				kab = (gn ^ (gs*Gr[j] - gr*Gs[j]))*(N[i]*wr*w[n]);
				ke[4*i+3][4*j  ] = kab.x;
				ke[4*i+3][4*j+1] = kab.y;
				ke[4*i+3][4*j+2] = kab.z;
			}
	}
}

//-----------------------------------------------------------------------------
//! calculates the equivalent nodal volumetric flow rates due to fluid flux

bool FEFluxSurface::FlowRate(FESurfaceElement& el, vector<double>& fe)
{
	int i, n;

	// nr integration points
	int nint = el.GaussPoints();

	// nr of element nodes
	int neln = el.Nodes();

	// fluid flux at nodes
	double *wn = el.pt();

	// nodal coordinates
	vec3d *rt = el.rt();

	double* Gr, *Gs;
	double* N;
	double* w  = el.GaussWeights();

	// fluid flux at integration points
	double wr;

	vec3d dxr, dxs;

	// volumetric flow rate
	double f;

	// repeat over integration points
	fe.zero();
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

			dxr.x += Gr[i]*rt[i].x;
			dxr.y += Gr[i]*rt[i].y;
			dxr.z += Gr[i]*rt[i].z;

			dxs.x += Gs[i]*rt[i].x;
			dxs.y += Gs[i]*rt[i].y;
			dxs.z += Gs[i]*rt[i].z;
		}

		f = (dxr ^ dxs).norm()*wr*w[n];

		for (i=0; i<neln; ++i)
		{
			fe[i] += N[i]*f;
		}
	}

	return true;
}

//-----------------------------------------------------------------------------
//! calculates the equivalent nodal volumetric flow rates due to fluid flux

bool FEFluxSurface::LinearFlowRate(FESurfaceElement& el, vector<double>& fe)
{
	int i, n;

	// nr integration points
	int nint = el.GaussPoints();

	// nr of element nodes
	int neln = el.Nodes();

	// fluid flux at nodes
	double *wn = el.pt();

	// nodal coordinates
	vec3d *r0 = el.r0();

	double* Gr, *Gs;
	double* N;
	double* w  = el.GaussWeights();

	// flux at integration points
	double wr;

	vec3d dxr, dxs;

	// volumetric flow rate
	double f;

	// repeat over integration points
	fe.zero();
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

			dxr.x += Gr[i]*r0[i].x;
			dxr.y += Gr[i]*r0[i].y;
			dxr.z += Gr[i]*r0[i].z;

			dxs.x += Gs[i]*r0[i].x;
			dxs.y += Gs[i]*r0[i].y;
			dxs.z += Gs[i]*r0[i].z;
		}

		f = (dxr ^ dxs).norm()*wr*w[n];

		for (i=0; i<neln; ++i)
		{
			fe[i] += N[i]*f;
		}
	}

	return true;
}

//-----------------------------------------------------------------------------

void FEFluxSurface::Serialize(FEM& fem, Archive& ar)
{
	if (ar.IsSaving())
	{
		int i;
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
			ar << fc.blinear << fc.face << fc.lc;
			ar << fc.s[0] << fc.s[1] << fc.s[2] << fc.s[3];
			ar << fc.bc;
		}
	}
	else
	{
		int i, m, mat;
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
			ar >> fc.blinear >> fc.face >> fc.lc;
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

				// calculate nodal pressures
				double* pt = el.pt();

				if (!fc.blinear)
				{
					double g = fem.GetLoadCurve(fc.lc)->Value();

					for (int j=0; j<el.Nodes(); ++j) pt[j] = g*fc.s[j]*dt;

					// get the element stiffness matrix
					int neln = el.Nodes();
					int ndof = neln*4;
					ke.Create(ndof, ndof);

					// calculate pressure stiffness
					FluxStiffness(el, ke);

					// TODO: the problem here is that the LM array that is returned by the UnpackElement
					// function does not give the equation numbers in the right order. For this reason we
					// have to create a new lm array and place the equation numbers in the right order.
					// What we really ought to do is fix the UnpackElement function so that it returns
					// the LM vector in the right order for poroelastic elements.
					vector<int> lm(ndof);
					for (int i=0; i<neln; ++i)
					{
						lm[4*i  ] = el.LM()[3*i];
						lm[4*i+1] = el.LM()[3*i+1];
						lm[4*i+2] = el.LM()[3*i+2];
						lm[4*i+3] = el.LM()[3*neln+i];
					}
					
					// assemble element matrix in global stiffness matrix
					psolver->AssembleStiffness(el.m_node, lm, ke);
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

			// calculate nodal fluid fluxes
			double* pt = el.pt();

			double g = fem.GetLoadCurve(fc.lc)->Value();

			for (int j=0; j<el.Nodes(); ++j) pt[j] = g*fc.s[j]*dt;

			int ndof = el.Nodes();
			fe.resize(ndof);

			if (fc.blinear) LinearFlowRate(el, fe); else FlowRate(el, fe);

			// add element force vector to global force vector
			int *lm = el.LM() ,J;
			for (int k=0; k<ndof; ++k)
			{
				J = lm[3*ndof+k];
				if (J >= 0) R[J] += fe[k];
			}
		}
	}
}

