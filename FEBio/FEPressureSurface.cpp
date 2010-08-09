#include "stdafx.h"
#include "FEPressureSurface.h"
#include "FESolidSolver.h"

//-----------------------------------------------------------------------------
//! calculates the stiffness contribution due to hydrostatic pressure

void FEPressureSurface::PressureStiffness(FESurfaceElement& el, matrix& ke, vector<double>& tn)
{
	int i, j, n;

	int nint = el.GaussPoints();
	int neln = el.Nodes();

	// traction at integration point
	double tr;
	
	vec3d dxr, dxs;

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

		tr = 0;
		dxr = dxs = vec3d(0,0,0);
		for (i=0; i<neln; ++i) 
		{
			tr += N[i]*tn[i];
			dxr += rt[i]*Gr[i];
			dxs += rt[i]*Gs[i];
		}
		
		// calculate stiffness component
		for (i=0; i<neln; ++i)
			for (j=0; j<neln; ++j)
			{
				kab = (dxr*(N[j]*Gs[i]-N[i]*Gs[j])
					   -dxs*(N[j]*Gr[i]-N[i]*Gr[j]))*w[n]*0.5*tr;

				ke[3*i  ][3*j  ] +=      0;
				ke[3*i  ][3*j+1] += -kab.z;
				ke[3*i  ][3*j+2] +=  kab.y;

				ke[3*i+1][3*j  ] +=  kab.z;
				ke[3*i+1][3*j+1] +=      0;
				ke[3*i+1][3*j+2] += -kab.x;

				ke[3*i+2][3*j  ] += -kab.y;
				ke[3*i+2][3*j+1] +=  kab.x;
				ke[3*i+2][3*j+2] +=      0;
			}
	}
}

//-----------------------------------------------------------------------------
//! calculates the equivalent nodal forces due to hydrostatic pressure

bool FEPressureSurface::PressureForce(FESurfaceElement& el, vector<double>& fe, vector<double>& tn)
{
	int i, n;

	// nr integration points
	int nint = el.GaussPoints();

	// nr of element nodes
	int neln = el.Nodes();

	// nodal coordinates
	vec3d *rt = el.rt();

	double* Gr, *Gs;
	double* N;
	double* w  = el.GaussWeights();

	// traction at integration points
	double tr;

	vec3d dxr, dxs;

	// force vector
	vec3d f;

	// repeat over integration points
	zero(fe);
	for (n=0; n<nint; ++n)
	{
		N  = el.H(n);
		Gr = el.Gr(n);
		Gs = el.Gs(n);

		tr = 0;
		dxr = dxs = vec3d(0,0,0);
		for (i=0; i<neln; ++i) 
		{
			tr += N[i]*tn[i];
			dxr += rt[i]*Gr[i];
			dxs += rt[i]*Gs[i];
		}

		f = (dxr ^ dxs)*tr*w[n];

		for (i=0; i<neln; ++i)
		{
			fe[3*i  ] += N[i]*f.x;
			fe[3*i+1] += N[i]*f.y;
			fe[3*i+2] += N[i]*f.z;
		}
	}

	return true;
}

//-----------------------------------------------------------------------------
//! calculates the equivalent nodal forces due to hydrostatic pressure

bool FEPressureSurface::LinearPressureForce(FESurfaceElement& el, vector<double>& fe, vector<double>& tn)
{
	int i, n;

	// nr integration points
	int nint = el.GaussPoints();

	// nr of element nodes
	int neln = el.Nodes();

	// nodal coordinates
	vec3d *r0 = el.r0();

	double* Gr, *Gs;
	double* N;
	double* w  = el.GaussWeights();

	// traction at integration points
	double tr;

	vec3d dxr, dxs;

	// force vector
	vec3d f;

	// repeat over integration points
	zero(fe);
	for (n=0; n<nint; ++n)
	{
		N  = el.H(n);
		Gr = el.Gr(n);
		Gs = el.Gs(n);

		tr = 0;
		dxr = dxs = vec3d(0,0,0);
		for (i=0; i<neln; ++i) 
		{
			tr += N[i]*tn[i];
			dxr += r0[i]*Gr[i];
			dxs += r0[i]*Gs[i];
		}

		f = (dxr ^ dxs)*tr*w[n];

		for (i=0; i<neln; ++i)
		{
			fe[3*i  ] += N[i]*f.x;
			fe[3*i+1] += N[i]*f.y;
			fe[3*i+2] += N[i]*f.z;
		}
	}

	return true;
}

//-----------------------------------------------------------------------------

void FEPressureSurface::Serialize(FEM& fem, Archive& ar)
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

		// pressure forces
		for (i=0; i<m_PC.size(); ++i)
		{
			FEPressureLoad& pc = m_PC[i];
			ar << pc.blinear << pc.face << pc.lc;
			ar << pc.s[0] << pc.s[1] << pc.s[2] << pc.s[3];
		}
	}
	else
	{
		int i, m, mat;
		for (i=0; i<(int) m_el.size(); ++i)
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

		// pressure forces
		for (i=0; i<(int) m_PC.size(); ++i)
		{
			FEPressureLoad& pc = m_PC[i];
			ar >> pc.blinear >> pc.face >> pc.lc;
			ar >> pc.s[0] >> pc.s[1] >> pc.s[2] >> pc.s[3];
		}

		// initialize surface data
		Init();
	}
}

void FEPressureSurface::StiffnessMatrix(FESolidSolver* psolver)
{
	FEM& fem = psolver->m_fem;

	matrix ke;

	int npr = m_PC.size();
	for (int m=0; m<npr; ++m)
	{
		FEPressureLoad& pc = m_PC[m];
		if (pc.IsActive())
		{
			// get the surface element
			FESurfaceElement& el = m_el[m];

			// skip rigid surface elements
			// TODO: do we really need to skip rigid elements?
			if (!el.IsRigid())
			{
				UnpackElement(el);

				// calculate nodal normal tractions
				int neln = el.Nodes();
				vector<double> tn(neln);

				if (!pc.blinear)
				{
					double g = fem.GetLoadCurve(pc.lc)->Value();

					// evaluate the prescribed traction.
					// note the negative sign. This is because this boundary condition uses the 
					// convention that a positive pressure is compressive
					for (int j=0; j<neln; ++j) tn[j] = -g*pc.s[j];

					// get the element stiffness matrix
					int ndof = 3*neln;
					ke.Create(ndof, ndof);

					// calculate pressure stiffness
					PressureStiffness(el, ke, tn);

					// assemble element matrix in global stiffness matrix
					psolver->AssembleStiffness(el.m_node, el.LM(), ke);
				}
			}
		}
	}
}

void FEPressureSurface::Residual(FESolidSolver* psolver, vector<double>& R)
{
	FEM& fem = psolver->m_fem;

	vector<double> fe;

	int npr = m_PC.size();
	for (int i=0; i<npr; ++i)
	{
		FEPressureLoad& pc = m_PC[i];
		if (pc.IsActive())
		{
			FESurfaceElement& el = m_el[i];
			UnpackElement(el);

			// calculate nodal normal tractions
			int neln = el.Nodes();
			vector<double> tn(neln);

			double g = fem.GetLoadCurve(pc.lc)->Value();

			// evaluate the prescribed traction.
			// note the negative sign. This is because this boundary condition uses the 
			// convention that a positive pressure is compressive
			for (int j=0; j<el.Nodes(); ++j) tn[j] = -g*pc.s[j];
			
			int ndof = 3*neln;
			fe.resize(ndof);

			if (pc.blinear) LinearPressureForce(el, fe, tn); else PressureForce(el, fe, tn);

			// add element force vector to global force vector
			psolver->AssembleResidual(el.m_node, el.LM(), fe, R);
		}
	}
}

//-----------------------------------------------------------------------------

void FEConstTractionSurface::Residual(FESolidSolver* psolver, vector<double>& R)
{
	FEM& fem = psolver->m_fem;

	vector<double> fe;

	int i, n;
	int npr = m_TC.size();
	for (int iel=0; iel<npr; ++iel)
	{
		FETractionLoad& pc = m_TC[iel];
		if (pc.IsActive())
		{
			FESurfaceElement& el = m_el[iel];
			UnpackElement(el);

			double g = fem.GetLoadCurve(pc.lc)->Value();

			int ndof = 3*el.Nodes();
			fe.resize(ndof);

			// nr integration points
			int nint = el.GaussPoints();

			// nr of element nodes
			int neln = el.Nodes();

			// nodal coordinates
			vec3d *r0 = el.r0();

			double* Gr, *Gs;
			double* N;
			double* w  = el.GaussWeights();

			vec3d dxr, dxs;

			// repeat over integration points
			zero(fe);
			for (n=0; n<nint; ++n)
			{
				N  = el.H(n);
				Gr = el.Gr(n);
				Gs = el.Gs(n);

				// calculate the traction at the integration point
				vec3d t = el.eval(pc.s, n)*g;

				// calculate the tangent vectors
				dxr = dxs = vec3d(0,0,0);
				for (i=0; i<neln; ++i) 
				{
					dxr.x += Gr[i]*r0[i].x;
					dxr.y += Gr[i]*r0[i].y;
					dxr.z += Gr[i]*r0[i].z;

					dxs.x += Gs[i]*r0[i].x;
					dxs.y += Gs[i]*r0[i].y;
					dxs.z += Gs[i]*r0[i].z;
				}

				vec3d f = t*((dxr ^ dxs).norm()*w[n]);

				for (i=0; i<neln; ++i)
				{
					fe[3*i  ] += N[i]*f.x;
					fe[3*i+1] += N[i]*f.y;
					fe[3*i+2] += N[i]*f.z;
				}
			}

			// add element force vector to global force vector
			psolver->AssembleResidual(el.m_node, el.LM(), fe, R);
		}
	}
}
