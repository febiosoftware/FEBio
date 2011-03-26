#include "stdafx.h"
#include "FEPressureLoad.h"
#include "FESolidSolver.h"

//-----------------------------------------------------------------------------
//! calculates the stiffness contribution due to hydrostatic pressure

void FEPressureLoad::PressureStiffness(FESurfaceElement& el, matrix& ke, vector<double>& tn)
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

bool FEPressureLoad::PressureForce(FESurfaceElement& el, vector<double>& fe, vector<double>& tn)
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

bool FEPressureLoad::LinearPressureForce(FESurfaceElement& el, vector<double>& fe, vector<double>& tn)
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

void FEPressureLoad::Serialize(FEM& fem, DumpFile& ar)
{
	if (ar.IsSaving())
	{
		m_surf.Serialize(fem, ar);

		// pressure forces
		size_t i;
		for (i=0; i<m_PC.size(); ++i)
		{
			LOAD& pc = m_PC[i];
			ar << pc.lc;
			ar << pc.s[0] << pc.s[1] << pc.s[2] << pc.s[3];
		}
	}
	else
	{
		int i, m, mat;
		m_surf.Serialize(fem, ar);
/*
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
*/

		// pressure forces
		for (i=0; i<(int) m_PC.size(); ++i)
		{
			LOAD& pc = m_PC[i];
			ar >> pc.lc;
			ar >> pc.s[0] >> pc.s[1] >> pc.s[2] >> pc.s[3];
		}

		// initialize surface data
		m_surf.Init();
	}
}

void FEPressureLoad::StiffnessMatrix(FESolver* psolver)
{
	FESolidSolver& solver = dynamic_cast<FESolidSolver&>(*psolver);
	FEM& fem = solver.m_fem;

	matrix ke;

	int npr = m_PC.size();
	for (int m=0; m<npr; ++m)
	{
		LOAD& pc = m_PC[m];
		// get the surface element
		FESurfaceElement& el = m_surf.Element(m);

		// skip rigid surface elements
		// TODO: do we really need to skip rigid elements?
		if (!el.IsRigid())
		{
			m_surf.UnpackElement(el);

			// calculate nodal normal tractions
			int neln = el.Nodes();
			vector<double> tn(neln);

			if (m_blinear == false)
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
				solver.AssembleStiffness(el.m_node, el.LM(), ke);
			}
		}
	}
}

//-----------------------------------------------------------------------------
void FEPressureLoad::Residual(FESolver* psolver, vector<double>& R)
{
	FESolidSolver& solver = dynamic_cast<FESolidSolver&>(*psolver);
	FEM& fem = solver.m_fem;

	vector<double> fe;

	int npr = m_PC.size();
	for (int i=0; i<npr; ++i)
	{
		LOAD& pc = m_PC[i];
		FESurfaceElement& el = m_surf.Element(i);
		m_surf.UnpackElement(el);

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

		if (m_blinear) LinearPressureForce(el, fe, tn); else PressureForce(el, fe, tn);

		// add element force vector to global force vector
		solver.AssembleResidual(el.m_node, el.LM(), fe, R);
	}
}
