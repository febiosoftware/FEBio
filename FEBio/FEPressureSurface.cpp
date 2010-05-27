#include "stdafx.h"
#include "FEPressureSurface.h"
#include "FESolidSolver.h"

//-----------------------------------------------------------------------------
//! calculates the stiffness contribution due to hydrostatic pressure

void FEPressureSurface::PressureStiffness(FESurfaceElement& el, matrix& ke)
{
	int i, j, n;

	int nint = el.GaussPoints();
	int neln = el.Nodes();

	// pressure at nodes
	double *pn = el.pt();

	// pressure at integration point
	double p;

	// gauss weights
	double* w = el.GaussWeights();

	// nodal coordinates
	vec3d* rt = el.rt();

	// jacobian
	double J[3][2];

	double t1, t2;
	double kab[3];

	ke.zero();

	double* N, *Gr, *Gs;

	// repeat over integration points
	for (n=0; n<nint; ++n)
	{
		N = el.H(n);
		Gr = el.Gr(n);
		Gs = el.Gs(n);

		// calculate pressure at integration point
		p = 0;
		for (i=0; i<neln; ++i) p += N[i]*pn[i];

		// calculate jacobian
		J[0][0] = J[0][1] = 0;
		J[1][0] = J[1][1] = 0;
		J[2][0] = J[2][1] = 0;
		for (i=0; i<neln; ++i)
		{
			J[0][0] += Gr[i]*rt[i].x; J[0][1] += Gs[i]*rt[i].x;
			J[1][0] += Gr[i]*rt[i].y; J[1][1] += Gs[i]*rt[i].y;
			J[2][0] += Gr[i]*rt[i].z; J[2][1] += Gs[i]*rt[i].z;
		}

		// calculate stiffness component
		for (i=0; i<neln; ++i)
			for (j=0; j<neln; ++j)
			{
				t1 = 0.5*(Gs[i]*N[j] - Gs[j]*N[i]);
				t2 = 0.5*(Gr[i]*N[j] - Gr[j]*N[i]);

				kab[0] = p*(J[0][0]*t1 - J[0][1]*t2)*w[n];
				kab[1] = p*(J[1][0]*t1 - J[1][1]*t2)*w[n];
				kab[2] = p*(J[2][0]*t1 - J[2][1]*t2)*w[n];

				ke[3*i  ][3*j  ] +=       0; //(0,0,0)*kab[0]+(0,0,1)*kab[1]+(0,0,2)*kab[2];
				ke[3*i  ][3*j+1] +=  kab[2]; //(0,1,0)*kab[0]+(0,1,1)*kab[1]+(0,1,2)*kab[2];
				ke[3*i  ][3*j+2] += -kab[1]; //(0,2,0)*kab[0]+(0,2,1)*kab[1]+(0,2,2)*kab[2];

				ke[3*i+1][3*j  ] += -kab[2]; //(1,0,0)*kab[0]+(1,0,1)*kab[1]+(1,0,2)*kab[2];
				ke[3*i+1][3*j+1] +=       0; //(1,1,0)*kab[0]+(1,1,1)*kab[1]+(1,1,2)*kab[2];
				ke[3*i+1][3*j+2] +=  kab[0]; //(1,2,0)*kab[0]+(1,2,1)*kab[1]+(1,2,2)*kab[2];

				ke[3*i+2][3*j  ] +=  kab[1]; //(2,0,0)*kab[0]+(2,0,1)*kab[1]+(2,0,2)*kab[2];
				ke[3*i+2][3*j+1] += -kab[0]; //(2,1,0)*kab[0]+(2,1,1)*kab[1]+(2,1,2)*kab[2];
				ke[3*i+2][3*j+2] +=       0; //(2,2,0)*kab[0]+(2,2,1)*kab[1]+(2,2,2)*kab[2];
			}
	}
}

//-----------------------------------------------------------------------------
//! calculates the equivalent nodal forces due to hydrostatic pressure

bool FEPressureSurface::PressureForce(FESurfaceElement& el, vector<double>& fe)
{
	int i, n;

	// nr integration points
	int nint = el.GaussPoints();

	// nr of element nodes
	int neln = el.Nodes();

	// pressure at nodes
	double *pn = el.pt();

	// nodal coordinates
	vec3d *rt = el.rt();

	double* Gr, *Gs;
	double* N;
	double* w  = el.GaussWeights();

	// pressure at integration points
	double pr;

	vec3d dxr, dxs;

	// force vector
	vec3d f;

	// repeat over integration points
	fe.zero();
	for (n=0; n<nint; ++n)
	{
		N  = el.H(n);
		Gr = el.Gr(n);
		Gs = el.Gs(n);

		pr = 0;
		dxr = dxs = vec3d(0,0,0);
		for (i=0; i<neln; ++i) 
		{
			pr += N[i]*pn[i];

			dxr.x += Gr[i]*rt[i].x;
			dxr.y += Gr[i]*rt[i].y;
			dxr.z += Gr[i]*rt[i].z;

			dxs.x += Gs[i]*rt[i].x;
			dxs.y += Gs[i]*rt[i].y;
			dxs.z += Gs[i]*rt[i].z;
		}

		f = (dxr ^ dxs)*pr*w[n];

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

bool FEPressureSurface::LinearPressureForce(FESurfaceElement& el, vector<double>& fe)
{
	int i, n;

	// nr integration points
	int nint = el.GaussPoints();

	// nr of element nodes
	int neln = el.Nodes();

	// pressure at nodes
	double *pn = el.pt();

	// nodal coordinates
	vec3d *r0 = el.r0();

	double* Gr, *Gs;
	double* N;
	double* w  = el.GaussWeights();

	// pressure at integration points
	double pr;

	vec3d dxr, dxs;

	// force vector
	vec3d f;

	// repeat over integration points
	fe.zero();
	for (n=0; n<nint; ++n)
	{
		N  = el.H(n);
		Gr = el.Gr(n);
		Gs = el.Gs(n);

		pr = 0;
		dxr = dxs = vec3d(0,0,0);
		for (i=0; i<neln; ++i) 
		{
			pr += N[i]*pn[i];

			dxr.x += Gr[i]*r0[i].x;
			dxr.y += Gr[i]*r0[i].y;
			dxr.z += Gr[i]*r0[i].z;

			dxs.x += Gs[i]*r0[i].x;
			dxs.y += Gs[i]*r0[i].y;
			dxs.z += Gs[i]*r0[i].z;
		}

		f = (dxr ^ dxs)*pr*w[n];

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

		// pressure forces
		for (i=0; i<m_PC.size(); ++i)
		{
			FEPressureLoad& pc = m_PC[i];
			ar << pc.blinear << pc.face << pc.lc;
			ar << pc.s[0] << pc.s[1] << pc.s[2] << pc.s[3];
			ar << pc.bc;
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

		// pressure forces
		for (i=0; i<m_PC.size(); ++i)
		{
			FEPressureLoad& pc = m_PC[i];
			ar >> pc.blinear >> pc.face >> pc.lc;
			ar >> pc.s[0] >> pc.s[1] >> pc.s[2] >> pc.s[3];
			ar >> pc.bc;
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
		if ((pc.bc == 0) && pc.IsActive())
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

				if (!pc.blinear)
				{
					double g = fem.GetLoadCurve(pc.lc)->Value();

					for (int j=0; j<el.Nodes(); ++j) pt[j] = -g*pc.s[j];

					// get the element stiffness matrix
					int ndof = 3*el.Nodes();
					ke.Create(ndof, ndof);

					// calculate pressure stiffness
					PressureStiffness(el, ke);

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
		if ((pc.bc == 0) && pc.IsActive())
		{
			FESurfaceElement& el = m_el[i];
			UnpackElement(el);

			// calculate nodal pressures
			double* pt = el.pt();

			double g = fem.GetLoadCurve(pc.lc)->Value();

			for (int j=0; j<el.Nodes(); ++j) pt[j] = -g*pc.s[j];

			int ndof = 3*el.Nodes();
			fe.resize(ndof);

			if (pc.blinear) LinearPressureForce(el, fe); else PressureForce(el, fe);

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
			fe.zero();
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
