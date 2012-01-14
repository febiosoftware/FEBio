#include "stdafx.h"
#include "FETractionLoad.h"
#include "FESolidSolver.h"

//-----------------------------------------------------------------------------
//! Calculate the residual for the traction load
void FETractionLoad::Residual(FENLSolver* psolver, vector<double>& R)
{
	FEModel& fem = psolver->GetFEModel();

	vector<double> fe;
	vector<int> lm;

	vec3d r0[4];

	int i, n;
	int npr = m_TC.size();
	for (int iel=0; iel<npr; ++iel)
	{
		LOAD& pc = m_TC[iel];
		FESurfaceElement& el = m_psurf->Element(iel);

		double g = fem.GetLoadCurve(pc.lc)->Value();

		int ndof = 3*el.Nodes();
		fe.resize(ndof);

		// nr integration points
		int nint = el.GaussPoints();

		// nr of element nodes
		int neln = el.Nodes();

		// nodal coordinates
		for (i=0; i<neln; ++i) r0[i] = m_psurf->GetMesh()->Node(el.m_node[i]).m_r0;

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

		// get the element's LM vector
		m_psurf->UnpackLM(el, lm);

		// add element force vector to global force vector
		psolver->AssembleResidual(el.m_node, lm, fe, R);
	}
}

//-----------------------------------------------------------------------------
void FETractionLoad::Serialize(DumpFile& ar)
{
	if (ar.IsSaving())
	{
		ar << (int) m_TC.size();
		for (int i=0; i < (int) m_TC.size(); ++i)
		{
			LOAD& d = m_TC[i];
			ar << d.lc;
			ar << d.s[0] << d.s[1] << d.s[2] << d.s[3];
		}
	}
	else
	{
		int n;
		ar >> n;
		m_TC.resize(n);
		for (int i=0; i<n; ++i)
		{
			LOAD& d = m_TC[i];
			ar >> d.lc;
			ar >> d.s[0] >> d.s[1] >> d.s[2] >> d.s[3];
		}
	}
}
