#include "stdafx.h"
#include "FEHeatFlux.h"
#include "FESolver.h"

//-----------------------------------------------------------------------------
//! Calculate the heat flux residual
void FEHeatFlux::Residual(FESolver* psolver, vector<double>& R)
{
	int i, j, n;
	FEM& fem = dynamic_cast<FEM&>(psolver->GetFEModel());
	vector<int> elm;

	int nfc = m_psurf->Elements();
	for (i=0; i<nfc; ++i)
	{
		LOAD& hf = HeatFlux(i);
		FESurfaceElement& el = m_psurf->Element(i);

		int ne = el.Nodes();
		int ni = el.GaussPoints();

		double g = fem.GetLoadCurve(hf.lc)->Value();

		// calculate nodal fluxes
		double qn[4];
		for (j=0; j<el.Nodes(); ++j) qn[j] = g*hf.s[j];

		vector<double> fe(ne);

		// nodal coordinates
		vec3d rt[4];
		for (j=0; j<ne; ++j) rt[j] = m_psurf->GetMesh()->Node(el.m_node[j]).m_rt;

		double* Gr, *Gs;
		double* N;
		double* w  = el.GaussWeights();

		// pressure at integration points
		double q;

		vec3d dxr, dxs;

		// get the element's LM vector
		m_psurf->UnpackLM(el, elm);

		vector<int> lm(ne);
		for (j=0; j<ne; ++j) lm[j] = elm[ne*10 + j];

		// force vector
		// repeat over integration points
		zero(fe);
		for (n=0; n<ni; ++n)
		{
			N  = el.H(n);
			Gr = el.Gr(n);
			Gs = el.Gs(n);

			q = 0;
			dxr = dxs = vec3d(0,0,0);
			for (j=0; j<ne; ++j) 
			{
				q += N[j]*qn[j];
				dxr.x += Gr[j]*rt[j].x;
				dxr.y += Gr[j]*rt[j].y;
				dxr.z += Gr[j]*rt[j].z;

				dxs.x += Gs[j]*rt[j].x;
				dxs.y += Gs[j]*rt[j].y;
				dxs.z += Gs[j]*rt[j].z;
			}
	
			double J = (dxr ^ dxs).norm();

			for (j=0; j<ne; ++j) fe[j] += N[j]*q*J*w[n];
		}

		// add element force vector to global force vector
		for (j=0; j<ne; ++j)
		{
			if (lm[j] >= 0) R[lm[j]] += fe[j];
		}
	}
}

//-----------------------------------------------------------------------------
void FEHeatFlux::Serialize(DumpFile &ar)
{
	if (ar.IsSaving())
	{
		ar << (int) m_FC.size();
		for (int i=0; i<(int) m_FC.size(); ++i)
		{
			LOAD& d = m_FC[i];
			ar << d.lc;
			ar << d.s[0] << d.s[1] << d.s[2] << d.s[3];
		}
	}
	else
	{
		int n;
		ar >> n;
		m_FC.resize(n);
		for (int i=0; i<n; ++i)
		{
			LOAD& d = m_FC[i];
			ar >> d.lc;
			ar >> d.s[0] >> d.s[1] >> d.s[2] >> d.s[3];
		}
	}
}
