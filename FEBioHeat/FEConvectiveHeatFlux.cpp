#include "FEConvectiveHeatFlux.h"
#include "FECore/FEModel.h"

//-----------------------------------------------------------------------------
//! residual
void FEConvectiveHeatFlux::Residual(FEGlobalVector& R)
{
	FEModel& fem = R.GetFEModel();
	vector<int> elm;

	int nfc = m_psurf->Elements();
	for (int i=0; i<nfc; ++i)
	{
		LOAD& hf = HeatFlux(i);
		FESurfaceElement& el = m_psurf->Element(i);

		int ne = el.Nodes();
		int ni = el.GaussPoints();

		// get ambient temperature
		double Tc = fem.GetLoadCurve(hf.lc)->Value();

		// calculate nodal fluxes
		double qn[4];
		for (int j=0; j<el.Nodes(); ++j) qn[j] = Tc*hf.s[j]*hf.hc;

		vector<double> fe(ne);

		// nodal coordinates
		vec3d rt[4];
		for (int j=0; j<ne; ++j) rt[j] = m_psurf->GetMesh()->Node(el.m_node[j]).m_rt;

		double* Gr, *Gs;
		double* N;
		double* w  = el.GaussWeights();

		// pressure at integration points
		double q;

		vec3d dxr, dxs;

		// get the element's LM vector
		m_psurf->UnpackLM(el, elm);

		vector<int> lm(ne);
		for (int j=0; j<ne; ++j) lm[j] = elm[ne*10 + j];

		// force vector
		// repeat over integration points
		zero(fe);
		for (int n=0; n<ni; ++n)
		{
			N  = el.H(n);
			Gr = el.Gr(n);
			Gs = el.Gs(n);

			q = 0;
			dxr = dxs = vec3d(0,0,0);
			for (int j=0; j<ne; ++j) 
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

			for (int j=0; j<ne; ++j) fe[j] += N[j]*q*J*w[n];
		}

		// add element force vector to global force vector
		for (int j=0; j<ne; ++j)
		{
			if (lm[j] >= 0) R[lm[j]] += fe[j];
		}
	}
}

//-----------------------------------------------------------------------------
//! stiffness matrix
void FEConvectiveHeatFlux::StiffnessMatrix(FESolver* psolver)
{
	FEModel& fem = psolver->GetFEModel();

	matrix ke;
	vector<int> elm;

	int npr = m_FC.size();
	for (int m=0; m<npr; ++m)
	{
		LOAD& fc = m_FC[m];

		// get the surface element
		FESurfaceElement& el = m_psurf->Element(m);

		// get the element stiffness matrix
		int neln = el.Nodes();
		int ndof = neln;
		ke.resize(ndof, ndof);

		// calculate pressure stiffness
		ElementStiffness(el, ke, fc.hc);

		// get the element's LM vector
		m_psurf->UnpackLM(el, elm);

		vector<int> lm(neln);
		for (int j=0; j<neln; ++j) lm[j] = elm[neln*10 + j];

		// assemble element matrix in global stiffness matrix
		psolver->AssembleStiffness(el.m_node, lm, ke);
	}
}

//-----------------------------------------------------------------------------
void FEConvectiveHeatFlux::ElementStiffness(FESurfaceElement& el, matrix& ke, double hc)
{
	// nr integration points
	int nint = el.GaussPoints();

	// nr of element nodes
	int neln = el.Nodes();

	// nodal coordinates
	vec3d rt[4];
	for (int j=0; j<neln; ++j) rt[j] = m_psurf->GetMesh()->Node(el.m_node[j]).m_rt;

	double* Gr, *Gs;
	double* N;
	double* w  = el.GaussWeights();

	vec3d dxr, dxs;
	ke.zero();

	// repeat over integration points
	for (int n=0; n<nint; ++n)
	{
		N  = el.H(n);
		Gr = el.Gr(n);
		Gs = el.Gs(n);

		dxr = dxs = vec3d(0,0,0);
		for (int i=0; i<neln; ++i) 
		{
			dxr += rt[i]*Gr[i];
			dxs += rt[i]*Gs[i];
		}

		double J = (dxr ^ dxs).norm();

		for (int i=0; i<neln; ++i)
		{
			for (int j=0; j<neln; ++j)
			{
				double kij = hc*N[i]*N[j]*J*w[n];
				ke[i][j] += kij;
			}
		}
	}
}

//-----------------------------------------------------------------------------
bool FEConvectiveHeatFlux::SetFacetAttribute(int nface, const char* szatt, const char* szval)
{
	LOAD& pc = HeatFlux(nface);
	if      (strcmp(szatt, "id") == 0) {}
	else if (strcmp(szatt, "lc") == 0) pc.lc = atoi(szval) - 1;
	else if (strcmp(szatt, "hc") == 0)
	{
		double hc = atof(szval);
		pc.hc = hc;
	}
	else if (strcmp(szatt, "scale") == 0)
	{
		double s = atof(szval);
		pc.s[0] = pc.s[1] = pc.s[2] = pc.s[3] = s;
		pc.s[4] = pc.s[5] = pc.s[6] = pc.s[7] = s;
	}
	else return false;

	return true;
}
