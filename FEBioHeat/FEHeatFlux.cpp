#include "FEHeatFlux.h"
#include "FECore/FEModel.h"

//-----------------------------------------------------------------------------
FEHeatFlux::LOAD::LOAD()
{ 
	s[0] = s[1] = s[2] = s[3] = s[4] = s[5] = s[6] = s[7] = s[8] = 1.0; 
	lc = -1; 
}

//-----------------------------------------------------------------------------
BEGIN_PARAMETER_LIST(FEHeatFlux, FESurfaceLoad)
	ADD_PARAMETER(m_flux, FE_PARAM_DOUBLE, "flux");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FEHeatFlux::FEHeatFlux(FEModel* pfem) : FESurfaceLoad(pfem)
{
	m_flux = 1.0;
}

//-----------------------------------------------------------------------------
//! allocate storage
void FEHeatFlux::Create(int n)
{ 
	m_FC.resize(n);
}

//-----------------------------------------------------------------------------
//! Calculate the heat flux residual
void FEHeatFlux::Residual(FEGlobalVector& R)
{
	int i, j, n;
	FEModel& fem = R.GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	const int dof_T = fem.GetDOFS().GetDOF("T");
	if (dof_T == -1) { assert(false); return; }
	
	vector<int> elm;
	int nfc = m_psurf->Elements();
	for (i=0; i<nfc; ++i)
	{
		LOAD& hf = HeatFlux(i);
		FESurfaceElement& el = m_psurf->Element(i);

		int ne = el.Nodes();
		int ni = el.GaussPoints();

		double g = m_flux;
		if (hf.lc >= 0) g *= fem.GetLoadCurve(hf.lc)->Value();

		// calculate nodal fluxes
		double qn[FEElement::MAX_NODES];
		for (j=0; j<el.Nodes(); ++j) qn[j] = g*hf.s[j];

		vector<double> fe(ne);

		// nodal coordinates
		vec3d rt[FEElement::MAX_NODES];
		for (j=0; j<ne; ++j) rt[j] = m_psurf->GetMesh()->Node(el.m_node[j]).m_rt;

		double* Gr, *Gs;
		double* N;
		double* w  = el.GaussWeights();

		// pressure at integration points
		double q;

		vec3d dxr, dxs;

		// get the element's LM vector
		vector<int> lm(ne);
		for (j=0; j<ne; ++j) lm[j] = mesh.Node(el.m_node[j]).m_ID[dof_T];

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
bool FEHeatFlux::SetFacetAttribute(int nface, const char* szatt, const char* szval)
{
	LOAD& pc = HeatFlux(nface);
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
void FEHeatFlux::Serialize(DumpStream &ar)
{
	FESurfaceLoad::Serialize(ar);
	
	if (ar.IsSaving())
	{
		ar << (int) m_FC.size();
		for (int i=0; i<(int) m_FC.size(); ++i)
		{
			LOAD& d = m_FC[i];
			ar << d.lc;
			ar << d.s[0] << d.s[1] << d.s[2] << d.s[3];
			ar << d.s[4] << d.s[5] << d.s[6] << d.s[7] << d.s[8];
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
			ar >> d.s[4] >> d.s[5] >> d.s[6] >> d.s[7] >> d.s[8];
		}
	}
}
