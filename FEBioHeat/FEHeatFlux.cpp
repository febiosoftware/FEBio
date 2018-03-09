#include "FEHeatFlux.h"
#include "FECore/FEModel.h"

//-----------------------------------------------------------------------------
BEGIN_PARAMETER_LIST(FEHeatFlux, FESurfaceLoad)
	ADD_PARAMETER(m_flux, FE_PARAM_DOUBLE    , "flux" );
	ADD_PARAMETER(m_FC  , FE_PARAM_DATA_ARRAY, "value");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FEHeatFlux::FEHeatFlux(FEModel* pfem) : FESurfaceLoad(pfem), m_FC(FE_DOUBLE)
{
	m_flux = 1.0;
}

//-----------------------------------------------------------------------------
//! allocate storage
void FEHeatFlux::SetSurface(FESurface* psurf)
{ 
	FESurfaceLoad::SetSurface(psurf);
	m_FC.Create(psurf, 1.0);
}

//-----------------------------------------------------------------------------
//! Calculate the heat flux residual
void FEHeatFlux::Residual(const FETimeInfo& tp, FEGlobalVector& R)
{
	FEModel& fem = R.GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	const int dof_T = fem.GetDOFS().GetDOF("T");
	if (dof_T == -1) { assert(false); return; }
	
	vector<int> elm;
	int nfc = m_psurf->Elements();
	for (int i=0; i<nfc; ++i)
	{
		FESurfaceElement& el = m_psurf->Element(i);

		int ne = el.Nodes();
		int ni = el.GaussPoints();

		// calculate nodal fluxes
		double qn[FEElement::MAX_NODES];
		for (int j=0; j<el.Nodes(); ++j) qn[j] = m_flux*m_FC.value<double>(i, j);

		vector<double> fe(ne);

		// nodal coordinates
		vec3d rt[FEElement::MAX_NODES];
		for (int j=0; j<ne; ++j) rt[j] = m_psurf->GetMesh()->Node(el.m_node[j]).m_rt;

		double* Gr, *Gs;
		double* N;
		double* w  = el.GaussWeights();

		// pressure at integration points
		double q;

		vec3d dxr, dxs;

		// get the element's LM vector
		vector<int> lm(ne);
		for (int j=0; j<ne; ++j) lm[j] = mesh.Node(el.m_node[j]).m_ID[dof_T];

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
