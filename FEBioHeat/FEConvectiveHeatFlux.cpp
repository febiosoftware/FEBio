#include "FEConvectiveHeatFlux.h"
#include "FECore/FEModel.h"

//-----------------------------------------------------------------------------
BEGIN_PARAMETER_LIST(FEConvectiveHeatFlux, FESurfaceLoad)
	ADD_PARAMETER(m_hc, FE_PARAM_DOUBLE, "hc");
	ADD_PARAMETER(m_Ta, FE_PARAM_DOUBLE, "Ta");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FEConvectiveHeatFlux::FEConvectiveHeatFlux(FEModel* pfem) : FESurfaceLoad(pfem), m_FC(FE_DOUBLE)
{
	m_hc = 0;
	m_Ta = 1.0;
	m_FC.set(1.0);
}

//-----------------------------------------------------------------------------
void FEConvectiveHeatFlux::SetSurface(FESurface* psurf)
{
	FESurfaceLoad::SetSurface(psurf);
	m_FC.Create(psurf);
}

//-----------------------------------------------------------------------------
//! residual
void FEConvectiveHeatFlux::Residual(const FETimeInfo& tp, FEGlobalVector& R)
{
	FEModel& fem = R.GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	const int dof_T = fem.GetDOFS().GetDOF("T");
	if (dof_T == -1) { assert(false); return; }

	int nfc = m_psurf->Elements();
	for (int i=0; i<nfc; ++i)
	{
		FESurfaceElement& el = m_psurf->Element(i);

		int ne = el.Nodes();
		int ni = el.GaussPoints();

		// calculate nodal fluxes
		double qn[FEElement::MAX_NODES];
		for (int j=0; j<el.Nodes(); ++j) qn[j] = m_Ta*m_FC.get<double>(i)*m_hc;

		vector<double> fe(ne);

		// nodal coordinates
		vec3d rt[FEElement::MAX_NODES];
		for (int j=0; j<ne; ++j) rt[j] = mesh.Node(el.m_node[j]).m_rt;

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

//-----------------------------------------------------------------------------
//! stiffness matrix
void FEConvectiveHeatFlux::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	const int dof_T = fem.GetDOFS().GetDOF("T");
	if (dof_T == -1) { assert(false); return; }

	matrix ke;

	int npr = m_FC.size();
	for (int m=0; m<npr; ++m)
	{
		// get the surface element
		FESurfaceElement& el = m_psurf->Element(m);

		// get the element stiffness matrix
		int neln = el.Nodes();
		int ndof = neln;
		ke.resize(ndof, ndof);

		// calculate pressure stiffness
		ElementStiffness(el, ke, m_hc);

		// get the element's LM vector
		vector<int> lm(neln);
		for (int j=0; j<neln; ++j) lm[j] = mesh.Node(el.m_node[j]).m_ID[dof_T];

		// assemble element matrix in global stiffness matrix
		LS.AssembleLHS(lm, ke);
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
	vec3d rt[FEElement::MAX_NODES];
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
