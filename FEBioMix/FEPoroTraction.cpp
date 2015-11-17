#include "FEPoroTraction.h"
#include "FECore/FEModel.h"

//-----------------------------------------------------------------------------
FEPoroNormalTraction::LOAD::LOAD()
{ 
	s[0] = s[1] = s[2] = s[3] = s[4] = s[5] = s[6] = s[7] = s[8] = 1.0; 
	lc = -1;
}

//-----------------------------------------------------------------------------
BEGIN_PARAMETER_LIST(FEPoroNormalTraction, FESurfaceLoad)
	ADD_PARAMETER(m_traction  , FE_PARAM_DOUBLE, "traction" );
	ADD_PARAMETER(m_blinear   , FE_PARAM_BOOL  , "linear"   );
	ADD_PARAMETER(m_beffective, FE_PARAM_BOOL  , "effective");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! constructor
FEPoroNormalTraction::FEPoroNormalTraction(FEModel* pfem) : FESurfaceLoad(pfem)
{ 
	m_traction = 1.0;
	m_blinear = false; 
	m_beffective = false; 
}

//-----------------------------------------------------------------------------
//! allocate storage
void FEPoroNormalTraction::Create(int n)
{ 
	m_PC.resize(n); 
}

//-----------------------------------------------------------------------------
//! \deprecated This function is only needed for the 1.2 file format which is obsolete
bool FEPoroNormalTraction::SetAttribute(const char* szatt, const char* szval)
{
	if (strcmp(szatt, "type") == 0)
	{
		if      (strcmp(szval, "linear"   ) == 0) SetLinear(true );
		else if (strcmp(szval, "nonlinear") == 0) SetLinear(false);
		else return false;
		return true;
	}
	else if (strcmp(szatt, "traction") == 0)
	{
		if      (strcmp(szval, "effective") == 0) SetEffective(true);
		else if (strcmp(szval, "total"    ) == 0) SetEffective(false);
		else if (strcmp(szval, "mixture"  ) == 0) SetEffective(false);
		else return false;
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
bool FEPoroNormalTraction::SetFacetAttribute(int nface, const char* szatt, const char* szval)
{
	LOAD& pc = NormalTraction(nface);
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
//! calculates the stiffness contribution due to normal traction
void FEPoroNormalTraction::TractionStiffness(FESurfaceElement& el, matrix& ke, vector<double>& tn, bool effective, bool bsymm)
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
	vec3d rt[FEElement::MAX_NODES];
	for (j=0; j<neln; ++j) rt[j] = m_psurf->GetMesh()->Node(el.m_node[j]).m_rt;

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
		if (!bsymm) {
			// non-symmetric
			for (i=0; i<neln; ++i)
				for (j=0; j<neln; ++j)
				{
					kab = (dxs*Gr[j] - dxr*Gs[j])*N[i]*w[n]*tr;
					
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
			
			// if prescribed traction is effective, add stiffness component
			if (effective)
			{
				for (i=0; i<neln; ++i)
					for (j=0; j<neln; ++j)
					{
						kab = (dxr ^ dxs)*w[n]*N[i]*N[j];
						
						ke[3*i  ][3*neln+j] += kab.x;
						ke[3*i+1][3*neln+j] += kab.y;
						ke[3*i+2][3*neln+j] += kab.z;
					}
			}
		} else {
			// symmetric
			for (i=0; i<neln; ++i)
				for (j=0; j<neln; ++j)
				{
					kab = ((dxs*Gr[j] - dxr*Gs[j])*N[i]-(dxs*Gr[i] - dxr*Gs[i])*N[j])*0.5*w[n]*tr;
					
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
			
			// if prescribed traction is effective, add stiffness component
			if (effective)
			{
				for (i=0; i<neln; ++i)
					for (j=0; j<neln; ++j)
					{
						kab = (dxr ^ dxs)*w[n]*0.5*N[i]*N[j];
						
						ke[3*i  ][3*neln+j] += kab.x;
						ke[3*i+1][3*neln+j] += kab.y;
						ke[3*i+2][3*neln+j] += kab.z;

						ke[3*i  ][3*neln+j] += kab.x;
						ke[3*i+1][3*neln+j] += kab.y;
						ke[3*i+2][3*neln+j] += kab.z;
					}
			}
		}

	}
}

//-----------------------------------------------------------------------------
//! calculates the equivalent nodal forces due to hydrostatic pressure

bool FEPoroNormalTraction::TractionForce(FESurfaceElement& el, vector<double>& fe, vector<double>& tn)
{
	int i, n;

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

bool FEPoroNormalTraction::LinearTractionForce(FESurfaceElement& el, vector<double>& fe, vector<double>& tn)
{
	int i, n;

	// nr integration points
	int nint = el.GaussPoints();

	// nr of element nodes
	int neln = el.Nodes();
	assert(neln <= 4);

	// nodal coordinates
	vec3d r0[FEElement::MAX_NODES];
	for (i=0; i<neln; ++i) r0[i] = m_psurf->GetMesh()->Node(el.m_node[i]).m_r0;

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

void FEPoroNormalTraction::Serialize(DumpFile& ar)
{
	FESurfaceLoad::Serialize(ar);

	if (ar.IsSaving())
	{
		ar << m_traction;
		ar << m_blinear << m_beffective;
		ar << (int) m_PC.size();
		for (int i=0; i<(int) m_PC.size(); ++i)
		{
			LOAD& pc = m_PC[i];
			ar << pc.lc;
			ar << pc.s[0] << pc.s[1] << pc.s[2] << pc.s[3];
			ar << pc.s[4] << pc.s[5] << pc.s[6] << pc.s[7] << pc.s[8];
		}
	}
	else
	{
		int n;
		ar >> m_traction;
		ar >> m_blinear >> m_beffective;
		ar >> n;
		m_PC.resize(n);
		for (int i=0; i<n; ++i)
		{
			LOAD& pc = m_PC[i];
			ar >> pc.lc;
			ar >> pc.s[0] >> pc.s[1] >> pc.s[2] >> pc.s[3];
			ar >> pc.s[4] >> pc.s[5] >> pc.s[6] >> pc.s[7] >> pc.s[8];
		}
	}
}

//-----------------------------------------------------------------------------
void FEPoroNormalTraction::StiffnessMatrix(FESolver* psolver)
{
	FEModel& fem = psolver->GetFEModel();

	matrix ke;

	vector<int> lm;

	int npr = m_PC.size();
	for (int m=0; m<npr; ++m)
	{
		LOAD& pc = m_PC[m];
		// get the surface element
		FESurfaceElement& el = m_psurf->Element(m);
		int neln = el.Nodes();

		// fluid pressure
		double pt[FEElement::MAX_NODES];
		for (int i=0; i<neln; ++i) pt[i] = m_psurf->GetMesh()->Node(el.m_node[i]).m_pt;
			
		// calculate nodal normal tractions
		vector<double> tn(neln);

		if (m_blinear == false)
		{
			double g = m_traction;
			if (pc.lc >= 0) g *= fem.GetLoadCurve(pc.lc)->Value();

			// evaluate the prescribed traction.
			for (int j=0; j<neln; ++j) tn[j] = g*pc.s[j];

			// if the prescribed traction is effective, evaluate the total traction
			if (m_beffective) for (int j=0; j<neln; ++j) tn[j] -= pt[j];
				
			// get the element stiffness matrix
			int ndof = (m_beffective ? 4*neln : 3*neln);
			ke.resize(ndof, ndof);

			// calculate pressure stiffness
			TractionStiffness(el, ke, tn, m_beffective, psolver->m_bsymm);

			// get the element's LM vector
			m_psurf->UnpackLM(el, lm);

			// assemble element matrix in global stiffness matrix
			psolver->AssembleStiffness(el.m_node, lm, ke);
		}
	}
}

//-----------------------------------------------------------------------------
void FEPoroNormalTraction::Residual(FEGlobalVector& R)
{
	FEModel& fem = R.GetFEModel();

	vector<double> fe;

	vector<int> lm;

	int npr = m_PC.size();
	for (int i=0; i<npr; ++i)
	{
		LOAD& pc = m_PC[i];
		FESurfaceElement& el = m_psurf->Element(i);
		int neln = el.Nodes();

		// fluid pressure
		double pt[FEElement::MAX_NODES];
		for (int j=0; j<neln; ++j) pt[j] = m_psurf->GetMesh()->Node(el.m_node[j]).m_pt;

		// calculate nodal normal tractions
		vector<double> tn(neln);

		double g = m_traction;
		if (pc.lc >= 0) g *= fem.GetLoadCurve(pc.lc)->Value();

		// evaluate the prescribed traction.
		for (int j=0; j<neln; ++j) tn[j] = g*pc.s[j];
		
		// if the prescribed traction is effective, evaluate the total traction
		if (m_beffective) for (int j=0; j<neln; ++j) tn[j] -= pt[j];

		int ndof = (m_beffective? 4*neln : 3*neln);
		fe.resize(ndof);

		if (m_blinear) LinearTractionForce(el, fe, tn); else TractionForce(el, fe, tn);

		// get the element's LM vector
		m_psurf->UnpackLM(el, lm);

		// add element force vector to global force vector
		R.Assemble(el.m_node, lm, fe);
	}
}
