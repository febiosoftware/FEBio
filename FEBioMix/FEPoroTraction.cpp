#include "FEPoroTraction.h"
#include "FECore/FEModel.h"

//-----------------------------------------------------------------------------
BEGIN_PARAMETER_LIST(FEPoroNormalTraction, FESurfaceLoad)
	ADD_PARAMETER(m_traction  , FE_PARAM_DOUBLE, "traction" );
	ADD_PARAMETER(m_blinear   , FE_PARAM_BOOL  , "linear"   );
    ADD_PARAMETER(m_bshellb , FE_PARAM_BOOL  , "shell_bottom");
	ADD_PARAMETER(m_beffective, FE_PARAM_BOOL  , "effective");
	ADD_PARAMETER(m_PC        , FE_PARAM_DATA_ARRAY, "value");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! constructor
FEPoroNormalTraction::FEPoroNormalTraction(FEModel* pfem) : FESurfaceLoad(pfem), m_PC(FE_DOUBLE)
{ 
	m_traction = 1.0;
	m_blinear = false; 
    m_bshellb = false;
	m_beffective = false;
	m_PC.set(1.0);

	// get the degrees of freedom
	m_dofX = pfem->GetDOFIndex("x");
	m_dofY = pfem->GetDOFIndex("y");
	m_dofZ = pfem->GetDOFIndex("z");
	m_dofP = pfem->GetDOFIndex("p");
    m_dofU = pfem->GetDOFIndex("u");
    m_dofV = pfem->GetDOFIndex("v");
    m_dofW = pfem->GetDOFIndex("w");
    m_dofQ = pfem->GetDOFIndex("q");
}

//-----------------------------------------------------------------------------
//! allocate storage
void FEPoroNormalTraction::SetSurface(FESurface* ps)
{ 
	FESurfaceLoad::SetSurface(ps);
	m_PC.Create(ps); 
}

//-----------------------------------------------------------------------------
void FEPoroNormalTraction::UnpackLM(FEElement& el, vector<int>& lm)
{
	FEMesh& mesh = *GetSurface().GetMesh();
	int N = el.Nodes();
	lm.resize(N*4);
    if (!m_bshellb) {
        for (int i=0; i<N; ++i)
        {
            int n = el.m_node[i];
            FENode& node = mesh.Node(n);
            vector<int>& id = node.m_ID;
            
            // first the displacement dofs
            lm[3*i  ] = id[m_dofX];
            lm[3*i+1] = id[m_dofY];
            lm[3*i+2] = id[m_dofZ];
            
            // now the pressure dofs
            lm[3*N+i] = id[m_dofP];
        }
    }
    else {
        for (int i=0; i<N; ++i)
        {
            int n = el.m_node[i];
            FENode& node = mesh.Node(n);
            vector<int>& id = node.m_ID;
            
            // first the displacement dofs
            lm[3*i  ] = id[m_dofU];
            lm[3*i+1] = id[m_dofV];
            lm[3*i+2] = id[m_dofW];
            
            // now the pressure dofs
            lm[3*N+i] = id[m_dofQ];
        }
    }
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
    if (m_bshellb) {
        for (int j=0; j<neln; ++j) {
            FENode& nd = m_psurf->GetMesh()->Node(el.m_node[j]);
            rt[j] -= nd.m_d0 + nd.get_vec3d(m_dofX, m_dofY, m_dofZ) - nd.get_vec3d(m_dofU, m_dofV, m_dofW);
            
        }
    }

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
        if (m_bshellb) tr = -tr;
		
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
    if (m_bshellb) {
        for (int j=0; j<neln; ++j) {
            FENode& nd = m_psurf->GetMesh()->Node(el.m_node[j]);
            rt[j] -= nd.m_d0 + nd.get_vec3d(m_dofX, m_dofY, m_dofZ) - nd.get_vec3d(m_dofU, m_dofV, m_dofW);
            
        }
    }

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
        if (m_bshellb) tr = -tr;

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
    if (m_bshellb) {
        for (int j=0; j<neln; ++j) {
            FENode& nd = m_psurf->GetMesh()->Node(el.m_node[j]);
            r0[j] -= nd.m_d0;
        }
    }

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
        if (m_bshellb) tr = -tr;

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
void FEPoroNormalTraction::StiffnessMatrix(const FETimeInfo& tp, FESolver* psolver)
{
	FEMesh& mesh = *GetSurface().GetMesh();

	matrix ke;

	vector<int> lm;

	int npr = m_PC.size();
	for (int m=0; m<npr; ++m)
	{
		// get the surface element
		FESurfaceElement& el = m_psurf->Element(m);
		int neln = el.Nodes();

		// fluid pressure
		double pt[FEElement::MAX_NODES];
        if (!m_bshellb) {
            for (int i=0; i<neln; ++i) pt[i] = mesh.Node(el.m_node[i]).get(m_dofP);
        }
        else {
            for (int i=0; i<neln; ++i) pt[i] = mesh.Node(el.m_node[i]).get(m_dofQ);
        }
			
		// calculate nodal normal tractions
		vector<double> tn(neln);

		if (m_blinear == false)
		{
			// evaluate the prescribed traction.
			for (int j=0; j<neln; ++j) tn[j] = m_traction*m_PC.get<double>(m);

			// if the prescribed traction is effective, evaluate the total traction
			if (m_beffective) for (int j=0; j<neln; ++j) tn[j] -= pt[j];
				
			// get the element stiffness matrix
			int ndof = (m_beffective ? 4*neln : 3*neln);
			ke.resize(ndof, ndof);

			// calculate pressure stiffness
			TractionStiffness(el, ke, tn, m_beffective, psolver->m_bsymm);

			// get the element's LM vector
			UnpackLM(el, lm);

			// assemble element matrix in global stiffness matrix
			psolver->AssembleStiffness(el.m_node, lm, ke);
		}
	}
}

//-----------------------------------------------------------------------------
void FEPoroNormalTraction::Residual(const FETimeInfo& tp, FEGlobalVector& R)
{
	FEMesh& mesh = *GetSurface().GetMesh();

	vector<double> fe;

	vector<int> lm;

	int npr = m_PC.size();
	for (int i=0; i<npr; ++i)
	{
		FESurfaceElement& el = m_psurf->Element(i);
		int neln = el.Nodes();

		// fluid pressure
		double pt[FEElement::MAX_NODES];
        if (!m_bshellb) {
            for (int j=0; j<neln; ++j) pt[j] = mesh.Node(el.m_node[j]).get(m_dofP);
        }
        else {
            for (int j=0; j<neln; ++j) pt[j] = mesh.Node(el.m_node[j]).get(m_dofQ);
        }

		// calculate nodal normal tractions
		vector<double> tn(neln);

		// evaluate the prescribed traction.
		for (int j=0; j<neln; ++j) tn[j] = m_traction*m_PC.get<double>(i);
		
		// if the prescribed traction is effective, evaluate the total traction
		if (m_beffective) for (int j=0; j<neln; ++j) tn[j] -= pt[j];

		int ndof = (m_beffective? 4*neln : 3*neln);
		fe.resize(ndof);

		if (m_blinear) LinearTractionForce(el, fe, tn); else TractionForce(el, fe, tn);

		// get the element's LM vector
		UnpackLM(el, lm);

		// add element force vector to global force vector
		R.Assemble(el.m_node, lm, fe);
	}
}
