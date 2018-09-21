#include "stdafx.h"
#include "FEPressureLoad.h"
#include "FECore/FEModel.h"

//-----------------------------------------------------------------------------
// Parameter block for pressure loads
BEGIN_PARAMETER_LIST(FEPressureLoad, FESurfaceLoad)
	ADD_PARAMETER(m_blinear , FE_PARAM_BOOL  , "linear"  );
    ADD_PARAMETER(m_bshellb , FE_PARAM_BOOL  , "shell_bottom");
	ADD_PARAMETER(m_pressure, FE_PARAM_DOUBLE_MAPPED, "pressure");
	ADD_PARAMETER(m_bsymm   , FE_PARAM_BOOL  , "symmetric_stiffness");
	ADD_PARAMETER(m_bstiff  , FE_PARAM_BOOL  , "pressure_stiffness");
END_PARAMETER_LIST()

//-----------------------------------------------------------------------------
//! constructor
FEPressureLoad::FEPressureLoad(FEModel* pfem) : FESurfaceLoad(pfem)
{ 
	m_blinear = false;
    m_bshellb = false;
	m_pressure.setValue(0.0);
	m_bsymm = true;
	m_bstiff = true;

	// get the degrees of freedom
	m_dofX = pfem->GetDOFIndex("x");
	m_dofY = pfem->GetDOFIndex("y");
	m_dofZ = pfem->GetDOFIndex("z");
    m_dofSX = pfem->GetDOFIndex("sx");
    m_dofSY = pfem->GetDOFIndex("sy");
    m_dofSZ = pfem->GetDOFIndex("sz");
}

//-----------------------------------------------------------------------------
//! allocate storage
void FEPressureLoad::SetSurface(FESurface* ps)
{
	FESurfaceLoad::SetSurface(ps);
	m_pressure.addDomain(ps);
}

//-----------------------------------------------------------------------------
//! calculates the stiffness contribution due to hydrostatic pressure
void FEPressureLoad::PressureStiffness(FESurfaceElement& el, matrix& ke)
{
	// choose the symmetric of unsymmetric formulation
	if (m_bsymm) 
		SymmetricPressureStiffness(el, ke);
	else
		UnsymmetricPressureStiffness(el, ke);
}

//-----------------------------------------------------------------------------
//! calculates the stiffness contribution due to hydrostatic pressure
void FEPressureLoad::SymmetricPressureStiffness(FESurfaceElement& el, matrix& ke)
{
	int nint = el.GaussPoints();
	int neln = el.Nodes();

	// gauss weights
	double* w = el.GaussWeights();

	// nodal coordinates
	FEMesh& mesh = *m_psurf->GetMesh();
	vec3d rt[FEElement::MAX_NODES];
	for (int j = 0; j < neln; ++j)
	{
		rt[j] = mesh.Node(el.m_node[j]).m_rt;
	}
    if (m_bshellb) {
        for (int j=0; j<neln; ++j) {
            FENode& nd = mesh.Node(el.m_node[j]);
			rt[j] -= nd.m_d0 + nd.get_vec3d(m_dofX, m_dofY, m_dofZ) - nd.get_vec3d(m_dofSX, m_dofSY, m_dofSZ);
        }
    }

	// repeat over integration points
	ke.zero();
	for (int n=0; n<nint; ++n)
	{
		FEMaterialPoint& pt = *el.GetMaterialPoint(n);

		double* N = el.H(n);
		double* Gr = el.Gr(n);
		double* Gs = el.Gs(n);

		// traction at integration point
		vec3d dxr(0,0,0), dxs(0,0,0);
		for (int i=0; i<neln; ++i)
		{
			dxr += rt[i]*Gr[i];
			dxs += rt[i]*Gs[i];
		}

		// evaluate pressure at this material point
		double P = -m_pressure.eval(pt);

        if (m_bshellb) P = -P;
		
		// calculate stiffness component
		for (int i=0; i<neln; ++i)
			for (int j=0; j<neln; ++j)
			{
				vec3d kab = (dxr*(N[j]*Gs[i]-N[i]*Gs[j])
					   -dxs*(N[j]*Gr[i]-N[i]*Gr[j]))*w[n]*0.5*P;

				ke.add(3*i, 3*j, mat3da(kab));
			}
	}
}

//-----------------------------------------------------------------------------
//! calculates the stiffness contribution due to hydrostatic pressure

void FEPressureLoad::UnsymmetricPressureStiffness(FESurfaceElement& el, matrix& ke)
{
	int nint = el.GaussPoints();
	int neln = el.Nodes();

	// gauss weights
	double* w = el.GaussWeights();

	// nodal coordinates
	FEMesh& mesh = *m_psurf->GetMesh();
	vec3d rt[FEElement::MAX_NODES];
	for (int j = 0; j < neln; ++j)
	{
		rt[j] = mesh.Node(el.m_node[j]).m_rt;
	}
	if (m_bshellb) {
        for (int j=0; j<neln; ++j) {
            FENode& nd = mesh.Node(el.m_node[j]);
			rt[j] -= nd.m_d0 + nd.get_vec3d(m_dofX, m_dofY, m_dofZ) - nd.get_vec3d(m_dofSX, m_dofSY, m_dofSZ);
        }
    }

	// repeat over integration points
	ke.zero();
	for (int n=0; n<nint; ++n)
	{
		FEMaterialPoint& pt = *el.GetMaterialPoint(n);

		double* N = el.H(n);
		double* Gr = el.Gr(n);
		double* Gs = el.Gs(n);

		// traction at integration point
		vec3d dxr(0,0,0), dxs(0,0,0);
		for (int i=0; i<neln; ++i) 
		{
			dxr += rt[i]*Gr[i];
			dxs += rt[i]*Gs[i];
		}

		// evaluate pressure at this material point
		double P = -m_pressure.eval(pt);

        if (m_bshellb) P = -P;
		
		// calculate stiffness component
		for (int i=0; i<neln; ++i)
			for (int j=0; j<neln; ++j)
			{
				vec3d Kab = (dxr*Gs[j] - dxs*Gr[j])*(P*N[i]*w[n]);
				ke.sub(3*i, 3*j, mat3da(Kab));
			}
	}
}

//-----------------------------------------------------------------------------
//! calculates the equivalent nodal forces due to hydrostatic pressure

void FEPressureLoad::PressureForce(FESurfaceElement& el, vector<double>& fe)
{
	// nr integration points
	int nint = el.GaussPoints();

	// nr of element nodes
	int neln = el.Nodes();

	// nodal coordinates
	FEMesh& mesh = *m_psurf->GetMesh();
	vec3d rt[FEElement::MAX_NODES];
	for (int j = 0; j < neln; ++j)
	{
		rt[j] = mesh.Node(el.m_node[j]).m_rt;
	}
	if (m_bshellb) {
		for (int j = 0; j<neln; ++j) {
			FENode& nd = mesh.Node(el.m_node[j]);
			rt[j] -= nd.m_d0 + nd.get_vec3d(m_dofX, m_dofY, m_dofZ) - nd.get_vec3d(m_dofSX, m_dofSY, m_dofSZ);
		}
	}

	// repeat over integration points
	zero(fe);
	double* w  = el.GaussWeights();
	for (int n=0; n<nint; ++n)
	{
		FEMaterialPoint& pt = *el.GetMaterialPoint(n);

		double* N  = el.H(n);
		double* Gr = el.Gr(n);
		double* Gs = el.Gs(n);

		// traction at integration points
		vec3d dxr(0,0,0), dxs(0,0,0);
		for (int i=0; i<neln; ++i)
		{
			dxr += rt[i]*Gr[i];
			dxs += rt[i]*Gs[i];
		}

		// evaluate pressure at this material point
		double P = -m_pressure.eval(pt);

        if (m_bshellb) P = -P;

		// force vector
		vec3d f = (dxr ^ dxs)*P*w[n];

		for (int i=0; i<neln; ++i)
		{
			fe[3*i  ] += N[i]*f.x;
			fe[3*i+1] += N[i]*f.y;
			fe[3*i+2] += N[i]*f.z;
		}
	}
}

//-----------------------------------------------------------------------------
//! calculates the equivalent nodal forces due to hydrostatic pressure

void FEPressureLoad::LinearPressureForce(FESurfaceElement& el, vector<double>& fe)
{
	// nr integration points
	int nint = el.GaussPoints();

	// nr of element nodes
	int neln = el.Nodes();

	// nodal coordinates
	FEMesh& mesh = *m_psurf->GetMesh();
	vec3d r0[FEElement::MAX_NODES];
	for (int i=0; i<neln; ++i) r0[i] = mesh.Node(el.m_node[i]).m_r0;
    if (m_bshellb) {
        for (int j=0; j<neln; ++j) {
            FENode& nd = mesh.Node(el.m_node[j]);
            r0[j] -= nd.m_d0;
        }
    }

	// repeat over integration points
	zero(fe);
	double* w  = el.GaussWeights();
	for (int n=0; n<nint; ++n)
	{
		FEMaterialPoint& pt = *el.GetMaterialPoint(n);

		double* N  = el.H(n);
		double* Gr = el.Gr(n);
		double* Gs = el.Gs(n);

		// traction at integration points
		vec3d dxr(0,0,0), dxs(0,0,0);
		for (int i=0; i<neln; ++i)
		{
			dxr += r0[i]*Gr[i];
			dxs += r0[i]*Gs[i];
		}

		// evaluate pressure at this material point
		double P = -m_pressure.eval(pt);

        if (m_bshellb) P = -P;

		vec3d f = (dxr ^ dxs)*P*w[n];

		for (int i=0; i<neln; ++i)
		{
			fe[3*i  ] += N[i]*f.x;
			fe[3*i+1] += N[i]*f.y;
			fe[3*i+2] += N[i]*f.z;
		}
	}
}

//-----------------------------------------------------------------------------

void FEPressureLoad::Serialize(DumpStream& ar)
{
	FESurfaceLoad::Serialize(ar);
}

//-----------------------------------------------------------------------------
void FEPressureLoad::UnpackLM(FEElement& el, vector<int>& lm)
{
	FEMesh& mesh = *GetSurface().GetMesh();
	int N = el.Nodes();
	lm.resize(N*3);
    if (!m_bshellb) {
        for (int i=0; i<N; ++i)
        {
            int n = el.m_node[i];
            FENode& node = mesh.Node(n);
            vector<int>& id = node.m_ID;
            
            lm[3*i  ] = id[m_dofX];
            lm[3*i+1] = id[m_dofY];
            lm[3*i+2] = id[m_dofZ];
        }
    }
    else {
        for (int i=0; i<N; ++i)
        {
            int n = el.m_node[i];
            FENode& node = mesh.Node(n);
            vector<int>& id = node.m_ID;
            
            lm[3*i  ] = id[m_dofSX];
            lm[3*i+1] = id[m_dofSY];
            lm[3*i+2] = id[m_dofSZ];
        }
    }
}

//-----------------------------------------------------------------------------
void FEPressureLoad::StiffnessMatrix(const FETimeInfo& tp, FESolver* psolver)
{
	// We only need the stiffness for nonlinear pressure forces
	if (m_blinear) return;

	// unless of course we don't want it at all
	if (m_bstiff == false) return;

	matrix ke;
	vector<int> lm;

	FESurface& surf = GetSurface();
	int npr = surf.Elements();
	for (int m=0; m<npr; ++m)
	{
		// get the surface element
		FESurfaceElement& el = m_psurf->Element(m);

		// calculate nodal normal tractions
		int neln = el.Nodes();

		// get the element stiffness matrix
		int ndof = 3*neln;
		ke.resize(ndof, ndof);

		// calculate pressure stiffness
		PressureStiffness(el, ke);

		// get the element's LM vector
		UnpackLM(el, lm);

		// assemble element matrix in global stiffness matrix
		psolver->AssembleStiffness(el.m_node, lm, ke);
	}
}

//-----------------------------------------------------------------------------
void FEPressureLoad::Residual(const FETimeInfo& tp, FEGlobalVector& R)
{
	vector<double> fe;
	vector<int> lm;

	FESurface& surf = GetSurface();
	int npr = surf.Elements();
	for (int i=0; i<npr; ++i)
	{
		FESurfaceElement& el = m_psurf->Element(i);

		// calculate nodal normal tractions
		int neln = el.Nodes();
	
		int ndof = 3*neln;
		fe.resize(ndof);

		if (m_blinear) LinearPressureForce(el, fe); else PressureForce(el, fe);

		// get the element's LM vector
		UnpackLM(el, lm);

		// add element force vector to global force vector
		R.Assemble(el.m_node, lm, fe);
	}
}
