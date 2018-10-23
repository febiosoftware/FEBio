#include "stdafx.h"
#include "FESurfaceTraction.h"
#include <FECore/FEModel.h>

BEGIN_FECORE_CLASS(FESurfaceTraction, FESurfaceLoad)
	ADD_PARAMETER(m_bshellb, "shell_bottom");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FESurfaceTraction::FESurfaceTraction(FEModel* fem) : FESurfaceLoad(fem)
{
	m_bshellb = false;
	m_blinear = false;

	// get the degrees of freedom
	m_dofX  = fem->GetDOFIndex("x");
	m_dofY  = fem->GetDOFIndex("y");
	m_dofZ  = fem->GetDOFIndex("z");
	m_dofSX = fem->GetDOFIndex("sx");
	m_dofSY = fem->GetDOFIndex("sy");
	m_dofSZ = fem->GetDOFIndex("sz");
}

//-----------------------------------------------------------------------------
void FESurfaceTraction::UnpackLM(FEElement& el, vector<int>& lm)
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
void FESurfaceTraction::GetNodalCoordinates(FESurfaceElement& el, vec3d* rt)
{
	FEMesh& mesh = *m_psurf->GetMesh();
	int neln = el.Nodes();
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
}

//-----------------------------------------------------------------------------
void FESurfaceTraction::GetReferenceNodalCoordinates(FESurfaceElement& el, vec3d* r0)
{
	FEMesh& mesh = *m_psurf->GetMesh();
	int neln = el.Nodes();
	for (int j = 0; j < neln; ++j)
	{
		r0[j] = mesh.Node(el.m_node[j]).m_r0;
	}
	if (m_bshellb) {
		for (int j = 0; j<neln; ++j) {
			FENode& nd = mesh.Node(el.m_node[j]);
			r0[j] -= nd.m_d0;
		}
	}
}

//-----------------------------------------------------------------------------
void FESurfaceTraction::Residual(const FETimeInfo& tp, FEGlobalVector& R)
{
	vector<double> fe;
	vector<int> lm;

	FESurface& surf = GetSurface();
	int npr = surf.Elements();
	for (int i = 0; i<npr; ++i)
	{
		FESurfaceElement& el = m_psurf->Element(i);

		// allocate element force vector
		int neln = el.Nodes();
		int ndof = 3 * neln;
		fe.resize(ndof);

		// Calculate the element residual
		ElementResidual(el, fe);
		
		// get the element's LM vector
		UnpackLM(el, lm);

		// add element force vector to global force vector
		R.Assemble(el.m_node, lm, fe);
	}
}

//-----------------------------------------------------------------------------
void FESurfaceTraction::ElementResidual(FESurfaceElement& el, std::vector<double>& fe)
{
	// nr integration points
	int nint = el.GaussPoints();

	// nr of element nodes
	int neln = el.Nodes();

	// nodal coordinates
	FEMesh& mesh = *m_psurf->GetMesh();
	vec3d re[FEElement::MAX_NODES];
	if (m_blinear)
		GetReferenceNodalCoordinates(el, re);
	else
		GetNodalCoordinates(el, re);

	// repeat over integration points
	zero(fe);
	double* w = el.GaussWeights();
	for (int n = 0; n<nint; ++n)
	{
		FESurfaceMaterialPoint& pt = static_cast<FESurfaceMaterialPoint&>(*el.GetMaterialPoint(n));

		double* N = el.H(n);
		double* Gr = el.Gr(n);
		double* Gs = el.Gs(n);

		// traction at integration points
		pt.dxr = vec3d(0, 0, 0);
		pt.dxs = vec3d(0, 0, 0);
		for (int i = 0; i<neln; ++i)
		{
			pt.dxr += re[i] * Gr[i];
			pt.dxs += re[i] * Gs[i];
		}
		double J = (pt.dxr ^ pt.dxs).norm();

		// evaluate traction at this material point
		vec3d t = Traction(pt);

		// flip sign if we are applying this to the bottom of a shell
		if (m_bshellb) t = -t;

		// force vector
		double wn = w[n] * J;
		for (int i = 0; i<neln; ++i)
		{
			fe[3*i    ] += N[i] * t.x * wn;
			fe[3*i + 1] += N[i] * t.y * wn;
			fe[3*i + 2] += N[i] * t.z * wn;
		}
	}
}

//-----------------------------------------------------------------------------
void FESurfaceTraction::StiffnessMatrix(const FETimeInfo& tp, FESolver* psolver)
{
	// Don't calculate stiffness for a linear load
	if (m_blinear) return;

	matrix ke;
	vector<int> lm;

	FESurface& surf = GetSurface();
	int npr = surf.Elements();
	for (int m = 0; m<npr; ++m)
	{
		// get the surface element
		FESurfaceElement& el = m_psurf->Element(m);

		// calculate nodal normal tractions
		int neln = el.Nodes();

		// get the element stiffness matrix
		int ndof = 3 * neln;
		ke.resize(ndof, ndof);

		// calculate element stiffness
		ElementStiffness(el, ke);

		// get the element's LM vector
		UnpackLM(el, lm);

		// assemble element matrix in global stiffness matrix
		psolver->AssembleStiffness(el.m_node, lm, ke);
	}
}
