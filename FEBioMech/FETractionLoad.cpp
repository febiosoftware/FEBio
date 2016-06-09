#include "stdafx.h"
#include "FETractionLoad.h"
#include "FECore/FEModel.h"

//=============================================================================
BEGIN_PARAMETER_LIST(FETractionLoad, FESurfaceLoad)
	ADD_PARAMETER(m_scale, FE_PARAM_DOUBLE    , "scale"   );
	ADD_PARAMETER(m_TC   , FE_PARAM_DATA_ARRAY, "traction");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! constructor
FETractionLoad::FETractionLoad(FEModel* pfem) : FESurfaceLoad(pfem), m_TC(FE_VEC3D)
{
	m_scale = 1.0;

	// get the degrees of freedom
	m_dofX = pfem->GetDOFIndex("x");
	m_dofY = pfem->GetDOFIndex("y");
	m_dofZ = pfem->GetDOFIndex("z");
}

//-----------------------------------------------------------------------------
//! allocate storage
void FETractionLoad::SetSurface(FESurface* ps)
{
	FESurfaceLoad::SetSurface(ps);
	m_TC.Create(ps); 
}

//-----------------------------------------------------------------------------
//! Calculate the residual for the traction load
void FETractionLoad::Residual(const FETimeInfo& tp, FEGlobalVector& R)
{
	vector<double> fe;
	vector<int> lm;

	vec3d r0[FEElement::MAX_NODES];

	FESurface& surf = *m_psurf;
	FEMesh& mesh = *surf.GetMesh();
	int NF = surf.Elements();
	for (int iel=0; iel<NF; ++iel)
	{
		FESurfaceElement& el = surf.Element(iel);

		int ndof = 3*el.Nodes();
		fe.resize(ndof);

		// nr integration points
		int nint = el.GaussPoints();

		// nr of element nodes
		int neln = el.Nodes();

		// nodal coordinates
		for (int i=0; i<neln; ++i) r0[i] = mesh.Node(el.m_node[i]).m_r0;

		// calculate the scaled traction for this element
		vec3d t = m_TC.get<vec3d>(iel)*m_scale;

		double* Gr, *Gs;
		double* N;
		double* w  = el.GaussWeights();

		// repeat over integration points
		zero(fe);
		for (int n=0; n<nint; ++n)
		{
			N  = el.H(n);
			Gr = el.Gr(n);
			Gs = el.Gs(n);


			// calculate the tangent vectors
			vec3d dxr(0,0,0), dxs(0,0,0);
			for (int i=0; i<neln; ++i) 
			{
				dxr.x += Gr[i]*r0[i].x;
				dxr.y += Gr[i]*r0[i].y;
				dxr.z += Gr[i]*r0[i].z;

				dxs.x += Gs[i]*r0[i].x;
				dxs.y += Gs[i]*r0[i].y;
				dxs.z += Gs[i]*r0[i].z;
			}

			vec3d f = t*((dxr ^ dxs).norm()*w[n]);

			for (int i=0; i<neln; ++i)
			{
				fe[3*i  ] += N[i]*f.x;
				fe[3*i+1] += N[i]*f.y;
				fe[3*i+2] += N[i]*f.z;
			}
		}

		// get the element's LM vector
		UnpackLM(el, lm);

		// add element force vector to global force vector
		R.Assemble(el.m_node, lm, fe);
	}
}

//-----------------------------------------------------------------------------
void FETractionLoad::UnpackLM(FEElement& el, vector<int>& lm)
{
	FEMesh& mesh = GetFEModel()->GetMesh();
	int N = el.Nodes();
	lm.resize(N*3);
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
