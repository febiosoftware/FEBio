#include "stdafx.h"
#include "FETractionLoad.h"
#include "FECore/FEModel.h"

//-----------------------------------------------------------------------------
FETractionLoad::LOAD::LOAD()
{
	s[0] = s[1] = s[2] = s[3] = s[4] = s[5] = s[6] = s[7] = s[8] = vec3d(0,0,0);
}

//=============================================================================
BEGIN_PARAMETER_LIST(FETractionLoad, FESurfaceLoad)
	ADD_PARAMETER(m_scale   , FE_PARAM_DOUBLE, "scale"   );
	ADD_PARAMETER(m_traction, FE_PARAM_VEC3D , "traction");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! constructor
FETractionLoad::FETractionLoad(FEModel* pfem) : FESurfaceLoad(pfem)
{
	m_traction = vec3d(0,0,0);
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

	int n = ps->Elements();
	m_TC.resize(n); 

	// TODO: This assumes the traction vector was read in before the surface
	for (int i=0; i<n; i++)
	{
		for (int j=0; j<9; ++j) m_TC[i].s[j] = m_traction;
	}
}

//-----------------------------------------------------------------------------
//! Calculate the residual for the traction load
void FETractionLoad::Residual(FEGlobalVector& R)
{
	FEModel& fem = R.GetFEModel();

	vector<double> fe;
	vector<int> lm;

	vec3d r0[FEElement::MAX_NODES];

	int i, n;
	int npr = m_TC.size();
	for (int iel=0; iel<npr; ++iel)
	{
		LOAD& pc = m_TC[iel];
		FESurfaceElement& el = m_psurf->Element(iel);

		double g = m_scale;

		int ndof = 3*el.Nodes();
		fe.resize(ndof);

		// nr integration points
		int nint = el.GaussPoints();

		// nr of element nodes
		int neln = el.Nodes();

		// nodal coordinates
		for (i=0; i<neln; ++i) r0[i] = m_psurf->GetMesh()->Node(el.m_node[i]).m_r0;

		double* Gr, *Gs;
		double* N;
		double* w  = el.GaussWeights();

		vec3d dxr, dxs;

		// repeat over integration points
		zero(fe);
		for (n=0; n<nint; ++n)
		{
			N  = el.H(n);
			Gr = el.Gr(n);
			Gs = el.Gs(n);

			// calculate the traction at the integration point
			vec3d t = el.eval(pc.s, n)*g;

			// calculate the tangent vectors
			dxr = dxs = vec3d(0,0,0);
			for (i=0; i<neln; ++i) 
			{
				dxr.x += Gr[i]*r0[i].x;
				dxr.y += Gr[i]*r0[i].y;
				dxr.z += Gr[i]*r0[i].z;

				dxs.x += Gs[i]*r0[i].x;
				dxs.y += Gs[i]*r0[i].y;
				dxs.z += Gs[i]*r0[i].z;
			}

			vec3d f = t*((dxr ^ dxs).norm()*w[n]);

			for (i=0; i<neln; ++i)
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
void FETractionLoad::Serialize(DumpStream& ar)
{
	FESurfaceLoad::Serialize(ar);

	if (ar.IsSaving())
	{
		ar << (int) m_TC.size();
		for (int i=0; i < (int) m_TC.size(); ++i)
		{
			LOAD& d = m_TC[i];
			ar << d.s[0] << d.s[1] << d.s[2] << d.s[3];
			ar << d.s[4] << d.s[5] << d.s[6] << d.s[7] << d.s[8];
		}
	}
	else
	{
		int n;
		ar >> n;
		m_TC.resize(n);
		for (int i=0; i<n; ++i)
		{
			LOAD& d = m_TC[i];
			ar >> d.s[0] >> d.s[1] >> d.s[2] >> d.s[3];
			ar >> d.s[4] >> d.s[5] >> d.s[6] >> d.s[7] >> d.s[8];
		}
	}
}

//-----------------------------------------------------------------------------
// \deprecated This is only needed for parsing the obsolete 1.2 format
bool FETractionLoad::SetFacetAttribute(int nface, const char* szatt, const char* szval)
{
	LOAD& tc = TractionLoad(nface);
	if      (strcmp(szatt, "id") == 0) {}
//	else if (strcmp(szatt, "lc") == 0) tc.lc = atoi(szval) - 1;
	else if (strcmp(szatt, "tx") == 0)
	{
		double tx = atof(szval);
		tc.s[0].x = tc.s[1].x = tc.s[2].x = tc.s[3].x = tx;
		tc.s[4].x = tc.s[5].x = tc.s[6].x = tc.s[7].x = tx;
		tc.s[8].x = tx;
	}
	else if (strcmp(szatt, "ty") == 0)
	{
		double ty = atof(szval);
		tc.s[0].y = tc.s[1].y = tc.s[2].y = tc.s[3].y = ty;
		tc.s[4].y = tc.s[5].y = tc.s[6].y = tc.s[7].y = ty;
		tc.s[8].y = ty;
	}
	else if (strcmp(szatt, "tz") == 0)
	{
		double tz = atof(szval);
		tc.s[0].z = tc.s[1].z = tc.s[2].z = tc.s[3].z = tz;
		tc.s[4].z = tc.s[5].z = tc.s[6].z = tc.s[7].z = tz;
		tc.s[8].z = tz;
	}
	else return false;

	return true;
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
