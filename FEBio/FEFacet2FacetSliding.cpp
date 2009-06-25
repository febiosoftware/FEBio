#include "stdafx.h"
#include "FEFacet2FacetSliding.h"
#include "fem.h"

//-----------------------------------------------------------------------------
// FEFacetSlidingSurface
//-----------------------------------------------------------------------------

void FEFacetSlidingSurface::Init()
{
	// initialize surface data first
	FESurface::Init();

	// count how many integration points we have
	int nint = 0, i;
	for (i=0; i<Elements(); ++i)
	{
		FESurfaceElement& el = Element(i);
		nint += el.GaussPoints();
	}

	// allocate data structures
	m_gap.create(nint);
	m_nu.create(nint);
	m_rs.create(nint);
	m_Lm.create(nint);
	m_pme.create(nint);

	// set intial values
	m_gap.zero();
	m_nu.zero();
	m_pme.set(0);
	m_Lm.zero();
}

//-----------------------------------------------------------------------------
//! Finds the (master) element that contains the projection of a (slave) node

FEElement* FEFacetSlidingSurface::FindMasterSegment(vec3d& x, vec3d& q, vec2d& r, bool& binit_nq, double tol)
{
	// get the mesh
	FEMesh& mesh = *m_pmesh;

	// see if we need to initialize the NQ structure
	if (binit_nq) m_NQ.Init();
	binit_nq = false;

	// let's find the closest master node
	int mn = m_NQ.Find(x);

	// mn is a local index, so get the global node number too
	int m = node[mn];

	// get the nodal position
	vec3d r0 = mesh.Node(m).m_rt;

	// now that we found the closest master node, lets see if we can find 
	// the best master element
	int N;

	// loop over all master elements that contain the node mn
	int nval = m_NEL.Valence(mn);
	FEElement** pe = m_NEL.ElementList(mn);
	for (int j=0; j<nval; ++j)
	{
		// get the master element
		FESurfaceElement& el = dynamic_cast<FESurfaceElement&> (*pe[j]);
		N = el.Nodes();

		// project the node on the element
		r[0] = 0;
		r[1] = 0;
		q = ProjectToSurface(el, x, r[0], r[1]);
		if (IsInsideElement(el, r[0], r[1], tol)) return pe[j];
	}

	// we did not find a master surface
	return 0;
}

//-----------------------------------------------------------------------------
// FEFacet2FacetSliding
//-----------------------------------------------------------------------------

FEFacet2FacetSliding::FEFacet2FacetSliding(FEM* pfem) : FEContactInterface(pfem), m_ss(&pfem->m_mesh), m_ms(&pfem->m_mesh)
{
	m_ntype = FE_FACET2FACET_SLIDING;
	static int ncount = 1;
	m_nID = ncount++;

	// default parameters
	m_epsn = 1.0;
}

//-----------------------------------------------------------------------------
//! Initialization routine
void FEFacet2FacetSliding::Init()
{
	// initialize surface data
	m_ss.Init();
	m_ms.Init();

	// project slave surface onto master surface
	ProjectSurface(m_ss, m_ms);
}

//-----------------------------------------------------------------------------
//! In this function we project the integration points to the master surface,
//! calculate the projection's natural coordinates and normal vector
//
void FEFacet2FacetSliding::ProjectSurface(FEFacetSlidingSurface &ss, FEFacetSlidingSurface &ms)
{
	bool bfirst = true;

	// keep a running counter of which integration 
	// point we are currently investigating
	int ni = 0;

	// get the mesh
	FEMesh* pm = ss.GetMesh();

	// loop over all slave elements
	for (int i=0; i<ss.Elements(); ++i)
	{
		// get the slave element
		FESurfaceElement& se = ss.Element(i);
		pm->UnpackElement(se);
		vec3d* re = se.rt();
		int nn = se.Nodes();

		// loop over all its integration points
		int nint = se.GaussPoints();
		for (int j=0; j<nint; ++j, ++ni)
		{
			// calculate the global coordinates of this integration point
			double* H = se.H(j);

			vec3d x(0,0,0), q;
			for (int k=0; k<nn; ++k) x += re[k]*H[k];

			// find the master segment this element belongs to
			ss.m_rs[ni] = vec2d(0,0);
			FESurfaceElement* pme = dynamic_cast<FESurfaceElement*>(ms.FindMasterSegment(x, q, ss.m_rs[ni], bfirst, 0.01));
			ss.m_pme[ni] = pme;

			if (pme)
			{
				double r = ss.m_rs[ni][0];
				double s = ss.m_rs[ni][1];

				FESurfaceElement& mel = *pme;

				// the slave normal is set to the master element normal
				ss.m_nu[ni] = ms.SurfaceNormal(mel, r, s);

				// calculate gap
				ss.m_gap[ni] = -ss.m_nu[ni]*(x - q);
			}
			else
			{
				ss.m_gap[ni] = 0;
				ss.m_Lm[ni] = 0;
			}
		}
	}
}

//-----------------------------------------------------------------------------
void FEFacet2FacetSliding::Update()
{
	// project slave surface to master surface
	ProjectSurface(m_ss, m_ms);
}

//-----------------------------------------------------------------------------
void FEFacet2FacetSliding::ShallowCopy(FEContactInterface &ci)
{
}

//-----------------------------------------------------------------------------
void FEFacet2FacetSliding::ContactForces(vector<double>& F)
{
	int i, j, k;
	vector<int> sLM, mLM, LM, en;
	vector<double> fe;

	// keep a running counter of integration points
	int ni = 0;

	// get the mesh
	FEMesh* pm = m_ss.GetMesh();

	// get the solver
	FESolver* psolver = m_pfem->m_pStep->m_psolver;

	double detJ[4], w[4], *Hs, Hm[4];

	// loop over all slave elements
	for (i=0; i<m_ss.Elements(); ++i)
	{
		FESurfaceElement& se = m_ss.Element(i);
		pm->UnpackElement(se);
		int nseln = se.Nodes();
		int nint = se.GaussPoints();

		// copy the LM vector
		sLM = se.LM();

		// nodal coordinates
		vec3d* r0 = se.r0();

		// we calculate all the metrics we need before we
		// calculate the nodal forces
		for (j=0; j<nint; ++j)
		{
			double* Gr = se.Gr(j);
			double* Gs = se.Gs(j);

			// calculate jacobian
			// note that we are integrating over the reference surface
			vec3d dxr, dxs;
			for (k=0; k<nseln; ++k)
			{
				dxr.x += Gr[k]*r0[k].x;
				dxr.y += Gr[k]*r0[k].y;
				dxr.z += Gr[k]*r0[k].z;

				dxs.x += Gs[k]*r0[k].x;
				dxs.y += Gs[k]*r0[k].y;
				dxs.z += Gs[k]*r0[k].z;
			}

			// jacobians
			detJ[j] = (dxr ^ dxs).norm();

			// integration weights
			w[j] = se.GaussWeights()[j];
		}

		// loop over all integration points
		for (int j=0; j<nint; ++j, ++ni)
		{
			// get the master element
			FESurfaceElement* pme = m_ss.m_pme[ni];
			if (pme)
			{
				FESurfaceElement& me = *pme;
				pm->UnpackElement(me);

				int nmeln = me.Nodes();
				mLM = me.LM();

				// calculate degrees of freedom
				int ndof = 3*(nseln + nmeln);

				// build the LM vector
				LM.create(ndof);
				for (k=0; k<nseln; ++k)
				{
					LM[3*k  ] = sLM[3*k  ];
					LM[3*k+1] = sLM[3*k+1];
					LM[3*k+2] = sLM[3*k+2];
				}

				for (k=0; k<nmeln; ++k)
				{
					LM[3*(k+nseln)  ] = mLM[3*k  ];
					LM[3*(k+nseln)+1] = mLM[3*k+1];
					LM[3*(k+nseln)+2] = mLM[3*k+2];
				}

				// build the en vector
				en.create(nseln+nmeln);
				for (k=0; k<nseln; ++k) en[k] = se.m_node[k];
				for (k=0; k<nmeln; ++k) en[k+nseln] = me.m_node[k];

				// calculate shape functions
				Hs = se.H(j);

				double r = m_ss.m_rs[ni][0];
				double s = m_ss.m_rs[ni][1];
				if (me.Nodes() == 4)
				{
					Hm[0] = 0.25*(1-r)*(1-s);
					Hm[1] = 0.25*(1+r)*(1-s);
					Hm[2] = 0.25*(1+r)*(1+s);
					Hm[3] = 0.25*(1-r)*(1+s);
				}

				// get normal vector
				vec3d nu = m_ss.m_nu[ni];

				// gap function
				double g = m_ss.m_gap[ni];
				
				// contact traction
				double tn = m_epsn*MBRACKET(g);

				// calculate the force vector
				fe.create(ndof);

				for (k=0; k<nseln; ++k)
				{
					fe[3*k  ] = Hs[k]*nu.x;
					fe[3*k+1] = Hs[k]*nu.y;
					fe[3*k+2] = Hs[k]*nu.z;
				}

				for (k=0; k<nmeln; ++k)
				{
					fe[3*(k+nseln)  ] = -Hm[k]*nu.x;
					fe[3*(k+nseln)+1] = -Hm[k]*nu.y;
					fe[3*(k+nseln)+2] = -Hm[k]*nu.z;
				}

				for (k=0; k<ndof; ++k) fe[k] *= tn*detJ[j]*w[j];

				// assemble the global residual
				psolver->AssembleResidual(en, LM, fe, F);
			}
		}
	}
}

//-----------------------------------------------------------------------------
void FEFacet2FacetSliding::ContactStiffness()
{
	int i, j, k, l;
	vector<int> sLM, mLM, LM, en;
	vector<double> N;
	matrix ke;

	// keep a running counter of integration points
	int ni = 0;

	// get the mesh
	FEMesh* pm = m_ss.GetMesh();

	// get the solver
	FESolver* psolver = m_pfem->m_pStep->m_psolver;

	double detJ[4], w[4], *Hs, Hm[4];

	// loop over all slave elements
	for (i=0; i<m_ss.Elements(); ++i)
	{
		FESurfaceElement& se = m_ss.Element(i);
		pm->UnpackElement(se);
		int nseln = se.Nodes();
		int nint = se.GaussPoints();

		// copy the LM vector
		sLM = se.LM();

		// nodal coordinates
		vec3d* r0 = se.r0();

		// we calculate all the metrics we need before we
		// calculate the nodal forces
		for (j=0; j<nint; ++j)
		{
			double* Gr = se.Gr(j);
			double* Gs = se.Gs(j);

			// calculate jacobian
			// note that we are integrating over the reference surface
			vec3d dxr, dxs;
			for (k=0; k<nseln; ++k)
			{
				dxr.x += Gr[k]*r0[k].x;
				dxr.y += Gr[k]*r0[k].y;
				dxr.z += Gr[k]*r0[k].z;

				dxs.x += Gs[k]*r0[k].x;
				dxs.y += Gs[k]*r0[k].y;
				dxs.z += Gs[k]*r0[k].z;
			}

			// jacobians
			detJ[j] = (dxr ^ dxs).norm();

			// integration weights
			w[j] = se.GaussWeights()[j];
		}

		// loop over all integration points
		for (int j=0; j<nint; ++j, ++ni)
		{
			// get the master element
			FESurfaceElement* pme = m_ss.m_pme[ni];
			if (pme)
			{
				FESurfaceElement& me = *pme;
				pm->UnpackElement(me);

				int nmeln = me.Nodes();
				mLM = me.LM();

				// calculate degrees of freedom
				int ndof = 3*(nseln + nmeln);

				// build the LM vector
				LM.create(ndof);
				for (k=0; k<nseln; ++k)
				{
					LM[3*k  ] = sLM[3*k  ];
					LM[3*k+1] = sLM[3*k+1];
					LM[3*k+2] = sLM[3*k+2];
				}

				for (k=0; k<nmeln; ++k)
				{
					LM[3*(k+nseln)  ] = mLM[3*k  ];
					LM[3*(k+nseln)+1] = mLM[3*k+1];
					LM[3*(k+nseln)+2] = mLM[3*k+2];
				}

				// build the en vector
				en.create(nseln+nmeln);
				for (k=0; k<nseln; ++k) en[k] = se.m_node[k];
				for (k=0; k<nmeln; ++k) en[k+nseln] = me.m_node[k];

				// calculate shape functions
				Hs = se.H(j);
				double r = m_ss.m_rs[ni][0];
				double s = m_ss.m_rs[ni][1];
				if (me.Nodes() == 4)
				{
					Hm[0] = 0.25*(1-r)*(1-s);
					Hm[1] = 0.25*(1+r)*(1-s);
					Hm[2] = 0.25*(1+r)*(1+s);
					Hm[3] = 0.25*(1-r)*(1+s);
				}
				else
				{
					Hm[0] = 1 - r - s;
					Hm[1] = r;
					Hm[2] = s;
				}

				// get normal vector
				vec3d nu = m_ss.m_nu[ni];

				// gap function
				double g = m_ss.m_gap[ni];
				
				// contact traction
				double tn = m_epsn*HEAVYSIDE(m_epsn*g);

				// calculate the N-vector
				N.create(ndof);

				for (k=0; k<nseln; ++k)
				{
					N[3*k  ] = Hs[k]*nu.x;
					N[3*k+1] = Hs[k]*nu.y;
					N[3*k+2] = Hs[k]*nu.z;
				}

				for (k=0; k<nmeln; ++k)
				{
					N[3*(k+nseln)  ] = -Hm[k]*nu.x;
					N[3*(k+nseln)+1] = -Hm[k]*nu.y;
					N[3*(k+nseln)+2] = -Hm[k]*nu.z;
				}

				// calculate the stiffness matrix
				ke.Create(ndof, ndof);
				for (k=0; k<ndof; ++k)
					for (l=0; l<ndof; ++l) ke[k][l] = tn*N[k]*N[l]*detJ[j]*w[j];

				// assemble the global residual
				psolver->AssembleStiffness(en, LM, ke);
			}
		}
	}
}

//-----------------------------------------------------------------------------
bool FEFacet2FacetSliding::Augment(int naug)
{
	return true;
}

//-----------------------------------------------------------------------------
void FEFacet2FacetSliding::Serialize(Archive &ar)
{

}
