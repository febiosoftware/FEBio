#include "stdafx.h"
#include "FESlidingInterface2.h"
#include "fem.h"

//-----------------------------------------------------------------------------
// FEContactSurface2
//-----------------------------------------------------------------------------
void FEContactSurface2::Init()
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
void FEContactSurface2::ShallowCopy(FEContactSurface2 &s)
{
	m_Lm  = s.m_Lm;
	m_gap = s.m_gap;
	m_pme.zero();
}

//-----------------------------------------------------------------------------

bool FEContactSurface2::IntersectTri(vec3d* y, vec3d r, vec3d n, double rs[2], double& g, double eps)
{
	vec3d e[2], E[2];
	mat2d G;
	e[0] = y[1]-y[0];
	e[1] = y[2]-y[0];
	vec3d m = e[0]^e[1]; m.unit();

	G[0][0] = e[0]*e[0]; G[0][1] = e[0]*e[1];
	G[1][0] = e[1]*e[0]; G[1][1] = e[1]*e[1];

	mat2d Gi = G.inverse();
	E[0] = e[0]*Gi[0][0] + e[1]*Gi[0][1];
	E[1] = e[0]*Gi[1][0] + e[1]*Gi[1][1];

	double d = n*m;
	if (d != 0)
	{
		g = m*(y[0] - r)/d;
		vec3d q = r + n*g;
		rs[0] = E[0]*(q - y[0]);
		rs[1] = E[1]*(q - y[0]);

		if ((rs[0] >= -eps) && (rs[1] >= -eps) && (rs[0]+rs[1] <= 1+eps)) return true;
	}

	return false;
}

//-----------------------------------------------------------------------------

bool FEContactSurface2::IntersectQuad(vec3d* y, vec3d r, vec3d n, double rs[2], double& g, double eps)
{
	vec3d x1[3], x2[3];
	x1[0] = y[1]; x2[0] = y[2];
	x1[1] = y[2]; x2[1] = y[3];
	x1[2] = y[0]; x2[2] = y[0];

	bool b = false;
	double rp, sp;

	if (IntersectTri(x1, r, n, rs, g, eps))
	{
		// we've intersected the first triangle
		b = true;
		rp =  1.0 - 2.0*rs[1];
		sp = -1.0 + 2.0*rs[0];
	}
	else if (IntersectTri(x2, r, n, rs, g, eps))
	{
		// we've intersected the sectond triangle
		b = true;
		rp = 1.0 - 2.0*rs[0];
		sp = 1.0 - 2.0*rs[1];
	}

	if (b)
	{
		// now that we've intersected one of the triangles
		// we try to find a more accurate intersection
		mat3d A;
		vec3d dx;
		vec3d F, F1, F2, F3;
		double H[4], H1[4], H2[4];

		double l1 = rp;
		double l2 = sp;
		double l3 = g;
		
		int nn = 0;
		int maxn = 5;
		do
		{
			H[0] = 0.25*(1 - l1)*(1 - l2);
			H[1] = 0.25*(1 + l1)*(1 - l2);
			H[2] = 0.25*(1 + l1)*(1 + l2);
			H[3] = 0.25*(1 - l1)*(1 + l2);

			H1[0] = -0.25*(1 - l2); H2[0] = -0.25*(1 - l1);
			H1[1] =  0.25*(1 - l2); H2[1] = -0.25*(1 - l1);
			H1[2] =  0.25*(1 + l2); H2[2] =  0.25*(1 - l1);
			H1[3] = -0.25*(1 + l2); H2[3] =  0.25*(1 - l1);

			F = r + n*l3 - y[0]*H[0] - y[1]*H[1] - y[2]*H[2] - y[3]*H[3];

			F1 = - y[0]*H1[0] - y[1]*H1[1] - y[2]*H1[2] - y[3]*H1[3];
			F2 = - y[0]*H2[0] - y[1]*H2[1] - y[2]*H2[2] - y[3]*H2[3];
			F3 = n;

			A[0][0] = F1.x; A[0][1] = F2.x; A[0][2] = F3.x;
			A[1][0] = F1.y; A[1][1] = F2.y; A[1][2] = F3.y;
			A[2][0] = F1.z; A[2][1] = F2.z; A[2][2] = F3.z;

			dx = -(A.inverse()*F);

			l1 += dx.x;
			l2 += dx.y;
			l3 += dx.z;

			++nn;
		}
		while ((dx.norm() > 0.001) && (nn < maxn));

		rs[0] = l1;
		rs[1] = l2;
		g     = l3;

		if ((rs[0] >= -1-eps) && (rs[0] <= 1+eps) && 
			(rs[1] >= -1-eps) && (rs[1] <= 1+eps)) return true;
	}

	return false;
}

//-----------------------------------------------------------------------------

bool FEContactSurface2::Intersect(FESurfaceElement& el, vec3d r, vec3d n, double rs[2], double& g, double eps)
{
	int N = el.Nodes();

	FEMesh& mesh = *m_pmesh;
	vec3d y[4];
	for (int i=0; i<N; ++i) y[i] = mesh.Node(el.m_node[i]).m_rt;

	if (N == 3) return IntersectTri(y, r, n, rs, g, eps);
	else if (N == 4) return IntersectQuad(y, r, n, rs, g, eps);
	return false;
}

//-----------------------------------------------------------------------------

FESurfaceElement* FEContactSurface2::FindIntersection(vec3d r, vec3d n, double rs[2])
{
	double g, gmin = 1e99;
	const double eps = 0.01;
	FESurfaceElement* pme = 0;
	for (int i=0; i<Elements(); ++i)
	{
		FESurfaceElement& el = Element(i);
		if (Intersect(el, r, n, rs, g, eps))
		{
			if (g < gmin)
			{
				pme = &el;
				gmin = g;
			}
		}	
	}

	return pme;
}

//-----------------------------------------------------------------------------
// FESlidingInterface2
//-----------------------------------------------------------------------------

FESlidingInterface2::FESlidingInterface2(FEM* pfem) : FEContactInterface(pfem), m_ss(&pfem->m_mesh), m_ms(&pfem->m_mesh)
{
	m_ntype = FE_CONTACT_SLIDING2;
	static int count = 1;
	m_nID = count++;

	// initial values
	m_knmult = 1;
	m_atol = 0.02;
	m_eps = 1;
	m_npass = 1;
}

//-----------------------------------------------------------------------------
void FESlidingInterface2::Init()
{
	// initialize surface data
	m_ss.Init();
	m_ms.Init();

	// project slave surface onto master surface
	ProjectSurface(m_ss, m_ms);

	if (m_npass == 2) ProjectSurface(m_ms, m_ss);
}

//-----------------------------------------------------------------------------
void FESlidingInterface2::ProjectSurface(FEContactSurface2& ss, FEContactSurface2& ms)
{
	FEMesh& mesh = m_pfem->m_mesh;
	FESurfaceElement* pme;
	vec3d r, nu;
	double rs[2];

	// loop over all integration points
	int n = 0;
	for (int i=0; i<ss.Elements(); ++i)
	{
		FESurfaceElement& el = ss.Element(i);
		mesh.UnpackElement(el);

		int nint = el.GaussPoints();
		for (int j=0; j<nint; ++j, ++n)
		{
			// calculate the global position of the integration point
			r = ss.Local2Global(el, j);

			// calculate the normal at this integration point
			nu = ss.SurfaceNormal(el, j);

			// find the intersection point with the master surface
			pme = ms.FindIntersection(r, nu, rs);

			ss.m_pme[n] = pme;
			ss.m_nu[n] = nu;
			ss.m_rs[n][0] = rs[0];
			ss.m_rs[n][1] = rs[1];
			if (pme)
			{
				// the node could potentially be in contact
				// find the global location of the intersection point
				vec3d q = ms.Local2Global(*pme, rs[0], rs[1]);

				// calculate the gap function
				// NOTE: this has the opposite sign compared
				// to Gerard's notes.
				ss.m_gap[n] = nu*(r - q);
			}
			else
			{
				// the node is not in contact
				ss.m_Lm[n] = 0;
				ss.m_gap[n] = 0;
			}
		}
	}
}

//-----------------------------------------------------------------------------
void FESlidingInterface2::Update()
{	
	ProjectSurface(m_ss, m_ms);
	if (m_npass == 2) ProjectSurface(m_ms, m_ss);
}

//-----------------------------------------------------------------------------
void FESlidingInterface2::ShallowCopy(FEContactInterface &ci)
{
	FESlidingInterface2& si = dynamic_cast<FESlidingInterface2&>(ci);
	m_ss.ShallowCopy(si.m_ss);
	m_ms.ShallowCopy(si.m_ms);
}

//-----------------------------------------------------------------------------
void FESlidingInterface2::ContactForces(vector<double> &F)
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

	for (int np=0; np<m_npass; ++np)
	{
		FEContactSurface2& ss = (np == 0? m_ss : m_ms);
		FEContactSurface2& ms = (np == 0? m_ms : m_ss);

		// loop over all slave elements
		ni = 0;
		for (i=0; i<ss.Elements(); ++i)
		{
			FESurfaceElement& se = ss.Element(i);
			pm->UnpackElement(se);
			int nseln = se.Nodes();
			int nint = se.GaussPoints();

			// copy the LM vector
			sLM = se.LM();

			// nodal coordinates
			vec3d* rt = se.rt();

			// we calculate all the metrics we need before we
			// calculate the nodal forces
			for (j=0; j<nint; ++j)
			{
				double* Gr = se.Gr(j);
				double* Gs = se.Gs(j);

				// calculate jacobian
				// note that we are integrating over the current surface
				vec3d dxr, dxs;
				for (k=0; k<nseln; ++k)
				{
					dxr.x += Gr[k]*rt[k].x;
					dxr.y += Gr[k]*rt[k].y;
					dxr.z += Gr[k]*rt[k].z;

					dxs.x += Gs[k]*rt[k].x;
					dxs.y += Gs[k]*rt[k].y;
					dxs.z += Gs[k]*rt[k].z;
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
				FESurfaceElement* pme = ss.m_pme[ni];
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

					double r = ss.m_rs[ni][0];
					double s = ss.m_rs[ni][1];
					if (me.Nodes() == 4)
					{
						Hm[0] = 0.25*(1-r)*(1-s);
						Hm[1] = 0.25*(1+r)*(1-s);
						Hm[2] = 0.25*(1+r)*(1+s);
						Hm[3] = 0.25*(1-r)*(1+s);
					}
					else if (me.Nodes() == 3)
					{
						Hm[0] = 1-r-s;
						Hm[1] = r;
						Hm[2] = s;
					}

					// get normal vector
					vec3d nu = ss.m_nu[ni];

					// gap function
					double g = ss.m_gap[ni];
					
					// lagrange multiplier
					double Lm = ss.m_Lm[ni];

					// contact traction
					double tn = Lm + m_eps*g;
					tn = MBRACKET(tn);

					// calculate the force vector
					fe.create(ndof);

					for (k=0; k<nseln; ++k)
					{
						fe[3*k  ] = -Hs[k]*nu.x;
						fe[3*k+1] = -Hs[k]*nu.y;
						fe[3*k+2] = -Hs[k]*nu.z;
					}

					for (k=0; k<nmeln; ++k)
					{
						fe[3*(k+nseln)  ] = Hm[k]*nu.x;
						fe[3*(k+nseln)+1] = Hm[k]*nu.y;
						fe[3*(k+nseln)+2] = Hm[k]*nu.z;
					}

					for (k=0; k<ndof; ++k) fe[k] *= tn*detJ[j]*w[j];

					// assemble the global residual
					psolver->AssembleResidual(en, LM, fe, F);
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------
void FESlidingInterface2::ContactStiffness()
{
	int i, j, k, l;
	vector<int> sLM, mLM, LM, en;
	double N[24], T1[24], T2[24], N1[24] = {0}, N2[24] = {0}, D1[24], D2[24], Nb1[24], Nb2[24];
	matrix ke;

	// keep a running counter of integration points
	int ni = 0;

	// get the mesh
	FEMesh* pm = m_ss.GetMesh();

	// get the solver
	FESolver* psolver = m_pfem->m_pStep->m_psolver;

	// see how many reformations we've had to do so far
	int nref = psolver->m_nref;

	// set higher order stiffness mutliplier
	double knmult = m_knmult;
	if (m_knmult < 0)
	{
		int ni = int(-m_knmult);
		if (nref >= ni)
		{
			knmult = 1; 
			m_pfem->m_log.printf("Higher order stiffness terms included.\n");
		}
		else knmult = 0;
	}

	double detJ[4], w[4], *Hs, Hm[4], Hmr[4], Hms[4];

	for (int np=0; np < m_npass; ++np)
	{
		FEContactSurface2& ss = (np == 0? m_ss : m_ms);
		FEContactSurface2& ms = (np == 0? m_ms : m_ss);

		// loop over all slave elements
		ni = 0;
		for (i=0; i<ss.Elements(); ++i)
		{
			FESurfaceElement& se = ss.Element(i);
			pm->UnpackElement(se);
			int nseln = se.Nodes();
			int nint = se.GaussPoints();

			// copy the LM vector
			sLM = se.LM();

			// nodal coordinates
			vec3d* rt = se.rt();

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
					dxr.x += Gr[k]*rt[k].x;
					dxr.y += Gr[k]*rt[k].y;
					dxr.z += Gr[k]*rt[k].z;

					dxs.x += Gs[k]*rt[k].x;
					dxs.y += Gs[k]*rt[k].y;
					dxs.z += Gs[k]*rt[k].z;
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
				FESurfaceElement* pme = ss.m_pme[ni];
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
					double r = ss.m_rs[ni][0];
					double s = ss.m_rs[ni][1];
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
					vec3d nu = ss.m_nu[ni];

					// gap function
					double g = ss.m_gap[ni];
					
					// lagrange multiplier
					double Lm = ss.m_Lm[ni];

					// contact traction
					double tn = Lm + m_eps*g;
					tn = MBRACKET(tn);

					double dtn = m_eps*HEAVYSIDE(Lm + m_eps*g);

					// calculate the N-vector
					for (k=0; k<nseln; ++k)
					{
						N[3*k  ] = -Hs[k]*nu.x;
						N[3*k+1] = -Hs[k]*nu.y;
						N[3*k+2] = -Hs[k]*nu.z;
					}

					for (k=0; k<nmeln; ++k)
					{
						N[3*(k+nseln)  ] = Hm[k]*nu.x;
						N[3*(k+nseln)+1] = Hm[k]*nu.y;
						N[3*(k+nseln)+2] = Hm[k]*nu.z;
					}

					// --- N O R M A L   S T I F F N E S S ---

					// create the stiffness matrix
					ke.Create(ndof, ndof);

					// add the first order term
					for (k=0; k<ndof; ++k)
						for (l=0; l<ndof; ++l) ke[k][l] = dtn*N[k]*N[l]*detJ[j]*w[j];

					// add the higher order terms (= tn*D(dg) )
					if (knmult > 0)
					{
					}

					// assemble the global residual
					psolver->AssembleStiffness(en, LM, ke);
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------
bool FESlidingInterface2::Augment(int naug)
{
	return true;
}

//-----------------------------------------------------------------------------
void FESlidingInterface2::Serialize(Archive &ar)
{

}
