#include "stdafx.h"
#include "FESlidingInterface2.h"
#include "fem.h"

//-----------------------------------------------------------------------------
// FEContactSurface2
//-----------------------------------------------------------------------------

FEContactSurface2::FEContactSurface2(FEM* pfem) : FESurface(&pfem->m_mesh)
{ 
	m_pfem = pfem; 
}

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

	// allocate biphasic stuff
	if (m_pfem->m_pStep->m_itype == FE_STATIC_PORO)
	{
		m_pg.create(nint);
		m_pg.zero();
	}
}

//-----------------------------------------------------------------------------
void FEContactSurface2::ShallowCopy(FEContactSurface2 &s)
{
	m_Lm  = s.m_Lm;
	m_gap = s.m_gap;
	m_pme.zero();
}

//-----------------------------------------------------------------------------
// This function calculates the intersection of a ray with a triangle
// and returns true if the ray intersects the triangle.
//
bool FEContactSurface2::IntersectTri(vec3d* y, vec3d r, vec3d n, double rs[2], double& g, double eps)
{
	vec3d e[2], E[2];
	mat2d G;

	// create base vectors on triangle
	e[0] = y[1]-y[0];
	e[1] = y[2]-y[0];

	// create triangle normal
	vec3d m = e[0]^e[1]; m.unit();

	double d = n*m;
	if (d != 0)
	{
		// distance from r to triangle
		g = m*(y[0] - r)/d;

		// intersection point with triangle
		vec3d q = r + n*g;

		// next, we decompose q into its components
		// in the triangle basis
		// we need to create the dual basis
		// first, we calculate the metric tensor
		G[0][0] = e[0]*e[0]; G[0][1] = e[0]*e[1];
		G[1][0] = e[1]*e[0]; G[1][1] = e[1]*e[1];

		// and its inverse
		mat2d Gi = G.inverse();

		// build dual basis
		E[0] = e[0]*Gi[0][0] + e[1]*Gi[0][1];
		E[1] = e[0]*Gi[1][0] + e[1]*Gi[1][1];

		// get the components
		rs[0] = E[0]*(q - y[0]);
		rs[1] = E[1]*(q - y[0]);

		// see if the intersection point is inside the triangle
		if ((rs[0] >= -eps) && (rs[1] >= -eps) && (rs[0]+rs[1] <= 1+eps)) return true;
	}

	return false;
}

//-----------------------------------------------------------------------------
//! This function calculates the intersection of a ray with a quad
//! and returns true if the ray intersected.
//!
bool FEContactSurface2::IntersectQuad(vec3d* y, vec3d r, vec3d n, double rs[2], double& g, double eps)
{
	// first we're going to see if the ray intersects the two subtriangles
	vec3d x1[3], x2[3];
	x1[0] = y[0]; x2[0] = y[2];
	x1[1] = y[1]; x2[1] = y[3];
	x1[2] = y[3]; x2[2] = y[1];

	bool b = false;
	double rp, sp;

	if (IntersectTri(x1, r, n, rs, g, eps))
	{
		// we've intersected the first triangle
		b = true;
		rp = -1.0 + 2.0*rs[0];
		sp = -1.0 + 2.0*rs[1];
	}
	else if (IntersectTri(x2, r, n, rs, g, eps))
	{
		// we've intersected the second triangle
		b = true;
		rp = 1.0 - 2.0*rs[0];
		sp = 1.0 - 2.0*rs[1];
	}

	// if one of the triangels was intersected,
	// we calculate a more accurate projection
	if (b)
	{
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
			// shape functions of quad
			H[0] = 0.25*(1 - l1)*(1 - l2);
			H[1] = 0.25*(1 + l1)*(1 - l2);
			H[2] = 0.25*(1 + l1)*(1 + l2);
			H[3] = 0.25*(1 - l1)*(1 + l2);

			// shape function derivatives
			H1[0] = -0.25*(1 - l2); H2[0] = -0.25*(1 - l1);
			H1[1] =  0.25*(1 - l2); H2[1] = -0.25*(1 + l1);
			H1[2] =  0.25*(1 + l2); H2[2] =  0.25*(1 + l1);
			H1[3] = -0.25*(1 + l2); H2[3] =  0.25*(1 - l1);

			// calculate residual
			F = r + n*l3 - y[0]*H[0] - y[1]*H[1] - y[2]*H[2] - y[3]*H[3];

			// residual derivatives
			F1 = - y[0]*H1[0] - y[1]*H1[1] - y[2]*H1[2] - y[3]*H1[3];
			F2 = - y[0]*H2[0] - y[1]*H2[1] - y[2]*H2[2] - y[3]*H2[3];
			F3 = n;

			// set up the tangent matrix
			A[0][0] = F1.x; A[0][1] = F2.x; A[0][2] = F3.x;
			A[1][0] = F1.y; A[1][1] = F2.y; A[1][2] = F3.y;
			A[2][0] = F1.z; A[2][1] = F2.z; A[2][2] = F3.z;

			// calculate solution increment
			dx = -(A.inverse()*F);

			// update solution
			l1 += dx.x;
			l2 += dx.y;
			l3 += dx.z;

			++nn;
		}
		while ((dx.norm() > 1e-7) && (nn < maxn));

		// store results
		rs[0] = l1;
		rs[1] = l2;
		g     = l3;

		// see if the point is inside the quad
		if ((rs[0] >= -1-eps) && (rs[0] <= 1+eps) && 
			(rs[1] >= -1-eps) && (rs[1] <= 1+eps)) return true;
	}

	return false;
}

//-----------------------------------------------------------------------------
//! This function calculates the intersection of a ray with a surface element.
//! It simply calls the tri or quad intersection function based on the type
//! of element.
//!
bool FEContactSurface2::Intersect(FESurfaceElement& el, vec3d r, vec3d n, double rs[2], double& g, double eps)
{
	int N = el.Nodes();

	// get the element nodes
	FEMesh& mesh = *m_pmesh;
	vec3d y[4];
	for (int i=0; i<N; ++i) y[i] = mesh.Node(el.m_node[i]).m_rt;

	// call the correct intersection function
	if (N == 3) return IntersectTri(y, r, n, rs, g, eps);
	else if (N == 4) return IntersectQuad(y, r, n, rs, g, eps);

	// if we get here, the ray did not intersect the element
	return false;
}

//-----------------------------------------------------------------------------
//! This function finds the element which is intersected by the ray (r,n).
//! It returns a pointer to the element, as well as the isoparametric coordinates
//! of the intersection point.
//!
FESurfaceElement* FEContactSurface2::FindIntersection(vec3d r, vec3d n, double rs[2], double eps)
{
	double g, gmin = 1e99, r2[2] = {rs[0], rs[1]};
	FESurfaceElement* pme = 0;

	// loop over all surface element
	for (int i=0; i<Elements(); ++i)
	{
		FESurfaceElement& el = Element(i);

		// see if the ray intersects this element
		if (Intersect(el, r, n, r2, g, eps))
		{
			// see if this is the best intersection found so far
			// TODO: should I put a limit on how small g can
			//       be to be considered a valid intersection?
			if (g < gmin)
			{
				// keep results
				pme = &el;
				gmin = g;
				rs[0] = r2[0];
				rs[1] = r2[1];
			}
		}	
	}

	// return the intersected element (or zero if none)
	return pme;
}

//-----------------------------------------------------------------------------
// FESlidingInterface2
//-----------------------------------------------------------------------------

FESlidingInterface2::FESlidingInterface2(FEM* pfem) : FEContactInterface(pfem), m_ss(pfem), m_ms(pfem)
{
	m_ntype = FE_CONTACT_SLIDING2;
	static int count = 1;
	m_nID = count++;

	// initial values
	m_knmult = 1;
	m_atol = 0.02;
	m_eps = 1;
	m_epsp = 1;
	m_npass = 1;
	m_stol = 0.01;
	m_bsymm = true;
	m_srad = 0.01;
}

//-----------------------------------------------------------------------------
void FESlidingInterface2::Init()
{
	// initialize surface data
	m_ss.Init();
	m_ms.Init();

	bool bporo = m_pfem->m_pStep->m_itype == FE_STATIC_PORO;

	// for now we enforce penalty method for biphasic contact
	if (bporo) m_blaugon = 0;

	// this contact implementation requires a non-symmetric stiffness matrix
	// so inform the FEM class
	if (!m_bsymm) m_pfem->SetSymmetryFlag(false);

	// make sure we are using full-Newton
	if (bporo && (m_pfem->m_pStep->m_psolver->m_maxups != 0))
	{
		m_pfem->m_pStep->m_psolver->m_maxups = 0;
		m_pfem->m_log.printbox("WARNING", "The biphasic contact algorithm does not work with BFGS yet.\nThe full-Newton method will be used instead.");
	}

	// update sliding interface data
	Update();
}

//-----------------------------------------------------------------------------
void FESlidingInterface2::ProjectSurface(FEContactSurface2& ss, FEContactSurface2& ms)
{
	FEMesh& mesh = m_pfem->m_mesh;
	FESurfaceElement* pme;
	vec3d r, nu;
	double rs[2];

	bool bporo = (m_pfem->m_pStep->m_itype == FE_STATIC_PORO);

	double ps[4], p1;

	double R = m_srad*mesh.GetBoundingBox().radius();

	// loop over all integration points
	int n = 0;
	for (int i=0; i<ss.Elements(); ++i)
	{
		FESurfaceElement& el = ss.Element(i);
		mesh.UnpackElement(el);

		int ne = el.Nodes();
		int nint = el.GaussPoints();

		// get the nodal pressures
		if (bporo)
		{
			for (int j=0; j<ne; ++j) ps[j] = el.pt()[j];
		}

		for (int j=0; j<nint; ++j, ++n)
		{
			// calculate the global position of the integration point
			r = ss.Local2Global(el, j);

			// get the pressure at the integration point
			if (bporo) p1 = el.eval(ps, j);

			// calculate the normal at this integration point
			nu = ss.SurfaceNormal(el, j);

			// first see if the old intersected face is still good enough
			pme = ss.m_pme[n];
			if (pme)
			{
				double g;

				// see if the ray intersects this element
				if (ms.Intersect(*pme, r, nu, rs, g, m_stol))
				{
					ss.m_rs[n][0] = rs[0];
					ss.m_rs[n][1] = rs[1];
				}
				else pme = 0;
			}

			// find the intersection point with the master surface
			if (pme == 0) pme = ms.FindIntersection(r, nu, rs, m_stol);

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
				double g = nu*(r - q);

				if (g < R)
				{
					ss.m_gap[n] = g;

					// calculate the pressure gap function
					if (bporo && ss.m_gap[n] >= 0)
					{
						mesh.UnpackElement(*pme);
						double p2 = pme->eval(pme->pt(), rs[0], rs[1]);
						ss.m_pg[n] = p1 - p2;
					}
					else ss.m_pg[n] = 0;
				}
				else
				{
					ss.m_Lm[n] = 0;
					ss.m_gap[n] = 0;
					ss.m_pme[n] = 0;
				}
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
	int i, j, n, id;

	// project the surfaces onto each other
	// this will update the gap functions as well
	ProjectSurface(m_ss, m_ms);
	if (m_npass == 2) ProjectSurface(m_ms, m_ss);

	// set poro flag
	bool bporo = (m_pfem->m_pStep->m_itype == FE_STATIC_PORO);

	// only continue if we are doing a poro-elastic simulation
	if (bporo == false) return;

	// now that the nodes have been projected, we need to figure out
	// if we need to modify the constraints on the pressure dofs
	for (int np=0; np<m_npass; ++np)
	{
		FEContactSurface2& ss = (np == 0? m_ss : m_ms);

		// first, mark all nodes as non free-draining (= pos. ID)
		// this is done by setting the dof's equation number
		// to a negative number
		for (i=0; i<ss.Nodes(); ++i) 
		{
			id = ss.Node(i).m_ID[6];
			if (id < -1) ss.Node(i).m_ID[6] = -id-2;
		}

		// keep a running counter of integration points.
		int ni = 0;

		// loop over all elements
		for (n=0; n<ss.Elements(); ++n)
		{
			FESurfaceElement& el = ss.Element(n);

			// get the nodal tractions at the integration points
			double ti[4], gi[4];
			int nint = el.GaussPoints();
			int neln = el.Nodes();
			for (i=0; i<nint; ++i, ++ni) 
			{
				gi[i] = ss.m_gap[ni];
				ti[i] = m_eps*MBRACKET(gi[i]);
			}

			// setup the (over-determined) system to find the nodal values
			// TODO: this is the same for all surface elements, so
			//       maybe we should do this only once
			matrix A;
			A.Create(nint, neln);
			for (i=0; i<nint; ++i)
			{
				double* H = el.H(i);
				for (j=0; j<neln; ++j) A[i][j] = H[j];
			}

			double tn[4];

			if (nint == neln)
			{
				matrix Ai = A.inverse();

				for (i=0; i<neln; ++i)
				{
					tn[i] = 0;
					for (j=0; j<nint; ++j) tn[i] += Ai[i][j]*ti[j];
				}
			}
			else
			{
				// TODO: I still need to do this case
				assert(false);
			}

			for (i=0; i<neln; ++i)
			{
				FENode& node = ss.Node(el.m_lnode[i]);
				id = node.m_ID[6];
				if ((id >= 0) && (tn[i] <= 0))
				{
					// fix the prescribed dof
					node.m_ID[6] = -id-2;

					// set the fluid pressure to zero
					node.m_pt = 0;
				}
			}
		}
	}

	// if we only did single pass, the dofs of the secondary surface 
	// have not been modified, so we modify them here.
	if (m_npass == 1)
	{

	}
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
	double detJ[4], w[4], *Hs, Hm[4];

	// get the mesh
	FEMesh* pm = m_ss.GetMesh();

	// get the solver
	FESolver* psolver = m_pfem->m_pStep->m_psolver;

	// set poro flag
	bool bporo = (m_pfem->m_pStep->m_itype == FE_STATIC_PORO);
	
	// get the poro-elasticity symmetry flag
	bool bsymm = m_pfem->m_bsym_poro;

	// if we're using the symmetric formulation
	// we need to multiply with the timestep
	double dt = (bsymm?m_pfem->m_pStep->m_dt:1);

	// loop over the nr of passes
	for (int np=0; np<m_npass; ++np)
	{
		// get slave and master surface
		FEContactSurface2& ss = (np == 0? m_ss : m_ms);
		FEContactSurface2& ms = (np == 0? m_ms : m_ss);

		// keep a running counter of integration points
		int ni = 0;

		// loop over all slave elements
		for (i=0; i<ss.Elements(); ++i)
		{
			// get the surface element
			FESurfaceElement& se = ss.Element(i);
			pm->UnpackElement(se);

			// get the nr of nodes and integration points
			int nseln = se.Nodes();
			int nint = se.GaussPoints();

			// copy the LM vector; we'll need it later
			sLM = se.LM();

			// we calculate all the metrics we need before we
			// calculate the nodal forces
			for (j=0; j<nint; ++j)
			{
				// get the base vectors
				vec3d g[2];
				ss.CoBaseVectors(se, j, g);

				// jacobians: J = |g0xg1|
				detJ[j] = (g[0] ^ g[1]).norm();

				// integration weights
				w[j] = se.GaussWeights()[j];
			}

			// loop over all integration points
			// note that we are integrating over the current surface
			for (j=0; j<nint; ++j, ++ni)
			{
				// get the master element
				FESurfaceElement* pme = ss.m_pme[ni];
				if (pme)
				{
					// get the master element
					FESurfaceElement& me = *pme;
					pm->UnpackElement(me);

					// get the nr of master element nodes
					int nmeln = me.Nodes();

					// copy LM vector
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
					for (k=0; k<nseln; ++k) en[k      ] = se.m_node[k];
					for (k=0; k<nmeln; ++k) en[k+nseln] = me.m_node[k];

					// get slave element shape functions
					Hs = se.H(j);

					// get master element shape functions
					double r = ss.m_rs[ni][0];
					double s = ss.m_rs[ni][1];
					me.shape_fnc(Hm, r, s);

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

					// do the biphasic stuff
					// TODO: I should only do this when the node is actually in contact
					//       in other words, when g >= 0
					if (bporo && (g >= 0))
					{
						// calculate nr of pressure dofs
						int ndof = nseln + nmeln;

						// calculate the flow rate
						double wn = m_epsp*ss.m_pg[ni];

						// fill the LM
						LM.create(ndof);
						for (k=0; k<nseln; ++k) LM[k        ] = sLM[3*nseln+k];
						for (k=0; k<nmeln; ++k) LM[k + nseln] = mLM[3*nmeln+k];

						// fill the force array
						fe.create(ndof);
						for (k=0; k<nseln; ++k) fe[k      ] =  Hs[k];
						for (k=0; k<nmeln; ++k) fe[k+nseln] = -Hm[k];

						// NOTE: note that dt is either the timestep
						//       or one, depending on whether we are 
						//       using the symmetric poro version or not
						for (k=0; k<ndof; ++k) fe[k] *= dt*wn*detJ[j]*w[j];

						// assemble residual
						psolver->AssembleResidual(en, LM, fe, F);
					}
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
	double detJ[4], w[4], *Hs, Hm[4], pt[4], dpr[4], dps[4];
	double N[24];
	matrix ke;

	// set the poro flag
	bool bporo = (m_pfem->m_pStep->m_itype == FE_STATIC_PORO);

	// get the poroelasticity symmetry flag
	bool bsymm = m_pfem->m_bsym_poro;

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

	// do single- or two-pass
	for (int np=0; np < m_npass; ++np)
	{
		// get the slave and master surface
		FEContactSurface2& ss = (np == 0? m_ss : m_ms);
		FEContactSurface2& ms = (np == 0? m_ms : m_ss);

		// keep a running counter of the integration points
		int ni = 0;

		// loop over all slave elements
		for (i=0; i<ss.Elements(); ++i)
		{
			// get ths slave element
			FESurfaceElement& se = ss.Element(i);
			pm->UnpackElement(se);

			// get nr of nodes and integration points
			int nseln = se.Nodes();
			int nint = se.GaussPoints();

			// copy the LM vector
			sLM = se.LM();

			// we calculate all the metrics we need before we
			// calculate the nodal forces
			for (j=0; j<nint; ++j)
			{
				// get the base vectors
				vec3d g[2];
				ss.CoBaseVectors(se, j, g);

				// jacobians: J = |g0xg1|
				detJ[j] = (g[0] ^ g[1]).norm();

				// integration weights
				w[j] = se.GaussWeights()[j];

				// pressure
				if (bporo)
				{
					pt[j] = se.eval(se.pt(), j);
					dpr[j] = se.eval_deriv1(se.pt(), j);
					dps[j] = se.eval_deriv2(se.pt(), j);
				}
			}

			// loop over all integration points
			for (j=0; j<nint; ++j, ++ni)
			{
				// get the master element
				FESurfaceElement* pme = ss.m_pme[ni];
				if (pme)
				{
					// --- S O L I D - S O L I D   C O N T A C T ---

					FESurfaceElement& me = *pme;
					pm->UnpackElement(me);

					// get the nr of master nodes
					int nmeln = me.Nodes();

					// copy the LM vector
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
					for (k=0; k<nseln; ++k) en[k      ] = se.m_node[k];
					for (k=0; k<nmeln; ++k) en[k+nseln] = me.m_node[k];

					// slave shape functions
					Hs = se.H(j);

					// master shape functions
					double r = ss.m_rs[ni][0];
					double s = ss.m_rs[ni][1];
					me.shape_fnc(Hm, r, s);

					// get slave normal vector
					vec3d nu = ss.m_nu[ni];

					// gap function
					double g = ss.m_gap[ni];
					
					// lagrange multiplier
					double Lm = ss.m_Lm[ni];

					// contact traction
					double tn = Lm + m_eps*g;
					tn = MBRACKET(tn);

					double dtn = m_eps*HEAVYSIDE(Lm + m_eps*g);

					// a. NxN-term
					//------------------------------------

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

					// create the stiffness matrix
					ke.Create(ndof, ndof);

					for (k=0; k<ndof; ++k)
						for (l=0; l<ndof; ++l) ke[k][l] = dtn*N[k]*N[l]*detJ[j]*w[j];

					// b. A-term
					//-------------------------------------

					for (k=0; k<nseln; ++k) N[k      ] =  Hs[k];
					for (k=0; k<nmeln; ++k) N[k+nseln] = -Hm[k];

					double* Gr = se.Gr(j);
					double* Gs = se.Gs(j);
					vec3d gs[2];
					ss.CoBaseVectors(se, j, gs);

					mat3d S1, S2;
					S1.skew(gs[0]);
					S2.skew(gs[1]);

					if (!m_bsymm)
					{
						for (k=0; k<nseln+nmeln; ++k)
							for (l=0; l<nseln; ++l)
							{
								ke[k*3  ][l*3  ] += tn*w[j]*N[k]*(-Gr[l]*S2[0][0] + Gs[l]*S1[0][0]);
								ke[k*3  ][l*3+1] += tn*w[j]*N[k]*(-Gr[l]*S2[0][1] + Gs[l]*S1[0][1]);
								ke[k*3  ][l*3+2] += tn*w[j]*N[k]*(-Gr[l]*S2[0][2] + Gs[l]*S1[0][2]);

								ke[k*3+1][l*3  ] += tn*w[j]*N[k]*(-Gr[l]*S2[1][0] + Gs[l]*S1[1][0]);
								ke[k*3+1][l*3+1] += tn*w[j]*N[k]*(-Gr[l]*S2[1][1] + Gs[l]*S1[1][1]);
								ke[k*3+1][l*3+2] += tn*w[j]*N[k]*(-Gr[l]*S2[1][2] + Gs[l]*S1[1][2]);

								ke[k*3+2][l*3  ] += tn*w[j]*N[k]*(-Gr[l]*S2[2][0] + Gs[l]*S1[2][0]);
								ke[k*3+2][l*3+1] += tn*w[j]*N[k]*(-Gr[l]*S2[2][1] + Gs[l]*S1[2][1]);
								ke[k*3+2][l*3+2] += tn*w[j]*N[k]*(-Gr[l]*S2[2][2] + Gs[l]*S1[2][2]);
							}
					}
					// c. M-term
					//---------------------------------------

					vec3d gm[2];
					ms.CoBaseVectors(me, r, s, gm);

					mat2d A, Ai;
					A[0][0] = gs[0]*gm[0]; A[0][1] = gs[0]*gm[1];
					A[1][0] = gs[1]*gm[0]; A[1][1] = gs[1]*gm[1];

					Ai = A.inverse();

					vec3d Gm[2];
					Gm[0] = gs[0]*Ai[0][0] + gs[1]*Ai[0][1];
					Gm[1] = gs[0]*Ai[1][0] + gs[1]*Ai[1][1];

					double Hmr[4], Hms[4];
					me.shape_deriv(Hmr, Hms, r, s);

					if (!m_bsymm)
					{
						for (k=0; k<nseln; ++k)
							for (l=0; l<nseln+nmeln; ++l)
							{
								ke[(k+nseln)*3  ][l*3  ] -= tn*detJ[j]*w[j]*(Hmr[k]*nu.x*Gm[0].x + Hms[k]*nu.x*Gm[1].x)*N[l];
								ke[(k+nseln)*3  ][l*3+1] -= tn*detJ[j]*w[j]*(Hmr[k]*nu.x*Gm[0].y + Hms[k]*nu.x*Gm[1].y)*N[l];
								ke[(k+nseln)*3  ][l*3+2] -= tn*detJ[j]*w[j]*(Hmr[k]*nu.x*Gm[0].z + Hms[k]*nu.x*Gm[1].z)*N[l];

								ke[(k+nseln)*3+1][l*3  ] -= tn*detJ[j]*w[j]*(Hmr[k]*nu.y*Gm[0].x + Hms[k]*nu.y*Gm[1].x)*N[l];
								ke[(k+nseln)*3+1][l*3+1] -= tn*detJ[j]*w[j]*(Hmr[k]*nu.y*Gm[0].y + Hms[k]*nu.y*Gm[1].y)*N[l];
								ke[(k+nseln)*3+1][l*3+2] -= tn*detJ[j]*w[j]*(Hmr[k]*nu.y*Gm[0].z + Hms[k]*nu.y*Gm[1].z)*N[l];

								ke[(k+nseln)*3+2][l*3  ] -= tn*detJ[j]*w[j]*(Hmr[k]*nu.z*Gm[0].x + Hms[k]*nu.z*Gm[1].x)*N[l];
								ke[(k+nseln)*3+2][l*3+1] -= tn*detJ[j]*w[j]*(Hmr[k]*nu.z*Gm[0].y + Hms[k]*nu.z*Gm[1].y)*N[l];
								ke[(k+nseln)*3+2][l*3+2] -= tn*detJ[j]*w[j]*(Hmr[k]*nu.z*Gm[0].z + Hms[k]*nu.z*Gm[1].z)*N[l];
							}
					}

					// assemble the global stiffness
					psolver->AssembleStiffness(en, LM, ke);

					// --- B I P H A S I C   S T I F F N E S S ---
					if (bporo && (g >= 0))
					{
						// the variable dt is either the timestep or one
						// depending on whether we are using the symmetric
						// poro version or not.
						double dt = (bsymm?m_pfem->m_pStep->m_dt:1);

						// --- S O L I D - P R E S S U R E   C O N T A C T ---

						int ndof = 4*(nseln+nmeln);
						LM.create(ndof);
						for (k=0; k<nseln; ++k)
						{
							LM[4*k  ] = sLM[3*k  ];			// x-dof
							LM[4*k+1] = sLM[3*k+1];			// y-dof
							LM[4*k+2] = sLM[3*k+2];			// z-dof
							LM[4*k+3] = sLM[3*nseln+k];		// p-dof
						}
						for (k=0; k<nmeln; ++k)
						{
							LM[4*(k+nseln)  ] = mLM[3*k  ];			// x-dof
							LM[4*(k+nseln)+1] = mLM[3*k+1];			// y-dof
							LM[4*(k+nseln)+2] = mLM[3*k+2];			// z-dof
							LM[4*(k+nseln)+3] = mLM[3*nmeln+k];		// p-dof
						}

						for (k=0; k<nseln; ++k) N[k      ] =  Hs[k];
						for (k=0; k<nmeln; ++k) N[k+nseln] = -Hm[k];

						ke.Create(ndof, ndof);
						ke.zero();

						double dpmr, dpms;
						dpmr = me.eval_deriv1(me.pt(), r, s);
						dpms = me.eval_deriv2(me.pt(), r, s);

						for (k=0; k<nseln+nmeln; ++k)
							for (l=0; l<nseln+nmeln; ++l)
							{
								ke[4*k + 3][4*l  ] += dt*w[j]*detJ[j]*m_epsp*N[k]*N[l]*(dpmr*Gm[0].x + dpms*Gm[1].x);
								ke[4*k + 3][4*l+1] += dt*w[j]*detJ[j]*m_epsp*N[k]*N[l]*(dpmr*Gm[0].y + dpms*Gm[1].y);
								ke[4*k + 3][4*l+2] += dt*w[j]*detJ[j]*m_epsp*N[k]*N[l]*(dpmr*Gm[0].z + dpms*Gm[1].z);
							}

						double H2r[8] = {0}, H2s[8] = {0};
						for (k=0; k<nmeln; ++k) 
						{
							H2r[k+nseln] = Hmr[k];
							H2s[k+nseln] = Hms[k];
						}

						double G2r[24], G2s[24];
						for (k=0; k<nseln+nmeln; ++k)
						{
							G2r[3*k  ] = N[k]*Gm[0].x;
							G2r[3*k+1] = N[k]*Gm[0].y;
							G2r[3*k+2] = N[k]*Gm[0].z;

							G2s[3*k  ] = N[k]*Gm[1].x;
							G2s[3*k+1] = N[k]*Gm[1].y;
							G2s[3*k+2] = N[k]*Gm[1].z;
						}

						double wn = m_epsp*ss.m_pg[ni];

						if (!m_bsymm)
						{
							for (k=0; k<nseln+nmeln; ++k)
								for (l=0; l<nseln+nmeln; ++l)
								{
									ke[4*k + 3][4*l  ] += dt*w[j]*detJ[j]*wn*(H2r[k]*G2r[3*l  ] + H2s[k]*G2s[3*l  ]);
									ke[4*k + 3][4*l+1] += dt*w[j]*detJ[j]*wn*(H2r[k]*G2r[3*l+1] + H2s[k]*G2s[3*l+1]);
									ke[4*k + 3][4*l+2] += dt*w[j]*detJ[j]*wn*(H2r[k]*G2r[3*l+2] + H2s[k]*G2s[3*l+2]);
								}

							for (k=0; k<nseln+nmeln; ++k)
								for (l=0; l<nseln+nmeln; ++l)
								{
									mat3d A;

									A[0][0] = -Gr[l]*S2[0][0] + Gs[l]*S1[0][0];
									A[0][1] = -Gr[l]*S2[0][1] + Gs[l]*S1[0][1];
									A[0][2] = -Gr[l]*S2[0][2] + Gs[l]*S1[0][2];

									A[1][0] = -Gr[l]*S2[1][0] + Gs[l]*S1[1][0];
									A[1][1] = -Gr[l]*S2[1][1] + Gs[l]*S1[1][1];
									A[1][2] = -Gr[l]*S2[1][2] + Gs[l]*S1[1][2];

									A[2][0] = -Gr[l]*S2[2][0] + Gs[l]*S1[2][0];
									A[2][1] = -Gr[l]*S2[2][1] + Gs[l]*S1[2][1];
									A[2][2] = -Gr[l]*S2[2][2] + Gs[l]*S1[2][2];

									ke[4*k + 3][4*l  ] += dt*w[j]*wn*N[k]*(A[0][0]*nu.x + A[0][1]*nu.y + A[0][2]*nu.z);
									ke[4*k + 3][4*l+1] += dt*w[j]*wn*N[k]*(A[1][0]*nu.x + A[1][1]*nu.y + A[1][2]*nu.z);
									ke[4*k + 3][4*l+2] += dt*w[j]*wn*N[k]*(A[2][0]*nu.x + A[2][1]*nu.y + A[2][2]*nu.z);
								}
	
							psolver->AssembleStiffness(en, LM, ke);
						}


						// --- P R E S S U R E - P R E S S U R E   C O N T A C T ---

						ndof = nseln+nmeln;

						for (k=0; k<nseln; ++k) N[k      ] =  Hs[k];
						for (k=0; k<nmeln; ++k) N[k+nseln] = -Hm[k];

						LM.create(ndof);
						for (k=0; k<nseln; ++k) LM[k      ] = sLM[3*nseln+k];
						for (k=0; k<nmeln; ++k) LM[k+nseln] = mLM[3*nmeln+k];

						// build the "element" stiffness
						ke.Create(ndof, ndof);
						for (k=0; k<ndof; ++k)
							for (l=0; l<ndof; ++l) ke[k][l] = -dt*m_epsp*w[j]*detJ[j]*N[k]*N[l];

						// assemble the global stiffness
						psolver->AssembleStiffness(en, LM, ke);
					}
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
