// FESlidingInterface.cpp: implementation of the FESlidingInterface class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FESlidingInterface.h"
#include "fem.h"
#include "FESolidSolver.h"
#include "log.h"

///////////////////////////////////////////////////////////////////////////////
// FESlidingInterface
///////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
//! constructor
FESlidingInterface::FESlidingInterface(FEM* pfem) : FEContactInterface(pfem), m_ss(&pfem->m_mesh), m_ms(&pfem->m_mesh)
{
	static int count = 1;
	m_ntype = FE_CONTACT_SLIDING;

	m_mu = 0;
	m_epsf = 0;

	m_naugmin = 0;
	m_naugmax = 10;

	m_gtol = 0;

	m_stol = 0.01;

	m_ktmult = 0;
	m_knmult = 1;

	m_breloc = true;

	m_nplc = -1;
	m_pplc = 0;
	m_nsegup = 0;	// always do segment updates
	m_bautopen = false;
	m_nID = count++;
};

//-----------------------------------------------------------------------------
//! Calculates the auto penalty factor

void FESlidingInterface::CalcAutoPenalty(FEContactSurface& s)
{
	int i, j, k, n, m, *fn, *en, *pl;

	FEMesh& mesh = *s.GetMesh();

	int LT_HEX[6][4] = {
		{ 0, 1, 5, 4 },
		{ 1, 2, 6, 5 },
		{ 2, 3, 7, 6 },
		{ 3, 0, 4, 7 },
		{ 3, 2, 1, 0 },
		{ 4, 5, 6, 7 }};

	int LT_PENTA[5][4] = {
		{ 0, 1, 4, 3 },
		{ 1, 2, 5, 4 },
		{ 2, 0, 3, 5 },
		{ 2, 1, 0, 0 },
		{ 3, 4, 5, 5 }};

	int LT_TET[4][4] = {
		{ 0, 1, 3, 3 },
		{ 1, 2, 3, 3 },
		{ 2, 0, 3, 3 },
		{ 2, 1, 0, 0 }};


	// create the node-element list.
	// That is, find for each node to which element it belongs
	FENodeElemList NE;
	NE.Create(mesh);

	FENodeElemList NF;
	NF.Create(s);

	// now, we need to find for each face, what element it belongs to
	bool bfound;
	FEElement *pe;
	bool b[4] = {true};
	for (i=0; i<s.Elements(); ++i)
	{
		FESurfaceElement& face = s.Element(i);

		// pick a node of this element
		fn = face.m_node;
		n = fn[0];

		// for this node, loop over all elements this node belongs to
		// and see if this face belongs to the element
		bfound = false;
		FEElement** pel = NE.ElementList(n);
		for (j=0; j<NE.Valence(n); ++j)
		{
			// get the element
			pe = pel[j];
			en = pe->m_node;

			// see what kind of element this is
			switch (pe->Type())
			{
			case FE_HEX:
			case FE_RIHEX:
			case FE_UDGHEX:
				{
					for (k=0; k<6; ++k)
					{
						pl = LT_HEX[k];
						b[1] = ((fn[1] == en[ pl[0] ]) || (fn[1] == en[pl[1]]) || (fn[1] == en[pl[2]]) || (fn[1] == en[pl[3]]));
						b[2] = ((fn[2] == en[ pl[0] ]) || (fn[2] == en[pl[1]]) || (fn[2] == en[pl[2]]) || (fn[2] == en[pl[3]]));
//						b[3] = ((fn[3] == en[ pl[0] ]) || (fn[3] == en[pl[1]]) || (fn[3] == en[pl[2]]) || (fn[3] == en[pl[3]]));

						if (b[1] && b[2]) { bfound = true; break; }
					}
				}
				break;
			case FE_PENTA:
				{
					for (k=0; k<5; ++k)
					{
						pl = LT_PENTA[k];
						b[1] = ((fn[1] == en[ pl[0] ]) || (fn[1] == en[pl[1]]) || (fn[1] == en[pl[2]]) || (fn[1] == en[pl[3]]));
						b[2] = ((fn[2] == en[ pl[0] ]) || (fn[2] == en[pl[1]]) || (fn[2] == en[pl[2]]) || (fn[2] == en[pl[3]]));
//						b[3] = ((fn[3] == en[ pl[0] ]) || (fn[3] == en[pl[1]]) || (fn[3] == en[pl[2]]) || (fn[3] == en[pl[3]]));

						if (b[1] && b[2]) { bfound = true; break; }
					}
				}
				break;
			case FE_TET:
				{
					for (k=0; k<4; ++k)
					{
						pl = LT_TET[k];
						b[1] = ((fn[1] == en[ pl[0] ]) || (fn[1] == en[pl[1]]) || (fn[1] == en[pl[2]]) || (fn[1] == en[pl[3]]));
						b[2] = ((fn[2] == en[ pl[0] ]) || (fn[2] == en[pl[1]]) || (fn[2] == en[pl[2]]) || (fn[2] == en[pl[3]]));
//						b[3] = ((fn[3] == en[ pl[0] ]) || (fn[3] == en[pl[1]]) || (fn[3] == en[pl[2]]) || (fn[3] == en[pl[3]]));

						if (b[1] && b[2]) { bfound = true; break; }
					}
				}
				break;
			case FE_SHELL_QUAD:
				{
					b[1] = ((fn[1] == en[0]) || (fn[1] == en[1]) || (fn[1] == en[2]) || (fn[1] == en[3]));
					b[2] = ((fn[2] == en[0]) || (fn[2] == en[1]) || (fn[2] == en[2]) || (fn[2] == en[3]));
//					b[3] = ((fn[3] == en[0]) || (fn[3] == en[1]) || (fn[3] == en[2]) || (fn[3] == en[3]));

					if (b[1] && b[2]) bfound = true;
				}
			case FE_SHELL_TRI:
				{
					b[1] = ((fn[1] == en[0]) || (fn[1] == en[1]) || (fn[1] == en[2]));
					b[2] = ((fn[2] == en[0]) || (fn[2] == en[1]) || (fn[2] == en[2]));
//					b[3] = ((fn[3] == en[0]) || (fn[3] == en[1]) || (fn[3] == en[2]));

					if (b[1] && b[2]) bfound = true;
				}
				break;
			}

			if (bfound) break;
		}

		// we found the element to which this face belongs
		if (bfound)
		{
			// get the element's material
			FESolidMaterial* pm = dynamic_cast<FESolidMaterial*>(m_pfem->GetMaterial(pe->GetMatID()));

			double K = pm->BulkModulus();

			// calculate the facet area
			double area = s.FaceArea(face);

			// calculate the volume element
			double vol = mesh.ElementVolume(*pe);

			// set the auto calculation factor
			double eps = K*area*area / vol;

			for (k=0; k<face.Nodes(); ++k)
			{
				m = face.m_lnode[k];
				s.eps[m] += eps;
			}
		}
		else assert(false);
	}

	for (i=0; i<s.Nodes(); ++i) s.eps[i] /= NF.Valence(i);
}

//-----------------------------------------------------------------------------
//! Initializes the sliding interface data

void FESlidingInterface::Init()
{
	// create the surfaces
	m_ss.Init();
	m_ms.Init();

	// update the master surface normals
//	m_ms.UpdateNormals();

	// project slave surface onto master surface
	ProjectSurface(m_ss, m_ms, m_breloc);

	if (m_bautopen)	CalcAutoPenalty(m_ss);

	// set penalty load curve
	if (m_nplc >= 0) m_pplc = m_pfem->GetLoadCurve(m_nplc);

	// for two-pass algorithms we repeat the previous
	// two steps with master and slave switched
	if (m_npass == 2)
	{
//		m_ss.UpdateNormals();
		ProjectSurface(m_ms, m_ss, true);

		if (m_bautopen) CalcAutoPenalty(m_ms);
	}
}

//-----------------------------------------------------------------------------
//!  Projects the slave surface onto the master surface.
//!  That is for each slave node we determine the closest
//!  master element and the projection of the slave node onto
//!  this master element.

// TODO: this function needs to identify the different types of node-slave contact:
//   1/ first contact
//   2/ crossing of element boundary
//	 3/ contact termination 
//			either by failure to find master segment or when g < tolerance

void FESlidingInterface::ProjectSurface(FEContactSurface& ss, FEContactSurface& ms, bool bupseg, bool bmove)
{
	int i;
	double r, s;
	FEElement* pme;
	bool bfirst = true;

	// spatial position of slave node
	vec3d x;

	// slave node projection
	vec3d q;

	// get the mesh
	FEMesh& mesh = *ss.GetMesh();

	// get the logfile
	Logfile& log = GetLogfile();

	// loop over all slave nodes
	for (i=0; i<ss.Nodes(); ++i)
	{
		// get the node
		FENode& node = ss.Node(i);

		// get the nodal position
		x = node.m_rt;

		// get the global node number
		int m = ss.node[i];

		// get the previous master element (if any)
		pme = ss.pme[i];

		// If the node is in contact, let's see if the node still is 
		// on the same master element
		if (pme != 0)
		{
			FESurfaceElement& mel = dynamic_cast<FESurfaceElement&>(*pme);

			r = ss.rs[i][0];
			s = ss.rs[i][1];

			q = ms.ProjectToSurface(mel, x, r, s);
			ss.rs[i][0] = r;
			ss.rs[i][1] = s;

			// we only check when we can update the segments
			// otherwise, we just stick with this element, even
			// if the node is no longer inside it.
			if (bupseg)
			{
				if (!ms.IsInsideElement(mel, r, s, m_stol) && bupseg)
				{
					// see if the node might have moved to another master element
					FEElement* pold = pme; 
					ss.rs[i] = vec2d(0,0);
					pme = ms.FindMasterSegment(x, q, ss.rs[i], bfirst, m_stol);

					if (pme == 0)
					{
						// nope, if has genuinly left contact
						int* n = pold->m_node;
//						log.printf("node %d has left element (%d, %d, %d, %d)\n", m+1, n[0]+1, n[1]+1, n[2]+1, n[3]+1);
					}
					else if (m_mu*m_epsf > 0)
					{
						// the node has moved to another master segment.
						// If friction is active we need to translate the frictional
						// data to the new master segment.
						FESurfaceElement& eo = dynamic_cast<FESurfaceElement&>(*pold);
						FESurfaceElement& en = dynamic_cast<FESurfaceElement&>(*pme );
						MapFrictionData(i, ss, ms, en, eo, q);
					}
				}
			}
		}
		else if (bupseg)
		{
			// get the master element
			// don't forget to initialize the search for the first node!
			ss.rs[i] = vec2d(0,0);
			pme = ms.FindMasterSegment(x, q, ss.rs[i], bfirst, m_stol);
			if (pme)
			{
				// the node has come into contact so make sure to initialize
				// the previous natural coordinates for friction.
				ss.rsp[i] = ss.rs[i];
			}
		}

		// if we found a master element, update the gap and normal data
		ss.pme[i] = pme;
		if (pme != 0)
		{
			FESurfaceElement& mel = dynamic_cast<FESurfaceElement&>(*pme);

			r = ss.rs[i][0];
			s = ss.rs[i][1];

			// if this is a new contact, copy the current coordinates
			// to the previous ones
			ss.M[i] = ss.Metric0(mel, r, s);

			// the slave normal is set to the master element normal
			ss.nu[i] = ss.SurfaceNormal(mel, r, s);

			// calculate gap
			ss.gap[i] = -(ss.nu[i]*(x - q)) + ss.off[i];
			if (bmove && (ss.gap[i]>0))
			{
				node.m_r0 = node.m_rt = q + ss.nu[i]*ss.off[i];
				ss.gap[i] = 0;
			}

			// TODO: what should we do if the gap function becomes
			// negative? setting the Lagrange multipliers to zero
			// might make the system unstable.
/*			if (ss.gap[i] < 0)
			{
				ss.Lm[i] = 0;
				ss.Lt[i][0] = 0;
				ss.Lt[i][1] = 0;
				ss.pme[i] = 0;
			}
*/		}
		else
		{
			// TODO: Is this a good criteria for out-of-contact?
			//		 perhaps this is not even necessary.
			// since the node is not in contact, we set the gap function 
			// and Lagrangian multiplier to zero
			ss.gap[i] = 0;
			ss.Lm[i]  = 0;
			ss.Lt[i][0] = ss.Lt[i][1] = 0;
		}
	}
}

//-----------------------------------------------------------------------------
//! updates sliding interface data

void FESlidingInterface::Update()
{
	static bool bfirst = true;

	// get the iteration number
	// we need this number to see if we can do segment updates or not
	int niter = m_pfem->m_pStep->m_psolver->m_niter;

	// should we do a segment update or not?
	// TODO: check what happens when m_nsegup == -1 and m_npass = 2;
	// We have to make sure that in this case, both surfaces get at least
	// one pass!
	bool bupdate = (bfirst || (m_nsegup == 0)? true : (niter <= m_nsegup));

	// update master surfaces
	m_ms.Update();

	// project slave surface onto master surface
	// this also calculates the nodal gap functions
	ProjectSurface(m_ss, m_ms, bupdate);

	if (m_npass == 2)
	{
		m_ss.Update();
		ProjectSurface(m_ms, m_ss, bupdate);
	}

	// set the first-entry-flag to false
	bfirst = false;
}

//-----------------------------------------------------------------------------

void FESlidingInterface::ContactForces(vector<double>& F)
{
	int j, k, l, m, n, np;
	int nseln, nmeln, ndof;

	// element contact force vector
	vector<double> fe;

	// the lm array for this force vector
	vector<int> lm;

	// the en array
	vector<int> en;

	// the elements LM vectors
	vector<int> sLM;
	vector<int> mLM;

	// get the mesh
	FEMesh& mesh = m_pfem->m_mesh;

	// get the solver
	FESolidSolver* psolver = dynamic_cast<FESolidSolver*>(m_pfem->m_pStep->m_psolver);

	vec3d* r0;
	double w[4];
	double* Gr, *Gs;
	double detJ[4];
	vec3d dxr, dxs;

	// do two-pass
	for (np=0; np<m_npass; ++np)
	{
		// pick the slave and master surfaces
		FEContactSurface& ss = (np==0? m_ss : m_ms);
		FEContactSurface& ms = (np==0? m_ms : m_ss);

		// loop over all slave facets
		int ne = ss.Elements();
		for (j=0; j<ne; ++j)
		{
			// get the slave element
			FESurfaceElement& sel = ss.Element(j);
			mesh.UnpackElement(sel);
			nseln = sel.Nodes();

			// get the LM array
			sLM = sel.LM();

			// nodal coordinates
			r0 = sel.r0();

			// we calculate all the metrics we need before we
			// calculate the nodal forces
			for (n=0; n<nseln; ++n)
			{
				Gr = sel.Gr(n);
				Gs = sel.Gs(n);

				// calculate jacobian
				// note that we are integrating over the reference surface
				dxr = dxs = vec3d(0,0,0);
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
				detJ[n] = (dxr ^ dxs).norm();

				// integration weights
				w[n] = sel.GaussWeights()[n];
			}

			// loop over slave element nodes (which are the integration points as well)
			// and calculate the contact nodal force
			for (n=0; n<nseln; ++n)
			{
				// get the local node number
				m = sel.m_lnode[n];

				// see if this node's constraint is active
				// that is, if it has a master element associated with it
				// TODO: is this a good way to test for an active constraint
				// The rigid wall criteria seems to work much better.
				if (ss.pme[m] != 0)
				{
					// This node is active and could lead to a non-zero
					// contact force.
					// get the master element
					FESurfaceElement& mel = dynamic_cast<FESurfaceElement&> (*ss.pme[m]);
					mesh.UnpackElement(mel, FE_UNPACK_LM);
					mLM = mel.LM();

					// calculate the degrees of freedom
					nmeln = mel.Nodes();
					ndof = 3*(nmeln+1);
					fe.create(ndof);

					// calculate the nodal force
					ContactNodalForce(m, ss, mel, fe);

					// multiply force with weights
					for (l=0; l<ndof; ++l) fe[l] *= detJ[n]*w[n];
					
					// fill the lm array
					lm.create(3*(nmeln+1));
					lm[0] = sLM[n*3  ];
					lm[1] = sLM[n*3+1];
					lm[2] = sLM[n*3+2];

					for (l=0; l<nmeln; ++l)
					{
						lm[3*(l+1)  ] = mLM[l*3  ];
						lm[3*(l+1)+1] = mLM[l*3+1];
						lm[3*(l+1)+2] = mLM[l*3+2];
					}

					// fill the en array
					en.create(nmeln+1);
					en[0] = sel.m_node[n];
					for (l=0; l<nmeln; ++l) en[l+1] = mel.m_node[l];

					// assemble into global force vector
					psolver->AssembleResidual(en, lm, fe, F);
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! Calculates the contact force on a slave node.
//! \param[in] m local node number
//! \param[out] fe force vector

void FESlidingInterface::ContactNodalForce(int m, FEContactSurface& ss, FESurfaceElement& mel, vector<double>& fe)
{
	int k, l;

	vec3d dxr, dxs;

	// normal force
	double tn, Ln;

	// gap function
	double gap;

	// tangents
	vec3d tau1, tau2;

	// master element nodes
	vec3d rtm[4];

	// shape function values
	double H[4], Hr[4], Hs[4];

	// contact vectors
	double N[15], N1[15], N2[15], T1[15], T2[15], D1[15], D2[15];
	double A[2][2], M[2][2], K[2][2];
	double detA;

	double eps, scale = Penalty();

	// get the mesh
	FEMesh& mesh = m_pfem->m_mesh;

	double Tt[2];

	int nmeln, ndof;

	// gap function
	gap = ss.gap[m];

	// normal penalty
	eps = ss.eps[m]*scale;

	// get slave node normal force
	Ln = ss.Lm[m];
	tn = Ln + eps*gap;
	tn = MBRACKET(tn);

	// get the slave node normal
	vec3d& nu = ss.nu[m];

	nmeln = mel.Nodes();
	ndof = 3*(1 + nmeln);

	// metric tensors
	mat2d Mk = ss.M[m];
	mat2d Mki = Mk.inverse();

	// get the master element node positions
	for (k=0; k<nmeln; ++k) rtm[k] = mel.rt()[k];

	// isoparametric coordinates of the projected slave node
	// onto the master element
	double r = ss.rs[m][0];
	double s = ss.rs[m][1];

	// get the coordinates at the previous step
	double rp = ss.rsp[m][0];
	double sp = ss.rsp[m][1];

	// get the master shape function values at this slave node
	if (nmeln == 4)
	{
		// quadrilateral
		H[0] = 0.25*(1-r)*(1-s);
		H[1] = 0.25*(1+r)*(1-s);
		H[2] = 0.25*(1+r)*(1+s);
		H[3] = 0.25*(1-r)*(1+s);
	}
	else if (nmeln == 3)
	{
		// triangle
		H[0] = 1 - r - s;
		H[1] = r;
		H[2] = s;
	}

	// --- N O R M A L   T R A C T I O N ---

	// calculate contact vectors for normal traction
	N[0] = nu.x;
	N[1] = nu.y;
	N[2] = nu.z;
	for (l=0; l<nmeln; ++l)
	{
		N[3*(l+1)  ] = -H[l]*nu.x;
		N[3*(l+1)+1] = -H[l]*nu.y;
		N[3*(l+1)+2] = -H[l]*nu.z;
	}

	// calculate force vector
	for (l=0; l<ndof; ++l) fe[l] = tn*N[l];

	// --- T A N G E N T I A L   T R A C T I O N ---
	if (m_mu*m_epsf > 0)
	{
		// Lagrangian traction
		double Lt[2];
		Lt[0] = ss.Lt[m][0];
		Lt[1] = ss.Lt[m][1];

		// calculate contact vector for tangential traction
		// only if both the friction coefficient and friction
		// penalty factor are non-zero

		// get the master shape function values at this slave node
		if (nmeln == 4)
		{
			Hr[0] = -0.25*(1-s); Hs[0] = -0.25*(1-r);
			Hr[1] =  0.25*(1-s); Hs[1] = -0.25*(1+r);
			Hr[2] =  0.25*(1+s); Hs[2] =  0.25*(1+r);
			Hr[3] = -0.25*(1+s); Hs[3] =  0.25*(1-r);
		}
		else if (nmeln == 3)
		{
			Hr[0] = -1; Hs[0] = -1;
			Hr[1] =  1; Hs[1] =  0;
			Hr[2] =  0; Hs[2] =  1;
		}

		// get the tangent vectors
		tau1 = tau2 = vec3d(0,0,0);
		for (k=0; k<nmeln; ++k)
		{
			tau1.x += Hr[k]*rtm[k].x;
			tau1.y += Hr[k]*rtm[k].y;
			tau1.z += Hr[k]*rtm[k].z;
		
			tau2.x += Hs[k]*rtm[k].x;
			tau2.y += Hs[k]*rtm[k].y;
			tau2.z += Hs[k]*rtm[k].z;
		}

		// set up the Ti vectors
		T1[0] = tau1.x; T2[0] = tau2.x;
		T1[1] = tau1.y; T2[1] = tau2.y;
		T1[2] = tau1.z; T2[2] = tau2.z;

		for (k=0; k<nmeln; ++k) 
		{
			T1[(k+1)*3  ] = -H[k]*tau1.x;
			T1[(k+1)*3+1] = -H[k]*tau1.y;
			T1[(k+1)*3+2] = -H[k]*tau1.z;

			T2[(k+1)*3  ] = -H[k]*tau2.x;
			T2[(k+1)*3+1] = -H[k]*tau2.y;
			T2[(k+1)*3+2] = -H[k]*tau2.z;
		}

		// set up the Ni vectors
		N1[0] = N2[0] = 0;
		N1[1] = N2[1] = 0;
		N1[2] = N2[2] = 0;

		for (k=0; k<nmeln; ++k) 
		{
			N1[(k+1)*3  ] = -Hr[k]*nu.x;
			N1[(k+1)*3+1] = -Hr[k]*nu.y;
			N1[(k+1)*3+2] = -Hr[k]*nu.z;

			N2[(k+1)*3  ] = -Hs[k]*nu.x;
			N2[(k+1)*3+1] = -Hs[k]*nu.y;
			N2[(k+1)*3+2] = -Hs[k]*nu.z;
		}

		// calculate metric tensor
		M[0][0] = tau1*tau1; M[0][1] = tau1*tau2; 
		M[1][0] = tau2*tau1; M[1][1] = tau2*tau2; 

		// calculate curvature tensor
		K[0][0] = 0; K[0][1] = 0;
		K[1][0] = 0; K[1][1] = 0;

		if (nmeln == 4)
		{
			const double h[4] = {0.25, -0.25, 0.25, -0.25};
			for (k=0; k<nmeln; ++k)
			{
				K[0][1] += (nu*rtm[k])*h[k];
				K[1][0] += (nu*rtm[k])*h[k];
			}
		}

		// setup A matrix
		A[0][0] = M[0][0] + gap*K[0][0];
		A[0][1] = M[0][1] + gap*K[0][1];
		A[1][0] = M[1][0] + gap*K[1][0];
		A[1][1] = M[1][1] + gap*K[1][1];

		detA = A[0][0]*A[1][1] - A[0][1]*A[1][0];

		// setup Di vectors
		for (k=0; k<ndof; ++k)
		{
			D1[k] = (1/detA)*(A[1][1]*(T1[k]+gap*N1[k]) - A[0][1]*(T2[k] + gap*N2[k]));
			D2[k] = (1/detA)*(A[0][0]*(T2[k]+gap*N2[k]) - A[0][1]*(T1[k] + gap*N1[k]));
		}

		// calculate friction tractions
		// a. calculate trial state
		Tt[0] = Lt[0] + m_epsf*(Mk[0][0]*(r - rp) + Mk[0][1]*(s - sp));
		Tt[1] = Lt[1] + m_epsf*(Mk[1][0]*(r - rp) + Mk[1][1]*(s - sp));

		double TMT = Tt[0]*(Mki[0][0]*Tt[0]+Mki[0][1]*Tt[1])+Tt[1]*(Mki[1][0]*Tt[0]+Mki[1][1]*Tt[1]);
		assert(TMT >= 0);

		double phi = sqrt(TMT) - m_mu*Ln;

		// b. return map
		if (phi > 0)
		{
			Tt[0] = m_mu*Ln*Tt[0]/sqrt(TMT);
			Tt[1] = m_mu*Ln*Tt[1]/sqrt(TMT);
		}

		// tangential force vector
		for (l=0; l<ndof; ++l) fe[l] -= (Tt[0]*D1[l] + Tt[1]*D2[l]);
	}
}

//-----------------------------------------------------------------------------

void FESlidingInterface::ContactStiffness()
{
	int j, k, l, n, m, np;
	int nseln, nmeln, ndof;

	matrix ke;

	vector<int> lm(15);
	vector<int> en(5);

	double *Gr, *Gs, w[4];
	vec3d *r0;

	double detJ[4];
	vec3d dxr, dxs;

	vector<int> sLM;
	vector<int> mLM;

	FEMesh& mesh = m_pfem->m_mesh;

	FESolidSolver* psolver = dynamic_cast<FESolidSolver*>(m_pfem->m_pStep->m_psolver);

	// do two-pass
	for (np=0; np<m_npass; ++np)
	{
		// get the master and slave surface
		FEContactSurface& ss = (np==0?m_ss:m_ms);	
		FEContactSurface& ms = (np==0?m_ms:m_ss);	

		// loop over all slave elements
		int ne = ss.Elements();
		for (j=0; j<ne; ++j)
		{
			// unpack the slave element
			FESurfaceElement& se = ss.Element(j);
			mesh.UnpackElement(se);
			nseln = se.Nodes();

			// get the LM array
			sLM = se.LM();

			// get the nodal coordinates
			r0 = se.r0();

			// get all the metrics we need 
			for (n=0; n<nseln; ++n)
			{
				Gr = se.Gr(n);
				Gs = se.Gs(n);

				// calculate jacobian
				dxr = dxs = vec3d(0,0,0);
				for (k=0; k<nseln; ++k)
				{
					dxr.x += Gr[k]*r0[k].x;
					dxr.y += Gr[k]*r0[k].y;
					dxr.z += Gr[k]*r0[k].z;

					dxs.x += Gs[k]*r0[k].x;
					dxs.y += Gs[k]*r0[k].y;
					dxs.z += Gs[k]*r0[k].z;
				}

				detJ[n] = (dxr ^ dxs).norm();
				w[n] = se.GaussWeights()[n];
			}

			// loop over all integration points (that is nodes)
			for (n=0; n<nseln; ++n)
			{
				m = se.m_lnode[n];

				// see if this node's constraint is active
				// that is, if it has a master element associated with it
				if (ss.pme[m] != 0)
				{
					// get the master element
					FESurfaceElement& me = dynamic_cast<FESurfaceElement&> (*ss.pme[m]);
					mesh.UnpackElement(me);

					// get the masters' LM array
					mLM = me.LM();
					nmeln = me.Nodes();
					ndof = 3*(nmeln+1);

					// calculate the stiffness matrix
					ke.Create(ndof, ndof);
					ContactNodalStiffness(m, ss, me, ke);

					// muliply with weights
					for (k=0; k<ndof; ++k)
						for (l=0; l<ndof; ++l) ke[k][l] *= detJ[n]*w[n];

					// fill the lm array
					lm[0] = sLM[n*3  ];
					lm[1] = sLM[n*3+1];
					lm[2] = sLM[n*3+2];

					for (k=0; k<nmeln; ++k)
					{
						lm[3*(k+1)  ] = mLM[k*3  ];
						lm[3*(k+1)+1] = mLM[k*3+1];
						lm[3*(k+1)+2] = mLM[k*3+2];
					}

					// create the en array
					en.create(nmeln+1);
					en[0] = se.m_node[n];
					for (k=0; k<nmeln; ++k) en[k+1] = me.m_node[k];
						
					// assemble stiffness matrix
					psolver->AssembleStiffness(en, lm, ke);
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------

void FESlidingInterface::ContactNodalStiffness(int m, FEContactSurface& ss, FESurfaceElement& mel, matrix& ke)
{
	int i, j, k, l;

	vector<int> lm(15);
	vector<int> en(5);

	vec3d dxr, dxs;
	double H[4], Hr[4], Hs[4];
	double N[15], T1[15], T2[15], N1[15], N2[15], D1[15], D2[15], Nb1[15], Nb2[15];

	// get the mesh
	FEMesh& mesh = m_pfem->m_mesh;

	// nr of element nodes and degrees of freedom 
	int nmeln = mel.Nodes();
	int ndof = 3*(1 + nmeln);

	// penalty factor
	double scale = Penalty();
	double eps = ss.eps[m]*scale;

	// get the master element node positions
	vec3d* rt = mel.rt();

	// slave node natural coordinates in master element
	double r = ss.rs[m][0];
	double s = ss.rs[m][1];

	// slave gap
	double gap = ss.gap[m];

	// lagrange multiplier
	double Lm = ss.Lm[m];

	// get slave node normal force
	double tn = Lm + eps*gap;
	tn = MBRACKET(tn);

	// get the slave node normal
	vec3d& nu = ss.nu[m];

	// get the master shape function values at this slave node
	if (nmeln == 4)
	{
		// quadrilateral
		H[0] = 0.25*(1-r)*(1-s);
		H[1] = 0.25*(1+r)*(1-s);
		H[2] = 0.25*(1+r)*(1+s);
		H[3] = 0.25*(1-r)*(1+s);

		Hr[0] = -0.25*(1-s); Hs[0] = -0.25*(1-r);
		Hr[1] =  0.25*(1-s); Hs[1] = -0.25*(1+r);
		Hr[2] =  0.25*(1+s); Hs[2] =  0.25*(1+r);
		Hr[3] = -0.25*(1+s); Hs[3] =  0.25*(1-r);
	}
	else if (nmeln == 3)
	{
		// triangle
		H[0] = 1 - r - s;
		H[1] = r;
		H[2] = s;

		Hr[0] = -1; Hs[0] = -1;
		Hr[1] =  1; Hs[1] =  0;
		Hr[2] =  0; Hs[2] =  1;
	}

	// get the tangent vectors
	vec3d tau1(0,0,0), tau2(0,0,0);
	for (k=0; k<nmeln; ++k)
	{
		tau1.x += Hr[k]*rt[k].x;
		tau1.y += Hr[k]*rt[k].y;
		tau1.z += Hr[k]*rt[k].z;
	
		tau2.x += Hs[k]*rt[k].x;
		tau2.y += Hs[k]*rt[k].y;
		tau2.z += Hs[k]*rt[k].z;
	}

	// set up the N vector
	N[0] = nu.x;
	N[1] = nu.y;
	N[2] = nu.z;

	for (k=0; k<nmeln; ++k) 
	{
		N[(k+1)*3  ] = -H[k]*nu.x;
		N[(k+1)*3+1] = -H[k]*nu.y;
		N[(k+1)*3+2] = -H[k]*nu.z;
	}
	
	// set up the Ti vectors
	T1[0] = tau1.x; T2[0] = tau2.x;
	T1[1] = tau1.y; T2[1] = tau2.y;
	T1[2] = tau1.z; T2[2] = tau2.z;

	for (k=0; k<nmeln; ++k) 
	{
		T1[(k+1)*3  ] = -H[k]*tau1.x;
		T1[(k+1)*3+1] = -H[k]*tau1.y;
		T1[(k+1)*3+2] = -H[k]*tau1.z;

		T2[(k+1)*3  ] = -H[k]*tau2.x;
		T2[(k+1)*3+1] = -H[k]*tau2.y;
		T2[(k+1)*3+2] = -H[k]*tau2.z;
	}

	// set up the Ni vectors
	N1[0] = N2[0] = 0;
	N1[1] = N2[1] = 0;
	N1[2] = N2[2] = 0;

	for (k=0; k<nmeln; ++k) 
	{
		N1[(k+1)*3  ] = -Hr[k]*nu.x;
		N1[(k+1)*3+1] = -Hr[k]*nu.y;
		N1[(k+1)*3+2] = -Hr[k]*nu.z;

		N2[(k+1)*3  ] = -Hs[k]*nu.x;
		N2[(k+1)*3+1] = -Hs[k]*nu.y;
		N2[(k+1)*3+2] = -Hs[k]*nu.z;
	}

	// calculate metric tensor
	mat2d M;
	M[0][0] = tau1*tau1; M[0][1] = tau1*tau2; 
	M[1][0] = tau2*tau1; M[1][1] = tau2*tau2; 

	// calculate reciprocal metric tensor
	mat2d Mi = M.inverse();

	// calculate curvature tensor
	double K[2][2] = {0};
	const double Grs[4] = {0.25, -0.25, 0.25, -0.25};
	if (nmeln == 4)
	{
		for (k=0; k<nmeln; ++k)
		{
			K[0][1] += (nu*rt[k])*Grs[k];
			K[1][0] += (nu*rt[k])*Grs[k];
		}
	}

	// setup A matrix A = M + gK
	double A[2][2];
	A[0][0] = M[0][0] + gap*K[0][0];
	A[0][1] = M[0][1] + gap*K[0][1];
	A[1][0] = M[1][0] + gap*K[1][0];
	A[1][1] = M[1][1] + gap*K[1][1];

	// calculate determinant of A
	double detA = A[0][0]*A[1][1] - A[0][1]*A[1][0];

	// setup Di vectors
	for (k=0; k<ndof; ++k)
	{
		D1[k] = (1/detA)*(A[1][1]*(T1[k]+gap*N1[k]) - A[0][1]*(T2[k] + gap*N2[k]));
		D2[k] = (1/detA)*(A[0][0]*(T2[k]+gap*N2[k]) - A[0][1]*(T1[k] + gap*N1[k]));
	}

	// setup Nbi vectors
	for (k=0; k<ndof; ++k)
	{
		Nb1[k] = N1[k] - K[0][1]*D2[k];
		Nb2[k] = N2[k] - K[0][1]*D1[k];
	}

	// --- N O R M A L   S T I F F N E S S ---
	double sum;
	for (k=0; k<ndof; ++k)
		for (l=0; l<ndof; ++l)
			{
				sum = 0;

				sum = Mi[0][0]*Nb1[k]*Nb1[l]+Mi[0][1]*(Nb1[k]*Nb2[l]+Nb2[k]*Nb1[l])+Mi[1][1]*Nb2[k]*Nb2[l];
				sum *= gap;
				sum -= D1[k]*N1[l]+D2[k]*N2[l]+N1[k]*D1[l]+N2[k]*D2[l];
				sum += K[0][1]*(D1[k]*D2[l]+D2[k]*D1[l]);
				sum *= tn*m_knmult;

				sum += eps*HEAVYSIDE(Lm+eps*gap)*N[k]*N[l];
	
				ke[k][l] = sum;
			}

	// --- T A N G E N T I A L   S T I F F N E S S ---
	// We only calculate the tangential stiffness if friction is enabled. We also
	// make sure that the gap >= 0, i.e. the node is actually in contact, otherwise
	// I've noticed that the solution can diverge quickly.
	if ((m_mu*m_epsf > 0) && (gap >=0))
	{
		// get the traction multipliers
		double Lt[2];
		Lt[0] = ss.Lt[m][0];
		Lt[1] = ss.Lt[m][1];

		// get the metric tensor and its inverse
		mat2d& Mk = ss.M[m];
		mat2d Mki = Mk.inverse();

		// get the previous isoparameteric coordinates
		double rp = ss.rsp[m][0];
		double sp = ss.rsp[m][1];

		// get the traction
		double Tt[2];
		// a. trial state
		Tt[0] = Lt[0] + Mk[0][0]*(r - rp) + Mk[0][1]*(s - sp);
		Tt[1] = Lt[1] + Mk[1][0]*(r - rp) + Mk[1][1]*(s - sp);

		double TMT = Tt[0]*(Mki[0][0]*Tt[0] + Mki[0][1]*Tt[1]) + Tt[1]*(Mki[1][0]*Tt[0] + Mki[1][1]*Tt[1]);

		// calculate the normalized traction vector
		vec3d pt = tau1*Tt[0] + tau2*Tt[1];
		pt.unit();
		
		// calculate the covariant version
		double Pt[2] = {Tt[0]/sqrt(TMT), Tt[1]/sqrt(TMT)};
		double Ptc[2];
		Ptc[0] = Mki[0][0]*Pt[0]+Mki[0][1]*Pt[1];
		Ptc[1] = Mki[1][0]*Pt[0]+Mki[1][1]*Pt[1];

		//b. return map
		bool bstick = true;
		double phi = sqrt(TMT) - m_mu*tn;
		if (phi > 0)
		{
			Tt[0] = m_mu*Tt[0]/sqrt(TMT);
			Tt[1] = m_mu*Tt[1]/sqrt(TMT);
			bstick = false;
		}

		// we need to define additional arrays for the tangential
		// contribution of the contact stiffness
		double T11[15] = {0}, T12[15] = {0}, T21[15] = {0}, T22[15] = {0};	// Tab matrix
		double N11[15] = {0}, N12[15] = {0}, N21[15] = {0}, N22[15] = {0};	// Nab matrix
		double P1[15] = {0}, P2[15] = {0};	// P arrays
		double Tb11[15], Tb21[15], Tb12[15], Tb22[15]; // Tbar matrix
		double Pb1[15], Pb2[15]; // Pbar array

		for (k=0; k<nmeln; ++k)
		{
			T11[3*(k+1)  ] = -Hr[k]*tau1.x;
			T11[3*(k+1)+1] = -Hr[k]*tau1.y;
			T11[3*(k+1)+2] = -Hr[k]*tau1.z;

			T12[3*(k+1)  ] = -Hs[k]*tau1.x;
			T12[3*(k+1)+1] = -Hs[k]*tau1.y;
			T12[3*(k+1)+2] = -Hs[k]*tau1.z;

			T21[3*(k+1)  ] = -Hr[k]*tau2.x;
			T21[3*(k+1)+1] = -Hr[k]*tau2.y;
			T21[3*(k+1)+2] = -Hr[k]*tau2.z;

			T22[3*(k+1)  ] = -Hs[k]*tau2.x;
			T22[3*(k+1)+1] = -Hs[k]*tau2.y;
			T22[3*(k+1)+2] = -Hs[k]*tau2.z;

			if (nmeln==4)
			{
				N12[3*(k+1)  ] = N21[3*(k+1)  ] = -0.25*nu.x;
				N12[3*(k+1)+1] = N21[3*(k+1)+1] = -0.25*nu.y;
				N12[3*(k+1)+2] = N21[3*(k+1)+2] = -0.25*nu.z;
			}

			P1[3*(k+1)  ] = -Hr[k]*pt.x;
			P1[3*(k+1)+1] = -Hr[k]*pt.y;
			P1[3*(k+1)+2] = -Hr[k]*pt.z;

			P2[3*(k+1)  ] = -Hs[k]*pt.x;
			P2[3*(k+1)+1] = -Hs[k]*pt.y;
			P2[3*(k+1)+2] = -Hs[k]*pt.z;
		}

		vec3d g12(0,0,0);
		if (nmeln == 4) g12 = rt[0]*Grs[0]+rt[1]*Grs[1]+rt[2]*Grs[2]+rt[3]*Grs[3];
		double gt1 = g12*tau1;
		double gt2 = g12*tau2;
		double gp = g12*pt;

		for (k=0; k<ndof; ++k)
		{
			Tb11[k] = T11[k] - gt1*D2[k];
			Tb12[k] = T12[k] - gt1*D1[k];

			Tb21[k] = T21[k] - gt2*D2[k];
			Tb22[k] = T22[k] - gt2*D1[k];

			Pb1[k] = P1[k] - gp*D2[k];
			Pb2[k] = P2[k] - gp*D1[k];
		}

		// raise the indices of A
		double Ac[2][2];
		for (k=0; k<2; ++k)
			for (l=0; l<2; ++l)
			{
				Ac[k][l] = 0;
				for (i=0; i<2; ++i)
					for (j=0; j<2; ++j) Ac[k][l] += Mki[k][i]*Mki[l][j]*A[i][j];
			}

		vec3d Hrs[2][2] = {{vec3d(0,0,0),vec3d(0,0,0)},{vec3d(0,0,0),vec3d(0,0,0)}};
		if (nmeln == 4)
		{
			Hrs[0][1] = Hrs[1][0] = rt[0]*Grs[0]+rt[1]*Grs[1]+rt[2]*Grs[2]+rt[3]*Grs[3];
		}

		// calculate stiffness matrix
		// NOTE: I think I need Mi (iso Mki) for KT1 and KT2. I only need Mk (iso M) for the direct stiffnesses
		double kij;
		for (i=0; i<ndof; ++i)
			for (j=0; j<ndof; ++j)
			{
				// KT1
				kij  = T11[i]*D1[j] + T12[i]*D2[j];
				kij += D1[i]*T11[j] + D2[i]*T12[j];
				kij -= (Hrs[0][1]*tau1)*D1[i]*D2[j] + (Hrs[1][0]*tau1)*D2[i]*D1[j];
				kij += Tb11[i]*D1[j] + Tb21[i]*D2[j];
				kij += D1[i]*Tb11[j] + D2[i]*Tb21[j];
				kij += gap*(N11[i]*D1[j] + N12[i]*D2[j] + D1[i]*N11[j] + D2[i]*N12[j]);
				kij -= N[i]*Nb1[j] - Nb1[i]*N[j];
				kij -= T1[i]*Mi[0][0]*Tb11[j] + T1[i]*Mi[0][1]*Tb21[j] + T2[i]*Mi[1][0]*Tb11[j] + T2[i]*Mi[1][1]*Tb21[j];
				kij -= Tb11[i]*Mi[0][0]*T1[j] + Tb21[i]*Mi[0][1]*T1[j] + Tb11[i]*Mi[1][0]*T2[j] + Tb21[i]*Mi[1][1]*T2[j];

				ke[i][j] += m_ktmult*(Tt[0]*Ac[0][0] + Tt[1]*Ac[1][0])*kij;

				// KT2
				kij  = T21[i]*D1[j] + T22[i]*D2[j];
				kij += D1[i]*T21[j] + D2[i]*T22[j];
				kij -= (Hrs[0][1]*tau2)*D1[i]*D2[j] + (Hrs[1][0]*tau2)*D2[i]*D1[j];
				kij += Tb12[i]*D1[j] + Tb22[i]*D2[j];
				kij += D1[i]*Tb12[j] + D2[i]*Tb22[j];
				kij += gap*(N21[i]*D1[j] + N22[i]*D2[j] + D1[i]*N21[j] + D2[i]*N22[j]);
				kij -= N[i]*Nb2[j] - Nb2[i]*N[j];
				kij -= T1[i]*Mi[0][0]*Tb12[j] + T1[i]*Mi[0][1]*Tb22[j] + T2[i]*Mi[1][0]*Tb12[j] + T2[i]*Mi[1][1]*Tb22[j];
				kij -= Tb12[i]*Mi[0][0]*T1[j] + Tb22[i]*Mi[0][1]*T1[j] + Tb12[i]*Mi[1][0]*T2[j] + Tb22[i]*Mi[1][1]*T2[j];

				ke[i][j] += m_ktmult*(Tt[0]*Ac[0][1] + Tt[1]*Ac[1][1])*kij;

				// kdirect
				if (bstick)
				{
					kij = Mk[0][0]*D1[i]*D1[j] + Mk[0][1]*D1[i]*D2[j] + Mk[1][0]*D2[i]*D1[j] + Mk[1][1]*D2[i]*D2[j];
					ke[i][j] += m_epsf*kij;
				}
				else
				{
					kij  = (1.0 - Ptc[0]*Pt[0])*(Mk[0][0]*D1[i]*D1[j]+Mk[0][1]*D1[i]*D2[j]);
					kij += (    - Ptc[0]*Pt[1])*(Mk[1][0]*D1[i]*D1[j]+Mk[1][1]*D1[i]*D2[j]);
					kij += (    - Ptc[1]*Pt[0])*(Mk[0][0]*D2[i]*D1[j]+Mk[0][1]*D2[i]*D2[j]);
					kij += (1.0 - Ptc[1]*Pt[1])*(Mk[1][0]*D2[i]*D1[j]+Mk[1][1]*D2[i]*D2[j]);
					
					ke[i][j] += m_ktmult*m_epsf*m_mu*tn/sqrt(TMT)*kij;
				}
			}
	}
}

//-----------------------------------------------------------------------------

bool FESlidingInterface::Augment(int naug)
{
	int i;
	double Ln;
	double Lt[2];
	bool bconv = true;
	mat2d Mi;

	static double normg0 = 0;

	// make sure we need to augment
	if (!m_blaugon) return true;

	// penalty factor
	double eps, scale = Penalty();


	// --- c a l c u l a t e   i n i t i a l   n o r m s ---
	// a. normal component
	double normL0 = 0;
	for (i=0; i<m_ss.Nodes(); ++i)	normL0 += m_ss.Lm[i]*m_ss.Lm[i];
	for (i=0; i<m_ms.Nodes(); ++i)	normL0 += m_ms.Lm[i]*m_ms.Lm[i];

	// b. tangential component
	if (m_mu*m_epsf > 0)
	{
		for (i=0; i<m_ss.Nodes(); ++i)
		{
			if (m_ss.pme[i])
			{
				Lt[0] = m_ss.Lt[i][0];
				Lt[1] = m_ss.Lt[i][1];
				mat2d& M = m_ss.M[i];
				Mi = M.inverse();
				normL0 += Lt[0]*(Mi[0][0]*Lt[0] + Mi[0][1]*Lt[1]) + Lt[1]*(Mi[1][0]*Lt[0] + Mi[1][1]*Lt[1]);
			}
		}

		for (i=0; i<m_ms.Nodes(); ++i)
		{
			if (m_ms.pme[i])
			{
				Lt[0] = m_ms.Lt[i][0];
				Lt[1] = m_ms.Lt[i][1];
				mat2d& M = m_ms.M[i];
				Mi = M.inverse();
				normL0 += Lt[0]*(Mi[0][0]*Lt[0] + Mi[0][1]*Lt[1]) + Lt[1]*(Mi[1][0]*Lt[0] + Mi[1][1]*Lt[1]);
			}
		}
	}
	normL0 = sqrt(normL0);

	// --- c a l c u l a t e   c u r r e n t   n o r m s ---
	// a. normal component
	double normL1 = 0;	// force norm
	double normg1 = 0;	// gap norm
	int N = 0;
	for (i=0; i<m_ss.Nodes(); ++i)
	{
		eps = m_ss.eps[i]*scale;

		// update Lagrange multipliers
		Ln = m_ss.Lm[i] + eps*m_ss.gap[i];
		Ln = MBRACKET(Ln);

		normL1 += Ln*Ln;

		if (m_ss.gap[i] > 0)
		{
			normg1 += m_ss.gap[i]*m_ss.gap[i];
			++N;
		}
	}	

	for (i=0; i<m_ms.Nodes(); ++i)
	{
		eps = m_ms.eps[i]*scale;

		// update Lagrange multipliers
		Ln = m_ms.Lm[i] + eps*m_ms.gap[i];
		Ln = MBRACKET(Ln);

		normL1 += Ln*Ln;
		if (m_ms.gap[i] > 0)
		{
			normg1 += m_ms.gap[i]*m_ms.gap[i];
			++N;
		}
	}
	if (N == 0) N=1;

	// b. tangential component
	if (m_mu*m_epsf > 0)
	{
		double r, s, rp, sp;
		for (i=0; i<m_ss.Nodes(); ++i)
		{
			if (m_ss.pme[i])
			{
				r = m_ss.rs[i][0];
				s = m_ss.rs[i][1];
				rp = m_ss.rsp[i][0];
				sp = m_ss.rsp[i][1];
				Ln = m_ss.Lm[i];
	
				mat2d& Mk = m_ss.M[i];
				Mi = Mk.inverse();

				Lt[0] = m_ss.Lt[i][0] + m_epsf*(Mk[0][0]*(r - rp) + Mk[0][1]*(s - sp));
				Lt[1] = m_ss.Lt[i][1] + m_epsf*(Mk[1][0]*(r - rp) + Mk[1][1]*(s - sp));

				double TMT = Lt[0]*(Mi[0][0]*Lt[0]+Mi[0][1]*Lt[1])+Lt[1]*(Mi[1][0]*Lt[0]+Mi[1][1]*Lt[1]);
				double phi = sqrt(TMT) - m_mu*Ln;

				// b. return map
				if (phi > 0)
				{
					Lt[0] = m_mu*Ln*Lt[0]/sqrt(TMT);
					Lt[1] = m_mu*Ln*Lt[1]/sqrt(TMT);
				}

				normL1 += Lt[0]*(Mi[0][0]*Lt[0] + Mi[0][1]*Lt[1]) + Lt[1]*(Mi[1][0]*Lt[0] + Mi[1][1]*Lt[1]);
			}
		}

		for (i=0; i<m_ms.Nodes(); ++i)
		{
			if (m_ms.pme[i])
			{
				r = m_ms.rs[i][0];
				s = m_ms.rs[i][1];
				rp = m_ms.rsp[i][0];
				sp = m_ms.rsp[i][1];
				Ln = m_ms.Lm[i];

				mat2d& Mk = m_ms.M[i];
				Mi = Mk.inverse();

				Lt[0] = m_ms.Lt[i][0] + m_epsf*(Mk[0][0]*(r - rp) + Mk[0][1]*(s - sp));
				Lt[1] = m_ms.Lt[i][1] + m_epsf*(Mk[1][0]*(r - rp) + Mk[1][1]*(s - sp));

				double TMT = Lt[0]*(Mi[0][0]*Lt[0]+Mi[0][1]*Lt[1])+Lt[1]*(Mi[1][0]*Lt[0]+Mi[1][1]*Lt[1]);
				double phi = sqrt(TMT) - m_mu*Ln;

				// b. return map
				if (phi > 0)
				{
					Lt[0] = m_mu*Ln*Lt[0]/sqrt(TMT);
					Lt[1] = m_mu*Ln*Lt[1]/sqrt(TMT);
				}

				normL1 += Lt[0]*(Mi[0][0]*Lt[0] + Mi[0][1]*Lt[1]) + Lt[1]*(Mi[1][0]*Lt[0] + Mi[1][1]*Lt[1]);
			}
		}
	}

	normL1 = sqrt(normL1);
	normg1 = sqrt(normg1 / N);

	if (naug == 0) normg0 = 0;

	// get the logfile
	Logfile& log = GetLogfile();

	// calculate and print convergence norms
	double lnorm = 0, gnorm = 0;
	if (normL1 != 0) lnorm = fabs(normL1 - normL0)/normL1; else lnorm = fabs(normL1 - normL0);
	if (normg1 != 0) gnorm = fabs(normg1 - normg0)/normg1; else gnorm = fabs(normg1 - normg0);

	log.printf(" sliding interface # %d\n", m_nID);
	log.printf("                        CURRENT        REQUIRED\n");
	log.printf("    normal force : %15le", lnorm);
	if (m_atol > 0) log.printf("%15le\n", m_atol); else log.printf("       ***\n");
	log.printf("    gap function : %15le", gnorm);
	if (m_gtol > 0) log.printf("%15le\n", m_gtol); else log.printf("       ***\n");

	// check convergence
	bconv = true;
	if ((m_atol > 0) && (lnorm > m_atol)) bconv = false;
	if ((m_gtol > 0) && (gnorm > m_gtol)) bconv = false;
	if (m_naugmin > naug) bconv = false;
	if (m_naugmax <= naug) bconv = true;
		
	if (bconv == false)
	{
		// we did not converge so update multipliers
		for (i=0; i<m_ss.Nodes(); ++i)
		{
			eps = m_ss.eps[i]*scale;

			// update Lagrange multipliers
			Ln = m_ss.Lm[i] + eps*m_ss.gap[i];
			m_ss.Lm[i] = MBRACKET(Ln);

			if ((m_mu*m_epsf > 0) && (m_ss.pme[i]))
			{
				// update the metrics
				FESurfaceElement& mel = dynamic_cast<FESurfaceElement&>(*m_ss.pme[i]);
				m_pfem->m_mesh.UnpackElement(mel);

				double r = m_ss.rs[i][0], s = m_ss.rs[i][1];
				double rp = m_ss.rsp[i][0], sp = m_ss.rsp[i][1];

				double Ln = m_ss.Lm[i];

				mat2d Mk = m_ss.M[i];
				mat2d Mki = Mk.inverse();
				
			
				// update traction multipliers
				// a. trial state
				Lt[0] = m_ss.Lt[i][0] + m_epsf*(Mk[0][0]*(r - rp) + Mk[0][1]*(s - sp));
				Lt[1] = m_ss.Lt[i][1] + m_epsf*(Mk[1][0]*(r - rp) + Mk[1][1]*(s - sp));

				double TMT = Lt[0]*(Mki[0][0]*Lt[0]+Mki[0][1]*Lt[1])+Lt[1]*(Mki[1][0]*Lt[0]+Mki[1][1]*Lt[1]);
				assert(TMT >= 0);

				double phi = sqrt(TMT) - m_mu*Ln;

				// b. return map
				if (phi > 0)
				{
					Lt[0] = m_mu*Ln*Lt[0]/sqrt(TMT);
					Lt[1] = m_mu*Ln*Lt[1]/sqrt(TMT);
				}

				m_ss.M[i] = m_ss.Metric0(mel, r, s);

				m_ss.Lt[i][0] = Lt[0];
				m_ss.Lt[i][1] = Lt[1];
			
			}
		}	

		for (i=0; i<m_ms.Nodes(); ++i)
		{
			eps = m_ms.eps[i]*scale;

			// update Lagrange multipliers
			Ln = m_ms.Lm[i] + eps*m_ms.gap[i];
			m_ms.Lm[i] = MBRACKET(Ln);

			if ((m_mu*m_epsf > 0) && (m_ms.pme[i]))
			{
				// update the metrics
				FESurfaceElement& mel = dynamic_cast<FESurfaceElement&>(*m_ms.pme[i]);
				m_pfem->m_mesh.UnpackElement(mel);

				double r = m_ms.rs[i][0], s = m_ms.rs[i][1];
				double rp = m_ms.rsp[i][0], sp = m_ms.rsp[i][1];

				double Ln = m_ms.Lm[i];

				mat2d Mk = m_ms.M[i];
				mat2d Mki = Mk.inverse();
				
				// update traction multipliers
				// a. trial state
				double Lt[2];
				Lt[0] = m_ms.Lt[i][0] + m_epsf*(Mk[0][0]*(r - rp) + Mk[0][1]*(s - sp));
				Lt[1] = m_ms.Lt[i][1] + m_epsf*(Mk[1][0]*(r - rp) + Mk[1][1]*(s - sp));

				double TMT = Lt[0]*(Mki[0][0]*Lt[0]+Mki[0][1]*Lt[1])+Lt[1]*(Mki[1][0]*Lt[0]+Mki[1][1]*Lt[1]);
				assert(TMT >= 0);

				double phi = sqrt(TMT) - m_mu*Ln;

				// b. return map
				if (phi > 0)
				{
					Lt[0] = m_mu*Ln*Lt[0]/sqrt(TMT);
					Lt[1] = m_mu*Ln*Lt[1]/sqrt(TMT);
				}

				m_ms.M[i] = m_ms.Metric0(mel, r, s);

				m_ms.Lt[i][0] = Lt[0];
				m_ms.Lt[i][1] = Lt[1];
			}
		}
	}

	if (bconv)
	{
		m_ss.rsp = m_ss.rs;
		m_ms.rsp = m_ms.rs;
	}

	// store the last gap norm
	normg0 = normg1;

	return bconv;
}

//-----------------------------------------------------------------------------
//! This function transforms friction data between two master segments

void FESlidingInterface::MapFrictionData(int inode, FEContactSurface& ss, FEContactSurface& ms, FESurfaceElement &en, FESurfaceElement &eo, vec3d &q)
{
	// first we find the projection of the old point on the new segment
	double r = ss.rs[inode][0];
	double s = ss.rs[inode][1];
	double rp = ss.rsp[inode][0], ro = rp;
	double sp = ss.rsp[inode][1], so = sp;
	vec3d xn = ms.PointOnSurface(eo, rp, sp);
	vec3d qn;
	qn = ms.ProjectToSurface(en, xn, rp, sp);
	ss.rsp[inode][0] = rp;
	ss.rsp[inode][1] = sp;

	// next, we transform the frictional traction
	// since these tractions are decomposed in the local 
	// element coordinate system, we have to do a coordinate transformation
	// note that this transformation needs to be done in curvilinear
	// coordinates since the base vectors may not be orthonormal. Also
	// note that we are doing this in the reference configuration
	vec3d to[2], tn[2];
	ms.ContraBaseVectors0(eo, ro, so, to);
	ms.CoBaseVectors0(en, r, s, tn);

	double Lt[2];
	Lt[0] = ss.Lt[inode][0];
	Lt[1] = ss.Lt[inode][1];
	
	vec3d t;
	t = to[0]*Lt[0] + to[1]*Lt[1];

	Lt[0] = t*tn[0];
	Lt[1] = t*tn[1];

	ss.Lt[inode][0] = Lt[0];
	ss.Lt[inode][1] = Lt[1];
}

//-----------------------------------------------------------------------------

void FESlidingInterface::Serialize(Archive& ar)
{
	int j, k, n, mat;

	if (ar.IsSaving())
	{
		ar << m_npass;
		ar << m_naugmax;
		ar << m_naugmin;
		ar << m_gtol;
		ar << m_atol;
		ar << m_blaugon;
		ar << m_ktmult;
		ar << m_knmult;
		ar << m_stol;
		ar << m_bautopen;
		ar << m_eps;
		ar << m_mu;
		ar << m_epsf;
		ar << m_nsegup;
		ar << m_nplc;

		for (j=0; j<2; ++j)
		{
			FEContactSurface& s = (j==0? m_ss : m_ms);

			int ne = s.Elements();
			ar << ne;

			for (k=0; k<ne; ++k)
			{
				FESurfaceElement& el = s.Element(k);
				ar << el.Type();
				ar << el.GetMatID() << el.m_nID << el.m_nrigid;
				ar << el.m_node;
				ar << el.m_lnode;
			}

			ar << s.gap;
			ar << s.nu;
			ar << s.rs;
			ar << s.rsp;
			ar << s.Lm;
			ar << s.M;
			ar << s.Lt;
			ar << s.off;
			ar << s.eps;
		}
	}
	else
	{
		ar >> m_npass;
		ar >> m_naugmax;
		ar >> m_naugmin;
		ar >> m_gtol;
		ar >> m_atol;
		ar >> m_blaugon;
		ar >> m_ktmult;
		ar >> m_knmult;
		ar >> m_stol;
		ar >> m_bautopen;
		ar >> m_eps;
		ar >> m_mu;
		ar >> m_epsf;
		ar >> m_nsegup;
		ar >> m_nplc;

		if (m_nplc >= 0) m_pplc = m_pfem->GetLoadCurve(m_nplc);

		for (j=0; j<2; ++j)
		{
			FEContactSurface& s = (j==0? m_ss : m_ms);

			int ne=0;
			ar >> ne;
			s.Create(ne);

			for (k=0; k<ne; ++k)
			{
				FESurfaceElement& el = s.Element(k);

				ar >> n;
				el.SetType(n);

				ar >> mat >> el.m_nID >> el.m_nrigid;
				ar >> el.m_node;
				ar >> el.m_lnode;

				el.SetMatID(mat);
			}

			// initialize surface
			s.Init();

			// read the contact data
			// Note that we do this after Init() since this data gets 
			// initialized to zero there
			ar >> s.gap;
			ar >> s.nu;
			ar >> s.rs;
			ar >> s.rsp;
			ar >> s.Lm;
			ar >> s.M;
			ar >> s.Lt;
			ar >> s.off;
			ar >> s.eps;
		}
	}
}
