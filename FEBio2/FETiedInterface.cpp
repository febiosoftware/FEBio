// FETiedInterface.cpp: implementation of the FETiedInterface class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FETiedInterface.h"
#include "fem.h"
#include "FESolidSolver.h"
#include "FECore/FENNQuery.h"
#include "log.h"

///////////////////////////////////////////////////////////////////////////////
// FETiedInterface
///////////////////////////////////////////////////////////////////////////////


//-----------------------------------------------------------------------------
//! constructor

FETiedInterface::FETiedInterface(FEM* pfem) : FEContactInterface(pfem), ss(&pfem->m_mesh), ms(&pfem->m_mesh)
{
	static int count = 1;
	m_ntype = FE_CONTACT_TIED;

	m_nplc = -1;
	m_pplc = 0;

	ss.SetSibling(&ms);
	ms.SetSibling(&ss);

	m_nID = count++;
}

//-----------------------------------------------------------------------------
//! initialization

void FETiedInterface::Init()
{
	// create the surfaces
	ss.Init();
	ms.Init();

	// project slave surface onto master surface
	ProjectSurface(ss, ms, false);

	// set penalty load curve
	if (m_nplc >= 0) m_pplc = m_pfem->GetLoadCurve(m_nplc);
}

//-----------------------------------------------------------------------------
//! update

void FETiedInterface::Update()
{
	int i, l, n;
	double r, s, H[4];
	vec3d T;
	vec3d y, q, rt;

	// get the mesh
	FEMesh& mesh = *ss.GetMesh();

	// update the gap functions
	FEElement* pme;
	for (i=0; i<ss.Nodes(); ++i)
	{
		pme = ss.m_pme[i];
		if (pme)
		{
			// get the nodal position
			vec3d rt = ss.Node(i).m_rt;

			// get the natural coordinates of the slave projection
			// onto the master element
			r = ss.rs[i][0];
			s = ss.rs[i][1];

			// calculate the shape function values
			n = pme->Nodes();
			if (n == 4)
			{

				H[0] = 0.25*(1 - r)*(1 - s);
				H[1] = 0.25*(1 + r)*(1 - s);
				H[2] = 0.25*(1 + r)*(1 + s);
				H[3] = 0.25*(1 - r)*(1 + s);
			}
			else if (n == 3)
			{
				H[0] = 1.0 - r - s;
				H[1] = r;
				H[2] = s;
			}
			else assert(false);

			// calculate the slave node projection
			q = vec3d(0,0,0);
			for (l=0; l<n; ++l)
			{
				y = mesh.Node( pme->m_node[l] ).m_rt;
				q += y*H[l];
			}

			// calculate the gap function
			ss.gap[i] = rt - q;
		}
	}
}

//-----------------------------------------------------------------------------
//! project surface

void FETiedInterface::ProjectSurface(FETiedContactSurface& ss, FETiedContactSurface& ms, bool bmove)
{
	int i, j, l, m, n, np1, nm1, a;

	int mn, ne;
	double r, s;

	int nval;
	FEElement** pe;
	FEElement* pme;

	double H[4];
	double l1, l2, N;

	vec3d q, r0;

	vec3d y[4];

	// plane normal
	vec3d np;
	// projection onto plane
	vec3d rp;

	// neirest neigbour query
	FENNQuery query(&ms);

	query.Init();

	// get the mesh
	FEMesh& mesh = *ss.GetMesh();

	const double eps = 0.01;

	// loop over all slave nodes
	for (i=0; i<ss.Nodes(); ++i)
	{
		FENode& node = ss.Node(i);

		a = ss.node[i];

		// get the nodal position
		vec3d rs = node.m_rt;

		// let's find the closest master node
		mn = query.Find(rs);

		// now that we found the closest master node, lets see if we can find 
		// the best master element
		pme = 0;

		// mn is a local index, so get the global node number too
		m = ms.node[mn];
		r0 = mesh.Node(m).m_rt;

		nval = ms.m_NEL.Valence(mn);
		pe = ms.m_NEL.ElementList(mn);
	
		for (j=0; j<nval; ++j)
		{
			// get the master element
			FESurfaceElement& el = dynamic_cast<FESurfaceElement&> (*pe[j]);

			N = el.Nodes();
			if (N == 4)
			{
				// find which node this is
				n = -1;
				if (el.m_node[0] == m) n = 0;
				else if (el.m_node[1] == m) n = 1;
				else if (el.m_node[2] == m) n = 2;
				else if (el.m_node[3] == m) n = 3;
				assert(n >= 0);

				np1 = (n+1)%4;
				nm1 = (n==0? 3 : n-1);

				// approximate element surface by a plane
				vec3d r1 = mesh.Node(el.m_node[np1]).m_rt - r0;
				vec3d r2 = mesh.Node(el.m_node[nm1]).m_rt - r0;

				l1 = r1.norm();
				l2 = r2.norm();

				np = r1^r2;
				np.unit();

				// project rs to the plane
				rp = rs - np*(np*(rs - r0));
	
				// get the plane coordinates
				r = ((rp - r0)*r1)/(l1*l1);
				s = ((rp - r0)*r2)/(l2*l2);
		
				if ((r >= -eps) && (s >= -eps))
				{
					// we have a winner !
					pme = pe[j];
					break;
				}
			}
			else if (N==3)
			{
				// approximate element surface by a plane
				vec3d rc = mesh.Node(el.m_node[0]).m_rt;
				vec3d e1 = mesh.Node(el.m_node[1]).m_rt - rc;
				vec3d e2 = mesh.Node(el.m_node[2]).m_rt - rc;
	
				l1 = e1.norm();
				l2 = e2.norm();
	
				np = e1^e2;
				np.unit();

				// project rs to the plane
				rp = rs - np*(np*(rs - rc));

				// get the plane coordinates
				r = ((rp - r0)*e1)/(l1*l1);
				s = ((rp - r0)*e2)/(l2*l2);
	
				if ((r >= -eps) && (s >= -eps) && (r+s<=1+eps))
				{
					// we have a winner !
					pme = pe[j];
					break;
				}
			}
		}

		ss.m_pme[i] = dynamic_cast<FESurfaceElement*>(pme);
		if (pme != 0)
		{
			// If we found a master element, we calculate the intersection
			// more accurately.
			FESurfaceElement& mel = dynamic_cast<FESurfaceElement&> (*pme);

			ne = mel.Nodes();

			if (ne == 4)
			{
				double fr, fs;

				// let's find the natural coordinates in the master element
				// for now we simply approximate the natural coordinate using the plane coordinates
				switch(n)
				{
				case 0:
					fr = -1 + 2*r;
					fs = -1 + 2*s;
					break;
				case 1:
					fr =  1 - 2*s;
					fs = -1 + 2*r;
					break;
				case 2:
					fr =  1 - 2*r;
					fs =  1 - 2*s;
					break;
				case 3:
					fr = -1 + 2*s;
					fs =  1 - 2*r;
					break;
				}

				// improve estimate of (r,s)
				r = fr;
				s = fs;
				q = ms.ProjectToSurface(mel, rs, fr, fs);

				if ((fr < -1-eps) || (fr > 1+eps) || (fs < -1-eps) || (fs > 1+eps))
				{
					// the projection has failed
					fr = r;
					fs = s;

					for (l=0; l<4; ++l)
						y[l] = mesh.Node( mel.m_node[l] ).m_rt;

					H[0] = 0.25*(1 - fr)*(1 - fs);
					H[1] = 0.25*(1 + fr)*(1 - fs);
					H[2] = 0.25*(1 + fr)*(1 + fs);
					H[3] = 0.25*(1 - fr)*(1 + fs);

					q = vec3d(0,0,0);
					for (l=0; l<4; ++l)	q += y[l]*H[l];
				}
				else
				{
					H[0] = 0.25*(1 - fr)*(1 - fs);
					H[1] = 0.25*(1 + fr)*(1 - fs);
					H[2] = 0.25*(1 + fr)*(1 + fs);
					H[3] = 0.25*(1 - fr)*(1 + fs);
				}

				ss.rs[i][0] = fr;
				ss.rs[i][1] = fs;
			}
			else if (ne == 3)
			{
				for (l=0; l<3; ++l)
					y[l] = mesh.Node( mel.m_node[l] ).m_rt;

				H[0] = 1.0 - r - s;
				H[1] = r;
				H[2] = s;

				q = vec3d(0,0,0);
				for (l=0; l<3; ++l)	q += y[l]*H[l];

				ss.rs[i][0] = r;
				ss.rs[i][1] = s;
			}

			// calculate gap
			ss.gap[i] = rs - q;
			if (bmove && (ss.gap[i].norm()>0))
			{
				node.m_r0 = node.m_rt = q;
				ss.gap[i] = vec3d(0,0,0);
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! This function calculates the contact forces for a tied interface.

void FETiedInterface::ContactForces(vector<double>& F)
{
	int j, k, l, m, n;
	int nseln, nmeln;

	double *Gr, *Gs;

	// jacobian
	double detJ;

	vec3d dxr, dxs;
	vec3d* rt, *r0;
	double* w;

	// natural coordinates of slave node in master element
	double r, s;

	// contact force
	vec3d tc;

	// shape function values
	double N[4];

	// element contact force vector
	vector<double> fe;

	// the lm array for this force vector
	vector<int> lm;

	// the en array
	vector<int> en;

	vector<int> sLM;
	vector<int> mLM;

	FEM& fem = dynamic_cast<FEM&>(*m_pfem);
	FEMesh& mesh = m_pfem->m_mesh;

	FESolidSolver* psolver = dynamic_cast<FESolidSolver*>(fem.m_pStep->m_psolver);
	
	// loop over all slave facets
	int ne = ss.Elements();
	for (j=0; j<ne; ++j)
	{
		// get the slave element
		FESurfaceElement& sel = ss.Element(j);
		ss.UnpackElement(sel);

		sLM = sel.LM();

		nseln = sel.Nodes();

		rt = sel.rt();
		r0 = sel.r0();
		w = sel.GaussWeights();

		// loop over slave element nodes (which are the integration points as well)
		for (n=0; n<nseln; ++n)
		{
			Gr = sel.Gr(n);
			Gs = sel.Gs(n);

			m = sel.m_lnode[n];

			// see if this node's constraint is active
			// that is, if it has a master element associated with it
			// TODO: is this a good way to test for an active constraint
			// The rigid wall criteria seems to work much better.
			if (ss.m_pme[m] != 0)
			{
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

				detJ = (dxr ^ dxs).norm();

				// get slave node contact force
				tc = ss.Lm[m] + ss.gap[m]*m_eps;

				// get the master element
				FESurfaceElement& mel = dynamic_cast<FESurfaceElement&> (*ss.m_pme[m]);
				ms.UnpackElement(mel, FE_UNPACK_LM);

				mLM = mel.LM();

				// we must unpack the slave element again
				ss.UnpackElement(sel);

				nmeln = mel.Nodes();

				// isoparametric coordinates of the projected slave node
				// onto the master element
				r = ss.rs[m][0];
				s = ss.rs[m][1];

				// get the master shape function values at this slave node
				if (nmeln == 4)
				{
					// quadrilateral
					N[0] = 0.25*(1-r)*(1-s);
					N[1] = 0.25*(1+r)*(1-s);
					N[2] = 0.25*(1+r)*(1+s);
					N[3] = 0.25*(1-r)*(1+s);
				}
				else if (nmeln == 3)
				{
					// triangle
					N[0] = 1 - r - s;
					N[1] = r;
					N[2] = s;
				}

				// calculate force vector
				fe.resize(3*(nmeln+1));
				fe[0] = -detJ*w[n]*tc.x;
				fe[1] = -detJ*w[n]*tc.y;
				fe[2] = -detJ*w[n]*tc.z;
				for (l=0; l<nmeln; ++l)
				{	
					fe[3*(l+1)  ] = detJ*w[n]*tc.x*N[l];
					fe[3*(l+1)+1] = detJ*w[n]*tc.y*N[l];
					fe[3*(l+1)+2] = detJ*w[n]*tc.z*N[l];
				}
	
				// fill the lm array
				lm.resize(3*(nmeln+1));
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
				en.resize(nmeln+1);
				en[0] = sel.m_node[n];
				for (l=0; l<nmeln; ++l) en[l+1] = mel.m_node[l];

				// assemble into global force vector
				psolver->AssembleResidual(en, lm, fe, F);
			}
		}
	}
}

//-----------------------------------------------------------------------------

void FETiedInterface::ContactStiffness()
{
	int j, k, l, n, m;
	int nseln, nmeln, ndof;

	matrix ke;

	vector<int> lm(15);
	vector<int> en(5);

	double *Gr, *Gs, *w;
	vec3d *rt, *r0;

	vec3d rtm[4];

	double detJ, r, s;
	vec3d dxr, dxs;
	double H[4];

	vec3d gap, Lm, tc;

	// penalty factor
	double eps = Penalty();

	// curvature tensor K
	double K[2][2] = {0};

//	double scale = -0.0035*m_fem.m_mesh.GetBoundingBox().radius();

	vector<int> sLM;
	vector<int> mLM;

	FEM& fem = dynamic_cast<FEM&>(*m_pfem);
	FEMesh& mesh = m_pfem->m_mesh;

	FESolidSolver* psolver = dynamic_cast<FESolidSolver*>(fem.m_pStep->m_psolver);

	// loop over all slave elements
	int ne = ss.Elements();
	for (j=0; j<ne; ++j)
	{
		FESurfaceElement& se = ss.Element(j);
		ss.UnpackElement(se);

		sLM = se.LM();

		nseln = se.Nodes();

		r0 = se.r0();
		rt = se.rt();
		w = se.GaussWeights();

		// loop over all integration points (that is nodes)
		for (n=0; n<nseln; ++n)
		{
			Gr = se.Gr(n);
			Gs = se.Gs(n);

			m = se.m_lnode[n];

			// see if this node's constraint is active
			// that is, if it has a master element associated with it
			if (ss.m_pme[m] != 0)
			{
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

				detJ = (dxr ^ dxs).norm();

				// get the master element
				FESurfaceElement& me = dynamic_cast<FESurfaceElement&> (*ss.m_pme[m]);
				ms.UnpackElement(me, FE_UNPACK_LM);

				mLM = me.LM();
				nmeln = me.Nodes();

				// get the master element node positions
				for (k=0; k<nmeln; ++k) rtm[k] = me.rt()[k];

				// we must unpack the slave element again
				ss.UnpackElement(se);

				// slave node natural coordinates in master element
				r = ss.rs[m][0];
				s = ss.rs[m][1];

				// slave gap
				gap = ss.gap[m];

				// lagrange multiplier
				Lm = ss.Lm[m];

				// get slave node normal force
				tc = ss.Lm[m] + ss.gap[m]*eps; //ss.T[m];

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

				// number of degrees of freedom
				ndof = 3*(1 + nmeln);

				// fill stiffness matrix
				ke.Create(ndof, ndof); ke.zero();
				ke[0][0] = w[n]*detJ*eps;
				ke[1][1] = w[n]*detJ*eps;
				ke[2][2] = w[n]*detJ*eps;
				for (k=0; k<nmeln; ++k)
				{
					ke[0][3+3*k  ] = -w[n]*detJ*eps*H[k];
					ke[1][3+3*k+1] = -w[n]*detJ*eps*H[k];
					ke[2][3+3*k+2] = -w[n]*detJ*eps*H[k];

					ke[3+3*k  ][0] = -w[n]*detJ*eps*H[k];
					ke[3+3*k+1][1] = -w[n]*detJ*eps*H[k];
					ke[3+3*k+2][2] = -w[n]*detJ*eps*H[k];
				}
				for (k=0; k<nmeln; ++k)
					for (l=0; l<nmeln; ++l)
					{
						ke[3+3*k  ][3+3*l  ] = w[n]*detJ*eps*H[k]*H[l];
						ke[3+3*k+1][3+3*l+1] = w[n]*detJ*eps*H[k]*H[l];
						ke[3+3*k+2][3+3*l+2] = w[n]*detJ*eps*H[k]*H[l];
					}

				// create lm array
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
				en.resize(nmeln+1);
				en[0] = se.m_node[n];
				for (k=0; k<nmeln; ++k) en[k+1] = me.m_node[k];
						
				// assemble stiffness matrix
				psolver->AssembleStiffness(en, lm, ke);
			}
		}
	}
}

//-----------------------------------------------------------------------------

bool FETiedInterface::Augment(int naug)
{
	int i;
	bool bconv = true;

	// penalty factor
	double eps = Penalty();

	double g;
	vec3d lm;

	// calculate initial norms
	double normL0 = 0;
	for (i=0; i<ss.Nodes(); ++i)
	{
		lm = ss.Lm[i];
		normL0 += lm*lm;
	}
	normL0 = sqrt(normL0);

	// update Lagrange multipliers and calculate current norms
	double normL1 = 0;
	double normgc = 0;
	int N = 0;
	for (i=0; i<ss.Nodes(); ++i)
	{
		lm = ss.Lm[i] + ss.gap[i]*eps;

		normL1 += lm*lm;
		if (ss.m_pme[i] != 0)
		{
			g = ss.gap[i].norm();
			normgc += g*g;
			++N;
		}
	}	
	if (N == 0) N=1;

	normL1 = sqrt(normL1);
	normgc = sqrt(normgc / N);

	// check convergence of constraints
	clog.printf(" tied interface # %d\n", m_nID);
	clog.printf("                        CURRENT        REQUIRED\n");
	double pctn = 0;
	if (fabs(normL1) > 1e-10) pctn = fabs((normL1 - normL0)/normL1);
	clog.printf("    normal force : %15le %15le\n", pctn, m_atol);
	clog.printf("    gap function : %15le       ***\n", normgc);
		
	if (pctn >= m_atol) 
	{
		bconv = false;
		for (i=0; i<ss.Nodes(); ++i)
		{
			// update Lagrange multipliers
			ss.Lm[i] = ss.Lm[i] + ss.gap[i]*eps;
		}	
	}

	return bconv;
}

//-----------------------------------------------------------------------------

void FETiedInterface::Serialize(DumpFile &ar)
{
	FEContactInterface::Serialize(ar);
	if (ar.IsSaving())
	{
		ar << m_eps;
		ar << m_atol;
		ar << m_nplc;
		ar << nse;
		ar << nme;

		ms.Serialize(ar);
		ss.Serialize(ar);
	}
	else
	{
		ar >> m_eps;
		ar >> m_atol;
		ar >> m_nplc;
		ar >> nse;
		ar >> nme;

		if (m_nplc >= 0) m_pplc = m_pfem->GetLoadCurve(m_nplc);

		ms.Serialize(ar);
		ss.Serialize(ar);
	}
}
