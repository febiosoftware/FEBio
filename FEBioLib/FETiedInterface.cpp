// FETiedInterface.cpp: implementation of the FETiedInterface class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FETiedInterface.h"
#include "FECore/FEModel.h"
#include "FECore/FENNQuery.h"
#include "log.h"

//-----------------------------------------------------------------------------
// Define sliding interface parameters
BEGIN_PARAMETER_LIST(FETiedInterface, FEContactInterface)
	ADD_PARAMETER(m_blaugon, FE_PARAM_BOOL  , "laugon"          ); 
	ADD_PARAMETER(m_atol   , FE_PARAM_DOUBLE, "tolerance"       );
	ADD_PARAMETER(m_eps    , FE_PARAM_DOUBLE, "penalty"         );
	ADD_PARAMETER(m_naugmin, FE_PARAM_INT   , "minaug"          );
	ADD_PARAMETER(m_naugmax, FE_PARAM_INT   , "maxaug"          );
	ADD_PARAMETER(m_stol   , FE_PARAM_DOUBLE, "search_tolerance");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! Constructor. Initialize default values.
FETiedInterface::FETiedInterface(FEModel* pfem) : FEContactInterface(pfem), ss(&pfem->GetMesh()), ms(&pfem->GetMesh())
{
	static int count = 1;
	m_ntype = FE_CONTACT_TIED;

	// define sibling relationships
	ss.SetSibling(&ms);
	ms.SetSibling(&ss);

	// initial parameter values
	m_blaugon = false;
	m_atol = 0.01;
	m_eps = 1.0;
	m_stol = 0.0001;
	m_naugmin = 0;
	m_naugmax = 10;

	// give this interface an ID (TODO: where is this actually used?)
	m_nID = count++;
}

//-----------------------------------------------------------------------------
//! Initialization. This function intializes the surfaces data and projects the
//! slave surface onto the master surface.
//! 
bool FETiedInterface::Init()
{
	// create the surfaces
	if (ss.Init() == false) return false;
	if (ms.Init() == false) return false;

	return true;
}

//-----------------------------------------------------------------------------
//! Interface activation
void FETiedInterface::Activate()
{
	// Don't forget to call base member!
	FEContactInterface::Activate();

	// project slave surface onto master surface
	ProjectSurface(ss, ms, false);
}

//-----------------------------------------------------------------------------
void FETiedInterface::ShallowCopy(FEContactInterface& ci)
{
	FETiedInterface& si = dynamic_cast<FETiedInterface&>(ci);
	ss.ShallowCopy(si.ss);
	ms.ShallowCopy(si.ms);
}

//-----------------------------------------------------------------------------
//! Update tied interface data. This function re-evaluates the gaps between
//! the slave node and their projections onto the master surface.
//!
void FETiedInterface::Update(int niter)
{
	// get the mesh
	FEMesh& mesh = *ss.GetMesh();

	// loop over all slave nodes
	for (int i=0; i<ss.Nodes(); ++i)
	{
		FEElement* pme = ss.m_pme[i];
		if (pme)
		{
			// get the current slave nodal position
			vec3d rt = ss.Node(i).m_rt;

			// get the natural coordinates of the slave projection
			// onto the master element
			double r = ss.m_rs[i][0];
			double s = ss.m_rs[i][1];

			// calculate the shape function values
			int ne = pme->Nodes();
			double H[4];
			if (ne == 4)
			{

				H[0] = 0.25*(1 - r)*(1 - s);
				H[1] = 0.25*(1 + r)*(1 - s);
				H[2] = 0.25*(1 + r)*(1 + s);
				H[3] = 0.25*(1 - r)*(1 + s);
			}
			else if (ne == 3)
			{
				H[0] = 1.0 - r - s;
				H[1] = r;
				H[2] = s;
			}
			else assert(false);

			// calculate the slave node projection
			vec3d q(0,0,0);
			for (int l=0; l<ne; ++l)
			{
				vec3d y = mesh.Node( pme->m_node[l] ).m_rt;
				q += y*H[l];
			}

			// calculate the gap function
			ss.m_gap[i] = rt - q;
		}
	}
}

//-----------------------------------------------------------------------------
//! project surface

void FETiedInterface::ProjectSurface(FETiedContactSurface& ss, FETiedContactSurface& ms, bool bmove)
{
	// loop over all slave nodes
	for (int i=0; i<ss.Nodes(); ++i)
	{
		// get the next node
		FENode& node = ss.Node(i);
		ss.m_pme[i] = 0;

		// get the nodal position of this slave node
		vec3d x = node.m_rt;

		// find the master element
		vec3d q; vec2d rs;
		FESurfaceElement* pme = ms.ClosestPointProjection(x, q, rs, (i==0), m_stol);
		if (pme)
		{
			// store the master element
			ss.m_pme[i] = pme;
			ss.m_rs[i][0] = rs[0];
			ss.m_rs[i][1] = rs[1];

			// calculate gap
			ss.m_gap[i] = x - q;

			// move the node if necessary
			if (bmove && (ss.m_gap[i].norm()>0))
			{
				node.m_r0 = node.m_rt = q;
				ss.m_gap[i] = vec3d(0,0,0);
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! This function calculates the contact forces for a tied interface.

void FETiedInterface::ContactForces(FEGlobalVector& R)
{
	int j, k, l, m, n;
	int nseln, nmeln;

	double *Gr, *Gs;

	// jacobian
	double detJ;

	vec3d dxr, dxs;
	vec3d rt[4], r0[4];
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

	// loop over all slave facets
	int ne = ss.Elements();
	for (j=0; j<ne; ++j)
	{
		// get the slave element
		FESurfaceElement& sel = ss.Element(j);

		// get the element's LM vector
		ss.UnpackLM(sel, sLM);

		nseln = sel.Nodes();

		for (int i=0; i<nseln; ++i)
		{
			r0[i] = ss.GetMesh()->Node(sel.m_node[i]).m_r0;
			rt[i] = ss.GetMesh()->Node(sel.m_node[i]).m_rt;
		}

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
				tc = ss.m_Lm[m] + ss.m_gap[m]*m_eps;

				// get the master element
				FESurfaceElement& mel = dynamic_cast<FESurfaceElement&> (*ss.m_pme[m]);
				ms.UnpackLM(mel, mLM);

				nmeln = mel.Nodes();

				// isoparametric coordinates of the projected slave node
				// onto the master element
				r = ss.m_rs[m][0];
				s = ss.m_rs[m][1];

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
				R.Assemble(en, lm, fe);
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! Calculate the stiffness matrix contribution.
void FETiedInterface::ContactStiffness(FENLSolver* psolver)
{
	int j, k, l, n, m;
	int nseln, nmeln, ndof;

	matrix ke;

	vector<int> lm(15);
	vector<int> en(5);

	double *Gr, *Gs, *w;
	vec3d rt[4], r0[4];

	vec3d rtm[4];

	double detJ, r, s;
	vec3d dxr, dxs;
	double H[4];

	vec3d gap, Lm, tc;

	// curvature tensor K
	double K[2][2] = {0};

//	double scale = -0.0035*m_fem.GetMesh().GetBoundingBox().radius();

	vector<int> sLM;
	vector<int> mLM;

	// loop over all slave elements
	int ne = ss.Elements();
	for (j=0; j<ne; ++j)
	{
		FESurfaceElement& se = ss.Element(j);

		// get the element's LM vector
		ss.UnpackLM(se, sLM);

		nseln = se.Nodes();

		for (int i=0; i<nseln; ++i)
		{
			r0[i] = ss.GetMesh()->Node(se.m_node[i]).m_r0;
			rt[i] = ss.GetMesh()->Node(se.m_node[i]).m_rt;
		}

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
				ms.UnpackLM(me, mLM);

				nmeln = me.Nodes();

				// get the master element node positions
				for (k=0; k<nmeln; ++k) rtm[k] = ms.GetMesh()->Node(me.m_node[k]).m_rt;

				// slave node natural coordinates in master element
				r = ss.m_rs[m][0];
				s = ss.m_rs[m][1];

				// slave gap
				gap = ss.m_gap[m];

				// lagrange multiplier
				Lm = ss.m_Lm[m];

				// get slave node normal force
				tc = ss.m_Lm[m] + ss.m_gap[m]*m_eps; //ss.T[m];

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
				ke.resize(ndof, ndof); ke.zero();
				ke[0][0] = w[n]*detJ*m_eps;
				ke[1][1] = w[n]*detJ*m_eps;
				ke[2][2] = w[n]*detJ*m_eps;
				for (k=0; k<nmeln; ++k)
				{
					ke[0][3+3*k  ] = -w[n]*detJ*m_eps*H[k];
					ke[1][3+3*k+1] = -w[n]*detJ*m_eps*H[k];
					ke[2][3+3*k+2] = -w[n]*detJ*m_eps*H[k];

					ke[3+3*k  ][0] = -w[n]*detJ*m_eps*H[k];
					ke[3+3*k+1][1] = -w[n]*detJ*m_eps*H[k];
					ke[3+3*k+2][2] = -w[n]*detJ*m_eps*H[k];
				}
				for (k=0; k<nmeln; ++k)
					for (l=0; l<nmeln; ++l)
					{
						ke[3+3*k  ][3+3*l  ] = w[n]*detJ*m_eps*H[k]*H[l];
						ke[3+3*k+1][3+3*l+1] = w[n]*detJ*m_eps*H[k]*H[l];
						ke[3+3*k+2][3+3*l+2] = w[n]*detJ*m_eps*H[k]*H[l];
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
//! Do an augmentation.
bool FETiedInterface::Augment(int naug)
{
	// make sure we need to augment
	if (!m_blaugon) return true;

	int i;

	// calculate initial norms
	double normL0 = 0;
	for (i=0; i<ss.Nodes(); ++i)
	{
		vec3d lm = ss.m_Lm[i];
		normL0 += lm*lm;
	}
	normL0 = sqrt(normL0);

	// update Lagrange multipliers and calculate current norms
	double normL1 = 0;
	double normgc = 0;
	int N = 0;
	for (i=0; i<ss.Nodes(); ++i)
	{
		vec3d lm = ss.m_Lm[i] + ss.m_gap[i]*m_eps;

		normL1 += lm*lm;
		if (ss.m_pme[i] != 0)
		{
			double g = ss.m_gap[i].norm();
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
		
	// check convergence
	bool bconv = true;
	if (pctn >= m_atol) bconv = false;
	if (naug < m_naugmin ) bconv = false;
	if (naug >= m_naugmax) bconv = true;

	if (bconv == false) 
	{
		for (i=0; i<ss.Nodes(); ++i)
		{
			// update Lagrange multipliers
			ss.m_Lm[i] = ss.m_Lm[i] + ss.m_gap[i]*m_eps;
		}	
	}

	return bconv;
}

//-----------------------------------------------------------------------------
//! Serialize the data to the archive.
void FETiedInterface::Serialize(DumpFile &ar)
{
	// store contact data
	FEContactInterface::Serialize(ar);

	// store contact surface data
	ms.Serialize(ar);
	ss.Serialize(ar);
}
