// FETiedInterface.cpp: implementation of the FETiedInterface class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FETiedInterface.h"
#include "FEStiffnessMatrix.h"
#include "FECore/FEModel.h"
#include "FECore/FEClosestPointProjection.h"
#include "FECore/log.h"

//-----------------------------------------------------------------------------
// Define sliding interface parameters
BEGIN_PARAMETER_LIST(FETiedInterface, FEContactInterface)
	ADD_PARAMETER(m_blaugon , FE_PARAM_BOOL  , "laugon"          ); 
	ADD_PARAMETER(m_atol    , FE_PARAM_DOUBLE, "tolerance"       );
	ADD_PARAMETER(m_eps     , FE_PARAM_DOUBLE, "penalty"         );
	ADD_PARAMETER(m_naugmin , FE_PARAM_INT   , "minaug"          );
	ADD_PARAMETER(m_naugmax , FE_PARAM_INT   , "maxaug"          );
	ADD_PARAMETER(m_stol    , FE_PARAM_DOUBLE, "search_tolerance");
	ADD_PARAMETER(m_boffset , FE_PARAM_BOOL  , "offset_shells"   );
	ADD_PARAMETER(m_Dmax    , FE_PARAM_DOUBLE, "max_distance"    );
	ADD_PARAMETER(m_bspecial, FE_PARAM_BOOL  , "special"         );
	ADD_PARAMETER(m_breloc  , FE_PARAM_BOOL  , "node_reloc"      );
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! Constructor. Initialize default values.
FETiedInterface::FETiedInterface(FEModel* pfem) : FEContactInterface(pfem), ss(pfem), ms(pfem)
{
	static int count = 1;
	SetID(count++);

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
	m_boffset = false;
	m_Dmax = 0.0;
	m_bspecial = true;
	m_breloc = false;
}

//-----------------------------------------------------------------------------
//! Initialization. This function intializes the surfaces data and projects the
//! slave surface onto the master surface.
//! 
bool FETiedInterface::Init()
{
	// set surface options
	ss.SetShellOffset(m_boffset);

	// create the surfaces
	if (ss.Init() == false) return false;
	if (ms.Init() == false) return false;

	return true;
}

//-----------------------------------------------------------------------------
//! build the matrix profile for use in the stiffness matrix
void FETiedInterface::BuildMatrixProfile(FEGlobalMatrix& K)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// get the DOFS
	const int dof_X = fem.GetDOFIndex("x");
	const int dof_Y = fem.GetDOFIndex("y");
	const int dof_Z = fem.GetDOFIndex("z");
	const int dof_RU = fem.GetDOFIndex("Ru");
	const int dof_RV = fem.GetDOFIndex("Rv");
	const int dof_RW = fem.GetDOFIndex("Rw");

	const int LMSIZE = 6*(FEElement::MAX_NODES+1);
	vector<int> lm(LMSIZE);

	for (int j=0; j<ss.Nodes(); ++j)
	{
		FESurfaceElement* pe = ss.m_pme[j];
		if (pe != 0)
		{
			FESurfaceElement& me = *pe;
			int* en = &me.m_node[0];

			int n = me.Nodes();
			lm.assign(LMSIZE, -1);

			lm[0] = ss.Node(j).m_ID[dof_X];
			lm[1] = ss.Node(j).m_ID[dof_Y];
			lm[2] = ss.Node(j).m_ID[dof_Z];
			lm[3] = ss.Node(j).m_ID[dof_RU];
			lm[4] = ss.Node(j).m_ID[dof_RV];
			lm[5] = ss.Node(j).m_ID[dof_RW];

			for (int k=0; k<n; ++k)
			{
				vector<int>& id = mesh.Node(en[k]).m_ID;
				lm[6*(k+1)  ] = id[dof_X];
				lm[6*(k+1)+1] = id[dof_Y];
				lm[6*(k+1)+2] = id[dof_Z];
				lm[6*(k+1)+3] = id[dof_RU];
				lm[6*(k+1)+4] = id[dof_RV];
				lm[6*(k+1)+5] = id[dof_RW];
			}

			K.build_add(lm);
		}
	}
}

//-----------------------------------------------------------------------------
//! Interface activation
void FETiedInterface::Activate()
{
	// Don't forget to call base member!
	FEContactInterface::Activate();

	// project slave surface onto master surface
	ProjectSurface(ss, ms, m_breloc);
}

//-----------------------------------------------------------------------------
void FETiedInterface::ShallowCopy(DumpStream& dmp, bool bsave)
{
	ss.ShallowCopy(dmp, bsave);
	ms.ShallowCopy(dmp, bsave);
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
		FESurfaceElement* pme = ss.m_pme[i];
		if (pme)
		{
			// get the current slave nodal position
			vec3d rt = ss.Node(i).m_rt;

			// get the natural coordinates of the slave projection
			// onto the master element
			double r = ss.m_rs[i][0];
			double s = ss.m_rs[i][1];

			// get the nodal coordinates
			int ne = pme->Nodes();
			vec3d y[FEElement::MAX_NODES];
			for (int l=0; l<ne; ++l) y[l] = mesh.Node( pme->m_node[l] ).m_rt;

			// calculate the slave node projection
			vec3d q = pme->eval(y, r, s);

			// calculate the master normal
			vec3d nu = ss.SurfaceNormal(*pme, r, s);

			// calculate the gap function
			// (taking possible offset into account)
			ss.m_gap[i] = (rt - q) - nu*ss.m_off[i];

			// calculate force
			ss.m_Tc[i] = ss.m_Lm[i] + ss.m_gap[i]*m_eps;
		}
	}
}

//-----------------------------------------------------------------------------
//! project surface

void FETiedInterface::ProjectSurface(FETiedContactSurface& ss, FETiedContactSurface& ms, bool bmove)
{
	// closest point projection method
	FEClosestPointProjection cpp(ms);
	cpp.SetTolerance(m_stol);
	cpp.HandleSpecialCases(m_bspecial);
	cpp.Init();

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
		FESurfaceElement* pme = cpp.Project(x, q, rs);
		if (pme)
		{
			// make sure we are within the max distance
			double D = (x - q).norm();
			if ((m_Dmax == 0.0) || (D <= m_Dmax))
			{
				// store the master element
				ss.m_pme[i] = pme;

				// store the natural coordinates of the projection on the master element
				ss.m_rs[i] = rs;

				// calculate the master normal
				vec3d nu = ms.SurfaceNormal(*pme, rs[0], rs[1]);

				// calculate gap
				ss.m_gap[i] = (x - q) - nu*ss.m_off[i];

				// move the node if necessary
				if (bmove && (ss.m_gap[i].norm()>0))
				{
					node.m_r0 = node.m_rt = q + nu*ss.m_off[i];
					ss.m_gap[i] = vec3d(0,0,0);
				}

				// calculate force
				ss.m_Tc[i] = ss.m_Lm[i] + ss.m_gap[i]*m_eps;
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! This function calculates the contact forces for a tied interface.

void FETiedInterface::ContactForces(FEGlobalVector& R)
{
	// shape function values
	double N[FEElement::MAX_NODES];

	// element contact force vector
	vector<double> fe;

	// the lm array for this force vector
	vector<int> lm;

	// the en array
	vector<int> en;

	vector<int> sLM;
	vector<int> mLM;

	// loop over all slave facets
	const int NE = ss.Elements();
	for (int j=0; j<NE; ++j)
	{
		// get the slave element
		FESurfaceElement& sel = ss.Element(j);

		// get the element's LM vector
		ss.UnpackLM(sel, sLM);

		int nseln = sel.Nodes();

		double* w = sel.GaussWeights();

		// loop over slave element nodes (which are the integration points as well)
		for (int n=0; n<nseln; ++n)
		{
			int m = sel.m_lnode[n];

			// see if this node's constraint is active
			// that is, if it has a master element associated with it
			// TODO: is this a good way to test for an active constraint
			// The rigid wall criteria seems to work much better.
			if (ss.m_pme[m] != 0)
			{
				// calculate jacobian
				double detJ = ss.jac0(sel, n);

				// get slave node contact force
				vec3d tc = ss.m_Lm[m] + ss.m_gap[m]*m_eps;

				// get the master element
				FESurfaceElement& mel = *ss.m_pme[m];
				ms.UnpackLM(mel, mLM);

				int nmeln = mel.Nodes();

				// isoparametric coordinates of the projected slave node
				// onto the master element
				double r = ss.m_rs[m][0];
				double s = ss.m_rs[m][1];

				// get the master shape function values at this slave node
				mel.shape_fnc(N, r, s);

				// calculate force vector
				fe.resize(3*(nmeln+1));
				fe[0] = -detJ*w[n]*tc.x;
				fe[1] = -detJ*w[n]*tc.y;
				fe[2] = -detJ*w[n]*tc.z;
				for (int l=0; l<nmeln; ++l)
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

				for (int l=0; l<nmeln; ++l)
				{
					lm[3*(l+1)  ] = mLM[l*3  ];
					lm[3*(l+1)+1] = mLM[l*3+1];
					lm[3*(l+1)+2] = mLM[l*3+2];
				}

				// fill the en array
				en.resize(nmeln+1);
				en[0] = sel.m_node[n];
				for (int l=0; l<nmeln; ++l) en[l+1] = mel.m_node[l];

				// assemble into global force vector
				R.Assemble(en, lm, fe);
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! Calculate the stiffness matrix contribution.
void FETiedInterface::ContactStiffness(FESolver* psolver)
{
	vector<int> sLM, mLM, lm, en;
	const int MN = FEElement::MAX_NODES;
	matrix ke;

	// shape functions
	double H[MN];

	// loop over all slave elements
	const int NE = ss.Elements();
	for (int i=0; i<NE; ++i)
	{
		// get the next element
		FESurfaceElement& se = ss.Element(i);
		int nseln = se.Nodes();

		// get the element's LM vector
		ss.UnpackLM(se, sLM);

		double* w = se.GaussWeights();

		// loop over all integration points (that is nodes)
		for (int n=0; n<nseln; ++n)
		{
			int m = se.m_lnode[n];

			// get the master element
			FESurfaceElement* pme = ss.m_pme[m];
			if (pme)
			{
				// get the master element
				FESurfaceElement& me = *pme;
				int nmeln = me.Nodes();
				ms.UnpackLM(me, mLM);

				// calculate jacobian
				double detJ = ss.jac0(se, n);

				// slave node natural coordinates in master element
				double r = ss.m_rs[m][0];
				double s = ss.m_rs[m][1];

				// get the master shape function values at this slave node
				me.shape_fnc(H, r, s);

				// number of degrees of freedom
				int ndof = 3*(1 + nmeln);

				// fill stiffness matrix
				ke.resize(ndof, ndof); ke.zero();
				ke[0][0] = w[n]*detJ*m_eps;
				ke[1][1] = w[n]*detJ*m_eps;
				ke[2][2] = w[n]*detJ*m_eps;
				for (int k=0; k<nmeln; ++k)
				{
					ke[0][3+3*k  ] = -w[n]*detJ*m_eps*H[k];
					ke[1][3+3*k+1] = -w[n]*detJ*m_eps*H[k];
					ke[2][3+3*k+2] = -w[n]*detJ*m_eps*H[k];

					ke[3+3*k  ][0] = -w[n]*detJ*m_eps*H[k];
					ke[3+3*k+1][1] = -w[n]*detJ*m_eps*H[k];
					ke[3+3*k+2][2] = -w[n]*detJ*m_eps*H[k];
				}
				for (int k=0; k<nmeln; ++k)
					for (int l=0; l<nmeln; ++l)
					{
						ke[3+3*k  ][3+3*l  ] = w[n]*detJ*m_eps*H[k]*H[l];
						ke[3+3*k+1][3+3*l+1] = w[n]*detJ*m_eps*H[k]*H[l];
						ke[3+3*k+2][3+3*l+2] = w[n]*detJ*m_eps*H[k]*H[l];
					}

				// create lm array
				lm.resize(3*(1+nmeln));
				lm[0] = sLM[n*3  ];
				lm[1] = sLM[n*3+1];
				lm[2] = sLM[n*3+2];

				for (int k=0; k<nmeln; ++k)
				{
					lm[3*(k+1)  ] = mLM[k*3  ];
					lm[3*(k+1)+1] = mLM[k*3+1];
					lm[3*(k+1)+2] = mLM[k*3+2];
				}

				// create the en array
				en.resize(nmeln+1);
				en[0] = se.m_node[n];
				for (int k=0; k<nmeln; ++k) en[k+1] = me.m_node[k];
						
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
	felog.printf(" tied interface # %d\n", GetID());
	felog.printf("                        CURRENT        REQUIRED\n");
	double pctn = 0;
	if (fabs(normL1) > 1e-10) pctn = fabs((normL1 - normL0)/normL1);
	felog.printf("    normal force : %15le %15le\n", pctn, m_atol);
	felog.printf("    gap function : %15le       ***\n", normgc);
		
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
