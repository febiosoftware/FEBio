// FERigidWallInterface.cpp: implementation of the FERigidWallInterface class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FERigidWallInterface.h"
#include "fem.h"
#include "FESolidSolver.h"
#include "FENNQuery.h"
#include "log.h"

///////////////////////////////////////////////////////////////////////////////
// FERigidWallInterface
///////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
//! constructor
FERigidWallInterface::FERigidWallInterface(FEM* pfem) : FEContactInterface(pfem), m_ss(&pfem->m_mesh)
{
	static int count = 1;
	m_ntype = FE_CONTACT_RIGIDWALL;

	m_mp = 0;

	m_nplc = -1;
	m_pplc = 0;
	m_nID = count++;
};


//-----------------------------------------------------------------------------
//! Initializes the rigid wall interface data

void FERigidWallInterface::Init()
{
	// create the surface
	m_ss.Init();

	// set the penalty load curve
	if (m_nplc >= 0) m_pplc = m_pfem->GetLoadCurve(m_nplc);

	// initialize rigid surface
	m_mp->Init();

	// project slave surface onto master surface
	ProjectSurface(m_ss);
}

//-----------------------------------------------------------------------------
//!  Projects the slave surface onto the master plane

void FERigidWallInterface::ProjectSurface(FEContactSurface& ss)
{
	vec3d r, q;

	// surface normal
	vec3d np;

	// loop over all slave nodes
	for (int i=0; i<m_ss.Nodes(); ++i)
	{
		// get the nodal position
		r = m_ss.Node(i).m_rt;

		// project this node onto the plane
		q = m_mp->Project(r);

		// get the local surface normal
		np = m_mp->Normal(q);

		// the slave normal is set to the master element normal
		m_ss.nu[i] = np;
	
		// calculate initial gap
		m_ss.gap[i] = -(np*(r - q)) + m_ss.off[i];
	}
}

//-----------------------------------------------------------------------------
//!  Updates rigid wall data

void FERigidWallInterface::Update()
{
	// project slave surface onto master surface
	ProjectSurface(m_ss);
}

//-----------------------------------------------------------------------------

void FERigidWallInterface::ContactForces(vector<double>& F)
{
	int j, k, m, n;
	int nseln;

	double *Gr, *Gs;

	// jacobian
	double detJ;

	vec3d dxr, dxs;
	vec3d* rt, *r0;
	double* w;

	// normal force
	double tn;

	// element contact force vector
	// note that this assumes that the element to which the integration node
	// connects is a four noded quadrilateral
	vector<double> fe(3);

	// the lm array for this force vector
	vector<int> lm(3);

	// the en array
	vector<int> en(1);

	vector<int> sLM;

	FEMesh& mesh = m_pfem->m_mesh;

	FESolidSolver* psolver = dynamic_cast<FESolidSolver*>(m_pfem->m_pStep->m_psolver);

	double eps = m_eps;
	
	// loop over all slave facets
	int ne = m_ss.Elements();
	for (j=0; j<ne; ++j)
	{
		// get the slave element
		FESurfaceElement& sel = m_ss.Element(j);
		mesh.UnpackElement(sel);

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
//			if (ss.Lm[m] >= 0)
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

				// get slave node normal force
				tn = m_ss.Lm[m] + eps*m_ss.gap[m];
				tn = MBRACKET(tn);

				// get the slave node normal
				vec3d& nu = m_ss.nu[m];

				// calculate force vector
				fe[0] = detJ*w[n]*tn*nu.x;
				fe[1] = detJ*w[n]*tn*nu.y;
				fe[2] = detJ*w[n]*tn*nu.z;
	
				// fill the lm array
				lm[0] = sLM[n*3  ];
				lm[1] = sLM[n*3+1];
				lm[2] = sLM[n*3+2];

				// fill the en array
				en[0] = sel.m_node[n];

				// assemble into global force vector
				psolver->AssembleResidual(en, lm, fe, F);
			}
		}
	}

}

//-----------------------------------------------------------------------------
//! This function calculates the stiffness contribution for the rigid wall
//! interface.
//! \todo I think there are a couple of stiffness terms missing in this formulation

void FERigidWallInterface::ContactStiffness()
{
	int j, k, l, n, m;
	int nseln, ndof;

	matrix ke;

	vector<int> lm(3);
	vector<int> en(1);

	double *Gr, *Gs, *w;
	vec3d *rt, *r0;

	double detJ, tn;
	vec3d dxr, dxs;
	double N[3];

	double gap, Lm;

	vector<int> sLM;

	FEMesh& mesh = m_pfem->m_mesh;

	FESolidSolver* psolver = dynamic_cast<FESolidSolver*>(m_pfem->m_pStep->m_psolver);

	// penalty value
	double eps = Penalty();

	// loop over all slave elements
	int ne = m_ss.Elements();
	for (j=0; j<ne; ++j)
	{
		FESurfaceElement& se = m_ss.Element(j);
		mesh.UnpackElement(se);

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
//			if (ss.Lm[m] >= 0)
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

				// slave gap
				gap = m_ss.gap[m];

				// lagrange multiplier
				Lm = m_ss.Lm[m];

				// get slave node normal force
				tn = m_ss.Lm[m] + eps*m_ss.gap[m];
				tn = MBRACKET(tn);

				// get the slave node normal
				vec3d& nu = m_ss.nu[m];

				// set up the N vector
				N[0] = nu.x;
				N[1] = nu.y;
				N[2] = nu.z;

				ndof = 3;

				// fill stiffness matrix
				// TODO: I don't think this is correct, since
				// if the rigid wall is not a plance, some terms
				// are probably missing
				ke.Create(ndof, ndof);
				for (k=0; k<ndof; ++k)
					for (l=0; l<ndof; ++l)
						ke[k][l] = w[n]*detJ*eps*HEAVYSIDE(Lm+eps*gap)*N[k]*N[l];

				// create lm array
				lm[0] = sLM[n*3  ];
				lm[1] = sLM[n*3+1];
				lm[2] = sLM[n*3+2];

				// create the en array
				en[0] = se.m_node[n];
						
				// assemble stiffness matrix
				psolver->AssembleStiffness(en, lm, ke);
			}
		}
	}
}

//-----------------------------------------------------------------------------

bool FERigidWallInterface::Augment(int naug)
{
	int i;
	double Lm;
	bool bconv = true;

	// penalty value
	double eps = Penalty();

	// calculate initial norms
	double normL0 = 0;
	for (i=0; i<m_ss.Nodes(); ++i) normL0 += m_ss.Lm[i]*m_ss.Lm[i];
	normL0 = sqrt(normL0);

	// update Lagrange multipliers and calculate current norms
	double normL1 = 0;
	double normgc = 0;
	int N = 0;
	for (i=0; i<m_ss.Nodes(); ++i)
	{
		// update Lagrange multipliers
		Lm = m_ss.Lm[i] + eps*m_ss.gap[i];
		Lm = MBRACKET(Lm);
		normL1 += Lm*Lm;
		if (m_ss.gap[i] > 0)
		{
			normgc += m_ss.gap[i]*m_ss.gap[i];
			++N;
		}
	}	
	if (N==0) N = 1;

	normL1 = sqrt(normL1);
	normgc = sqrt(normgc / N);

	// get the logfile
	Logfile& log = GetLogfile();

	// check convergence of constraints
	log.printf(" rigid wall interface # %d\n", m_nID);
	log.printf("                        CURRENT        REQUIRED\n");
	double pctn = 0;
	if (fabs(normL1) > 1e-10) pctn = fabs((normL1 - normL0)/normL1);
	log.printf("    normal force : %15le %15le\n", pctn, m_atol);
	log.printf("    gap function : %15le       ***\n", normgc);
		
	if (pctn >= m_atol)
	{
		bconv = false;
		for (i=0; i<m_ss.Nodes(); ++i)
		{
			// update Lagrange multipliers
			Lm = m_ss.Lm[i] + eps*m_ss.gap[i];
			m_ss.Lm[i] = MBRACKET(Lm);
		}	
	}

	return bconv;
}

//-----------------------------------------------------------------------------

void FERigidWallInterface::Serialize(Archive &ar)
{
	int k, n, mat;
	if (ar.IsSaving())
	{
		ar << m_eps;
		ar << m_atol;
		ar << m_nplc;

		FEContactSurface& s = m_ss;

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
		ar << s.Lm;
		ar << s.off;
		
		// plane data
		if (dynamic_cast<FEPlane*>(m_mp))
		{
			FEPlane* pp = dynamic_cast<FEPlane*>(m_mp);
			ar << FE_RIGID_PLANE;
			ar << pp->m_nplc;
			double* a = pp->GetEquation();
			ar << a[0] << a[1] << a[2] << a[3];
		}
		else if (dynamic_cast<FERigidSphere*>(m_mp))
		{
			FERigidSphere* ps = dynamic_cast<FERigidSphere*>(m_mp);
			ar << FE_RIGID_SPHERE;
			ar << ps->m_rc;
			ar << ps->m_R;
			ar << ps->m_nplc[0] << ps->m_nplc[1] << ps->m_nplc[2];
		}
	}
	else
	{
		ar >> m_eps;
		ar >> m_atol;
		ar >> m_nplc;

		if (m_nplc >= 0) m_pplc = m_pfem->GetLoadCurve(m_nplc);

		FEContactSurface& s = m_ss;

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

		ar >> s.gap;
		ar >> s.nu;
		ar >> s.rs;
		ar >> s.Lm;
		ar >> s.off;

		// plane data
		int ntype;
		ar >> ntype;
		switch (ntype)
		{
		case FE_RIGID_PLANE:
			{
				SetMasterSurface(new FEPlane(m_pfem));
				FEPlane& pl = dynamic_cast<FEPlane&>(*m_mp);
				ar >> pl.m_nplc;
				if (pl.m_nplc >= 0) pl.m_pplc = m_pfem->GetLoadCurve(pl.m_nplc);
				double* a = pl.GetEquation();
				ar >> a[0] >> a[1] >> a[2] >> a[4];
			}
			break;
		case FE_RIGID_SPHERE:
			{
				SetMasterSurface(new FERigidSphere(m_pfem));
				FERigidSphere& s = dynamic_cast<FERigidSphere&>(*m_mp);
				ar >> s.m_rc;
				ar >> s.m_R;
				ar >> s.m_nplc[0] >> s.m_nplc[1] >> s.m_nplc[2];
			}
			break;
		default:
			assert(false);
		}
	}
}
