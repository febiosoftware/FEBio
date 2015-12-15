// FERigidWallInterface.cpp: implementation of the FERigidWallInterface class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FERigidWallInterface.h"
#include "FECore/FENNQuery.h"
#include "FEElasticShellDomain.h"
#include "FEStiffnessMatrix.h"
#include "FECore/log.h"

//-----------------------------------------------------------------------------
// Define sliding interface parameters
BEGIN_PARAMETER_LIST(FERigidWallInterface, FEContactInterface)
	ADD_PARAMETER(m_blaugon, FE_PARAM_BOOL  , "laugon"      ); 
	ADD_PARAMETER(m_atol   , FE_PARAM_DOUBLE, "tolerance"   );
	ADD_PARAMETER(m_eps    , FE_PARAM_DOUBLE, "penalty"     );
END_PARAMETER_LIST();

///////////////////////////////////////////////////////////////////////////////
// FERigidWallSurface
///////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
//! Creates a surface for use with a sliding interface. All surface data
//! structures are allocated.
//! Note that it is assumed that the element array is already created
//! and initialized.

bool FERigidWallSurface::Init()
{
	int i, j, n;

	// always intialize base class first!
	if (FESurface::Init() == false) return false;

	// get the number of nodes
	int nn = Nodes();

	// allocate other surface data
	m_gap.assign(nn, 0);		// gap funtion
	m_nu.resize(nn);		// node normal 
	m_pme.assign(nn, static_cast<FESurfaceElement*>(0));		// penetrated master element
	m_rs.resize(nn);		// natural coords of projected slave node on master element
	m_rsp.resize(nn);
	m_Lm.assign(nn, 0);
	m_M.resize(nn);
	m_Lt.resize(nn);
	m_off.assign(nn, 0.0);
	m_eps.assign(nn, 1.0);

	// we calculate the gap offset values
	// This value is used to take the shell thickness into account
	// note that we force rigid shells to have zero thickness
	FEMesh& m = *m_pMesh;
	vector<double> tag(m.Nodes());
	zero(tag);
	for (int nd=0; nd<m.Domains(); ++nd)
	{
		FEElasticShellDomain* psd = dynamic_cast<FEElasticShellDomain*>(&m.Domain(nd));
		if (psd)
		{
			for (i=0; i<psd->Elements(); ++i)
			{
				FEShellElement& el = psd->Element(i);
				n = el.Nodes();
				for (j=0; j<n; ++j) tag[el.m_node[j]] = 0.5*el.m_h0[j];
			}
		}
	}
	for (i=0; i<nn; ++i) m_off[i] = tag[NodeIndex(i)];

	return true;
}


//-----------------------------------------------------------------------------
//! 

vec3d FERigidWallSurface::traction(int inode)
{
	vec3d t(0,0,0);
	FESurfaceElement* pe = m_pme[inode];
	if (pe)
	{
		FESurfaceElement& el = *pe;
		double Tn = m_Lm[inode];
		double T1 = m_Lt[inode][0];
		double T2 = m_Lt[inode][1];
		double r = m_rs[inode][0];
		double s = m_rs[inode][1];

		vec3d tn = m_nu[inode]*Tn, tt;
		vec3d e[2];
		ContraBaseVectors0(el, r, s, e);
		tt = e[0]*T1 + e[1]*T2;
		t = tn + tt;
	}

	return t;
}

//-----------------------------------------------------------------------------

void FERigidWallSurface::UpdateNormals()
{
	int i, j, jp1, jm1;
	int N = Nodes();
	int NE = Elements();
	for (i=0; i<N; ++i) m_nu[i] = vec3d(0,0,0);
	vec3d y[FEElement::MAX_NODES], e1, e2;

	for (i=0; i<NE; ++i)
	{
		FESurfaceElement& el = Element(i);
		int ne = el.Nodes();
		for (j=0; j<ne; ++j) y[j] = Node(el.m_lnode[j]).m_rt;

		for (j=0; j<ne; ++j)
		{
			jp1 = (j+1)%ne;
			jm1 = (j+ne-1)%ne;

			e1 = y[jp1] - y[j];
			e2 = y[jm1] - y[j];

			m_nu[el.m_lnode[j]] -= e1 ^ e2;						
		}
	}

	for (i=0; i<N; ++i) m_nu[i].unit();
}

//-----------------------------------------------------------------------------
void FERigidWallSurface::ShallowCopy(DumpStream& dmp, bool bsave)
{
	if (bsave)
	{
		dmp << m_Lm << m_gap << m_Lt;
	}
	else
	{
		dmp >> m_Lm >> m_gap >> m_Lt;
		zero(m_pme);
	}
}

//-----------------------------------------------------------------------------
void FERigidWallSurface::Serialize(DumpFile &ar)
{
	FESurface::Serialize(ar);
	if (ar.IsSaving())
	{
		ar << m_gap;
		ar << m_nu;
		ar << m_rs;
		ar << m_rsp;
		ar << m_Lm;
		ar << m_M;
		ar << m_Lt;
		ar << m_off;
		ar << m_eps;
	}
	else
	{
		ar >> m_gap;
		ar >> m_nu;
		ar >> m_rs;
		ar >> m_rsp;
		ar >> m_Lm;
		ar >> m_M;
		ar >> m_Lt;
		ar >> m_off;
		ar >> m_eps;
	}
}

//-----------------------------------------------------------------------------
void FERigidWallSurface::UnpackLM(FEElement& el, vector<int>& lm)
{
	int N = el.Nodes();
	lm.resize(N*3);
	for (int i=0; i<N; ++i)
	{
		int n = el.m_node[i];
		FENode& node = m_pMesh->Node(n);
		vector<int>& id = node.m_ID;

		lm[3*i  ] = id[DOF_X];
		lm[3*i+1] = id[DOF_Y];
		lm[3*i+2] = id[DOF_Z];
	}
}

///////////////////////////////////////////////////////////////////////////////
// FERigidWallInterface
///////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
//! constructor
FERigidWallInterface::FERigidWallInterface(FEModel* pfem) : FEContactInterface(pfem), m_ss(&pfem->GetMesh())
{
	static int count = 1;
	SetID(count++);

	m_mp = 0;

	m_eps = 0;
	m_atol = 0;
};

//-----------------------------------------------------------------------------
void FERigidWallInterface::ShallowCopy(DumpStream& dmp, bool bsave)
{
	m_ss.ShallowCopy(dmp, bsave);
}

//-----------------------------------------------------------------------------
//! Initializes the rigid wall interface data

bool FERigidWallInterface::Init()
{
	// create the surface
	if (m_ss.Init() == false) return false;

	// initialize rigid surface
	m_mp->Init();

	return true;
}

//-----------------------------------------------------------------------------
//! get the number of properties
int FERigidWallInterface::Properties()
{
	return 2;
}

//-----------------------------------------------------------------------------
//! get a specific property
FECoreBase* FERigidWallInterface::GetProperty(int i)
{
	switch (i)
	{
	case 0: if (m_mp == 0) SetMasterSurface(new FEPlane(GetFEModel())); break;
	case 1: if (m_mp == 0) SetMasterSurface(new FERigidSphere(GetFEModel())); break;
	}
	return m_mp;
}

//-----------------------------------------------------------------------------
//! find a property index ( returns <0 for error)
int FERigidWallInterface::FindPropertyIndex(const char* szname)
{
	if (strcmp(szname, "plane" ) == 0) return 0;
	if (strcmp(szname, "sphere") == 0) return 1;
	return -1;
}

//-----------------------------------------------------------------------------
//! set a property (returns false on error)
bool FERigidWallInterface::SetProperty(int i, FECoreBase* pm)
{
	return false;
}

//-----------------------------------------------------------------------------
//! build the matrix profile for use in the stiffness matrix
void FERigidWallInterface::BuildMatrixProfile(FEStiffnessMatrix& K)
{
	vector<int> lm(6);
	for (int j=0; j<m_ss.Nodes(); ++j)
	{
		if (m_ss.m_gap[j] >= 0)
		{
			lm[0] = m_ss.Node(j).m_ID[DOF_X];
			lm[1] = m_ss.Node(j).m_ID[DOF_Y];
			lm[2] = m_ss.Node(j).m_ID[DOF_Z];
			lm[3] = m_ss.Node(j).m_ID[DOF_RU];
			lm[4] = m_ss.Node(j).m_ID[DOF_RV];
			lm[5] = m_ss.Node(j).m_ID[DOF_RW];

			K.build_add(lm);
		}
	}
}

//-----------------------------------------------------------------------------
void FERigidWallInterface::Activate()
{
	// don't forget to call the base class
	FEContactInterface::Activate();

	// project slave surface onto master surface
	ProjectSurface(m_ss);
}

//-----------------------------------------------------------------------------
//!  Projects the slave surface onto the master plane

void FERigidWallInterface::ProjectSurface(FERigidWallSurface& ss)
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
		m_ss.m_nu[i] = np;
	
		// calculate initial gap
		m_ss.m_gap[i] = -(np*(r - q)) + m_ss.m_off[i];
	}
}

//-----------------------------------------------------------------------------
//!  Updates rigid wall data

void FERigidWallInterface::Update(int niter)
{
	// project slave surface onto master surface
	ProjectSurface(m_ss);
}

//-----------------------------------------------------------------------------

void FERigidWallInterface::ContactForces(FEGlobalVector& R)
{
	int j, k, m, n;
	int nseln;

	double *Gr, *Gs;

	// jacobian
	double detJ;

	vec3d dxr, dxs;
	vec3d rt[FEElement::MAX_NODES], r0[FEElement::MAX_NODES];
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

	// penalty value
	double pen = m_eps, eps;
	
	// loop over all slave facets
	int ne = m_ss.Elements();
	for (j=0; j<ne; ++j)
	{
		// get the slave element
		FESurfaceElement& sel = m_ss.Element(j);

		// get the element's LM vector
		m_ss.UnpackLM(sel, sLM);

		nseln = sel.Nodes();

		for (int i=0; i<nseln; ++i)
		{
			r0[i] = m_ss.GetMesh()->Node(sel.m_node[i]).m_r0;
			rt[i] = m_ss.GetMesh()->Node(sel.m_node[i]).m_rt;
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
				eps = pen*m_ss.m_eps[m];
				tn = m_ss.m_Lm[m] + eps*m_ss.m_gap[m];
				tn = MBRACKET(tn);

				// get the slave node normal
				vec3d& nu = m_ss.m_nu[m];

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
				R.Assemble(en, lm, fe);
			}
		}
	}

}

//-----------------------------------------------------------------------------
//! This function calculates the stiffness contribution for the rigid wall
//! interface.
//! \todo I think there are a couple of stiffness terms missing in this formulation

void FERigidWallInterface::ContactStiffness(FESolver* psolver)
{
	int j, k, l, n, m;
	int nseln, ndof;

	matrix ke;

	vector<int> lm(3);
	vector<int> en(1);

	double *Gr, *Gs, *w;
	vec3d rt[FEElement::MAX_NODES], r0[FEElement::MAX_NODES];

	double detJ, tn;
	vec3d dxr, dxs;
	double N[3];

	double gap, Lm;

	vector<int> sLM;

	// penalty value
	double pen = m_eps, eps;

	// loop over all slave elements
	int ne = m_ss.Elements();
	for (j=0; j<ne; ++j)
	{
		FESurfaceElement& se = m_ss.Element(j);

		// get the element's LM vector
		m_ss.UnpackLM(se, sLM);

		nseln = se.Nodes();

		for (int i=0; i<nseln; ++i)
		{
			r0[i] = m_ss.GetMesh()->Node(se.m_node[i]).m_r0;
			rt[i] = m_ss.GetMesh()->Node(se.m_node[i]).m_rt;
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
				gap = m_ss.m_gap[m];

				// lagrange multiplier
				Lm = m_ss.m_Lm[m];

				// get slave node normal force
				eps = pen*m_ss.m_eps[m];

				tn = m_ss.m_Lm[m] + eps*m_ss.m_gap[m];
				tn = MBRACKET(tn);

				// get the slave node normal
				vec3d& nu = m_ss.m_nu[m];

				// set up the N vector
				N[0] = nu.x;
				N[1] = nu.y;
				N[2] = nu.z;

				ndof = 3;

				// fill stiffness matrix
				// TODO: I don't think this is correct, since
				// if the rigid wall is not a plance, some terms
				// are probably missing
				ke.resize(ndof, ndof);
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
	// make sure we need to augment
	if (!m_blaugon) return true;

	int i;
	double Lm;
	bool bconv = true;

	// penalty value
	double pen = m_eps, eps;

	// calculate initial norms
	double normL0 = 0;
	for (i=0; i<m_ss.Nodes(); ++i) normL0 += m_ss.m_Lm[i]*m_ss.m_Lm[i];
	normL0 = sqrt(normL0);

	// update Lagrange multipliers and calculate current norms
	double normL1 = 0;
	double normgc = 0;
	int N = 0;
	for (i=0; i<m_ss.Nodes(); ++i)
	{
		// update Lagrange multipliers
		eps = pen*m_ss.m_eps[i];

		Lm = m_ss.m_Lm[i] + eps*m_ss.m_gap[i];
		Lm = MBRACKET(Lm);
		normL1 += Lm*Lm;
		if (m_ss.m_gap[i] > 0)
		{
			normgc += m_ss.m_gap[i]*m_ss.m_gap[i];
			++N;
		}
	}	
	if (N==0) N = 1;

	normL1 = sqrt(normL1);
	normgc = sqrt(normgc / N);

	// check convergence of constraints
	felog.printf(" rigid wall interface # %d\n", GetID());
	felog.printf("                        CURRENT        REQUIRED\n");
	double pctn = 0;
	if (fabs(normL1) > 1e-10) pctn = fabs((normL1 - normL0)/normL1);
	felog.printf("    normal force : %15le %15le\n", pctn, m_atol);
	felog.printf("    gap function : %15le       ***\n", normgc);
		
	if (pctn >= m_atol)
	{
		bconv = false;
		for (i=0; i<m_ss.Nodes(); ++i)
		{
			// update Lagrange multipliers
			eps = pen*m_ss.m_eps[i];

			Lm = m_ss.m_Lm[i] + eps*m_ss.m_gap[i];
			m_ss.m_Lm[i] = MBRACKET(Lm);
		}	
	}

	return bconv;
}

//-----------------------------------------------------------------------------

void FERigidWallInterface::Serialize(DumpFile &ar)
{
	FEContactInterface::Serialize(ar);

	if (ar.IsSaving())
	{
		ar << m_eps;
		ar << m_atol;

		m_ss.Serialize(ar);
		
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
		}
	}
	else
	{
		ar >> m_eps;
		ar >> m_atol;

		m_ss.Serialize(ar);

		// plane data
		int ntype;
		ar >> ntype;
		switch (ntype)
		{
		case FE_RIGID_PLANE:
			{
				SetMasterSurface(new FEPlane(GetFEModel()));
				FEPlane& pl = dynamic_cast<FEPlane&>(*m_mp);
				ar >> pl.m_nplc;
				if (pl.m_nplc >= 0) pl.m_pplc = GetFEModel()->GetLoadCurve(pl.m_nplc);
				double* a = pl.GetEquation();
				ar >> a[0] >> a[1] >> a[2] >> a[3];
			}
			break;
		case FE_RIGID_SPHERE:
			{
				SetMasterSurface(new FERigidSphere(GetFEModel()));
				FERigidSphere& s = dynamic_cast<FERigidSphere&>(*m_mp);
				ar >> s.m_rc;
				ar >> s.m_R;
			}
			break;
		default:
			assert(false);
		}
	}
}
