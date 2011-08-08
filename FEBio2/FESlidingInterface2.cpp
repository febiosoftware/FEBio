#include "stdafx.h"
#include "fem.h"
#include "FESlidingInterface2.h"
#include "FESolidSolver.h"
#include "FEPoroElastic.h"
#include "FEBiphasic.h"
#include "log.h"

//-----------------------------------------------------------------------------
// Register the class with the framework
REGISTER_FEBIO_CLASS(FESlidingInterface2, FEContactInterface, "sliding2");

//-----------------------------------------------------------------------------
// Define sliding interface parameters
BEGIN_PARAMETER_LIST(FESlidingInterface2, FEContactInterface)
	ADD_PARAMETER(m_blaugon  , FE_PARAM_BOOL  , "laugon"             );
	ADD_PARAMETER(m_atol     , FE_PARAM_DOUBLE, "tolerance"          );
	ADD_PARAMETER(m_gtol     , FE_PARAM_DOUBLE, "gaptol"             );
	ADD_PARAMETER(m_ptol     , FE_PARAM_DOUBLE, "ptol"               );
	ADD_PARAMETER(m_epsn     , FE_PARAM_DOUBLE, "penalty"            );
	ADD_PARAMETER(m_bautopen , FE_PARAM_BOOL  , "auto_penalty"       );
	ADD_PARAMETER(m_btwo_pass, FE_PARAM_BOOL  , "two_pass"           );
	ADD_PARAMETER(m_knmult   , FE_PARAM_INT   , "knmult"             );
	ADD_PARAMETER(m_stol     , FE_PARAM_DOUBLE, "search_tol"         );
	ADD_PARAMETER(m_epsp     , FE_PARAM_DOUBLE, "pressure_penalty"   );
	ADD_PARAMETER(m_bsymm    , FE_PARAM_BOOL  , "symmetric_stiffness");
	ADD_PARAMETER(m_srad     , FE_PARAM_DOUBLE, "search_radius"      );
	ADD_PARAMETER(m_nsegup   , FE_PARAM_INT   , "seg_up"             );
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
// FESlidingSurface2
//-----------------------------------------------------------------------------

FESlidingSurface2::FESlidingSurface2(FEModel* pfem) : FEContactSurface(&pfem->m_mesh)
{ 
	m_bporo = false;
	m_pfem = pfem; 
}

//-----------------------------------------------------------------------------
void FESlidingSurface2::Init()
{
	// initialize surface data first
	FEContactSurface::Init();

	// count how many integration points we have
	int nint = 0, i;
	for (i=0; i<Elements(); ++i)
	{
		FESurfaceElement& el = Element(i);
		nint += el.GaussPoints();
	}

	// allocate data structures
	m_gap.assign(nint, 0.0);
	m_nu.resize(nint);
	m_rs.resize(nint);
	m_Lmd.assign(nint, 0.0);
	m_Lmp.assign(nint, 0.0);
	m_pme.assign(nint, static_cast<FESurfaceElement*>(0));
	m_epsn.assign(nint, 1.0);
	m_epsp.assign(nint, 1.0);
	m_Ln.assign(nint, 0.0);

	m_nn.assign(Nodes(), 0);

	m_nei.resize(Elements());
	nint = 0;
	for (i=0; i<Elements(); ++i)
	{
		m_nei[i] = nint;
		nint += Element(i).GaussPoints();
	}

	// set intial values
	zero(m_nu);

	// determine biphasic status
	i = 0;
	while (!m_bporo && (i<Elements())) {
		// get the surface element
		FESurfaceElement& se = Element(i);

		// get the solid element this surface element belongs to
		FESolidElement* pe = dynamic_cast<FESolidElement*>(m_pMesh->FindElementFromID(se.m_nelem));
		if (pe)
		{
			// get the material
			FEMaterial* pm = dynamic_cast<FEMaterial*>(m_pfem->GetMaterial(pe->GetMatID()));
			
			// see if this is a poro-elastic element
			FEPoroElastic* poro = dynamic_cast<FEPoroElastic*>(pm);
			FEBiphasic* biph = dynamic_cast<FEBiphasic*> (pm);
			if (poro || biph) m_bporo = true;
		}
		++i;
	}
		
	// allocate biphasic stuff
	if (m_bporo)
	{
		m_pg.assign(nint, 0);
	}
}

//-----------------------------------------------------------------------------
void FESlidingSurface2::ShallowCopy(FESlidingSurface2 &s)
{
	m_Lmd = s.m_Lmd;
	m_gap = s.m_gap;
	m_Ln = s.m_Ln;
	zero(m_pme);
	m_bporo = s.m_bporo;

	if (m_bporo)
	{
		m_pg  = s.m_pg;
		m_Lmp = s.m_Lmp;
	}
}

//-----------------------------------------------------------------------------
//! This function calculates the node normal. Due to the piecewise continuity
//! of the surface elements this normal is not uniquely defined so in order to
//! obtain a unique normal the normal is averaged for each node over all the 
//! element normals at the node

void FESlidingSurface2::UpdateNodeNormals()
{
	int N = Nodes(), i, j, ne, jp1, jm1;
	vec3d y[4], n;

	// zero nodal normals
	zero(m_nn);

	// loop over all elements
	for (i=0; i<Elements(); ++i)
	{
		FESurfaceElement& el = Element(i);
		ne = el.Nodes();

		// get the nodal coordinates
		for (j=0; j<ne; ++j) y[j] = Node(el.m_lnode[j]).m_rt;

		// calculate the normals
		for (j=0; j<ne; ++j)
		{
			jp1 = (j+1)%ne;
			jm1 = (j+ne-1)%ne;
			n = (y[jp1] - y[j]) ^ (y[jm1] - y[j]);
			m_nn[el.m_lnode[j]] += n;
		}
	}

	// normalize all vectors
	for (i=0; i<N; ++i) m_nn[i].unit();
}

//-----------------------------------------------------------------------------
void FESlidingSurface2::Serialize(DumpFile& ar)
{
	// We need to store the m_bporo flag first 
	// since we need it before we initialize the surface data
	if (ar.IsSaving())
	{
		ar << m_bporo;
	}
	else
	{
		ar >> m_bporo;
	}

	// Next, we can serialize the base-class data
	FEContactSurface::Serialize(ar);

	// And finally, we serialize the surface data
	if (ar.IsSaving())
	{
		ar << m_gap;
		ar << m_nu;
		ar << m_rs;
		ar << m_Lmd;
		ar << m_Lmp;
		ar << m_nei;
		ar << m_epsn;
		ar << m_epsp;
		ar << m_nn;
		ar << m_pg;
		ar << m_Ln;

		int ne = (int) m_pme.size();
		ar << ne;
		for (int i=0; i<ne; ++i)
		{
			FESurfaceElement* pe = m_pme[i];
			if (pe) ar << pe->m_lid; else ar << -1;
		}
	}
	else
	{
		ar >> m_gap;
		ar >> m_nu;
		ar >> m_rs;
		ar >> m_Lmd;
		ar >> m_Lmp;
		ar >> m_nei;
		ar >> m_epsn;
		ar >> m_epsp;
		ar >> m_nn;
		ar >> m_pg;
		ar >> m_Ln;

		assert(m_pSibling);

		int ne = 0, id;
		ar >> ne;
		assert(ne == m_pme.size());
		for (int i=0; i<ne; ++i)
		{
			ar >> id;
			if (id < 0) m_pme[i] = 0; 
			else 
			{
				m_pme[i] = &m_pSibling->Element(id);
				assert(m_pme[i]->m_lid == id);
			}
		}
	}
}

//-----------------------------------------------------------------------------
// FESlidingInterface2
//-----------------------------------------------------------------------------

FESlidingInterface2::FESlidingInterface2(FEModel* pfem) : FEContactInterface(pfem), m_ss(pfem), m_ms(pfem)
{
	m_ntype = FE_CONTACT_SLIDING2;
	static int count = 1;
	m_nID = count++;

	// initial values
	m_knmult = 1;
	m_atol = 0.1;
	m_epsn = 1;
	m_epsp = 1;
	m_btwo_pass = false;
	m_stol = 0.01;
	m_bsymm = true;
	m_srad = 0.1;
	m_gtol = -1;	// we use augmentation tolerance by default
	m_ptol = -1;	// we use augmentation tolerance by default
	m_bautopen = false;
	m_nsegup = 0;

	m_naugmin = 0;
	m_naugmax = 10;

	m_ss.SetSibling(&m_ms);
	m_ms.SetSibling(&m_ss);
}

//-----------------------------------------------------------------------------

FESlidingInterface2::~FESlidingInterface2()
{
}

//-----------------------------------------------------------------------------
void FESlidingInterface2::Init()
{
	// initialize surface data
	m_ss.Init();
	m_ms.Init();

	FEM& fem = dynamic_cast<FEM&>(*m_pfem);
	
	bool bporo = (m_ss.m_bporo || m_ms.m_bporo);

	// this contact implementation requires a non-symmetric stiffness matrix
	// so inform the FEM class
	if (!m_bsymm) 
	{
		// request a non-symmetric stiffness matrix
		fem.SetSymmetryFlag(false);

		// make sure we are using full-Newton
		if (bporo && (fem.m_pStep->m_psolver->m_bfgs.m_maxups != 0))
		{
			fem.m_pStep->m_psolver->m_bfgs.m_maxups = 0;
			clog.printbox("WARNING", "The non-symmetric biphasic contact algorithm does not work with BFGS yet.\nThe full-Newton method will be used instead.");
		}
	}

	// calculate the penalty
	if (m_bautopen) 
	{
		CalcAutoPenalty(m_ss);
		CalcAutoPenalty(m_ms);
		CalcAutoPressurePenalty(m_ss);
		CalcAutoPressurePenalty(m_ms);
	}

	// update sliding interface data
	Update();
}

//-----------------------------------------------------------------------------
void FESlidingInterface2::CalcAutoPenalty(FESlidingSurface2& s)
{
	// get the mesh
	FEMesh& m = m_pfem->m_mesh;

	// loop over all surface elements
	int ni = 0;
	for (int i=0; i<s.Elements(); ++i)
	{
		// get the surface element
		FESurfaceElement& el = s.Element(i);

		// find the element this face belongs to
		FEElement* pe = m.FindElementFromID(el.m_nelem);
		assert(pe);

		// get the area of the surface element
		double A = s.FaceArea(el);

		// get the volume of the volume element
		double V = m.ElementVolume(*pe);

		// calculate a modulus
		double K = BulkModulus(el, s);

		// calculate penalty
		double eps = K*A/V;

		// assign to integation points of surface element
		int nint = el.GaussPoints();
		for (int j=0; j<nint; ++j, ++ni) s.m_epsn[ni] = eps;
	}
}

//-----------------------------------------------------------------------------
void FESlidingInterface2::CalcAutoPressurePenalty(FESlidingSurface2& s)
{
	// get the mesh
	FEMesh& m = m_pfem->m_mesh;

	// loop over all surface elements
	int ni = 0;
	for (int i=0; i<s.Elements(); ++i)
	{
		// get the surface element
		FESurfaceElement& el = s.Element(i);

		// find the element this face belongs to
		FEElement* pe = m.FindElementFromID(el.m_nelem);
		assert(pe);

		// get the area of the surface element
		double A = s.FaceArea(el);

		// get the volume of the volume element
		double V = m.ElementVolume(*pe);

		// calculate a modulus
		double k = AutoPressurePenalty(el, s);

		// calculate penalty
		double eps = k*A/V;

		// assign to integation points of surface element
		int nint = el.GaussPoints();
		for (int j=0; j<nint; ++j, ++ni) s.m_epsp[ni] = eps;
	}
}

//-----------------------------------------------------------------------------

double FESlidingInterface2::AutoPressurePenalty(FESurfaceElement& el, FESlidingSurface2& s)
{
	// get the mesh
	FEMesh& m = m_pfem->m_mesh;

	double eps = 0;

	// get the solid element this surface element belongs to
	FESolidElement* pe = dynamic_cast<FESolidElement*>(m.FindElementFromID(el.m_nelem));
	if (pe)
	{
		// get the material
		FEMaterial* pm = dynamic_cast<FEMaterial*>(m_pfem->GetMaterial(pe->GetMatID()));

		// see if this is a poro-elastic element
		FEPoroElastic* poro = dynamic_cast<FEPoroElastic*>(pm);
		FEBiphasic* biph = dynamic_cast<FEBiphasic*> (pm);
		if (poro || biph)
		{
			// get a material point
			FEMaterialPoint& mp = *pe->m_State[0];
			FEElasticMaterialPoint& ept = *(mp.ExtractData<FEElasticMaterialPoint>());

			// setup the material point
			ept.F = mat3dd(1.0);
			ept.J = 1;
			ept.s.zero();

			// if this is a poroelastic element, then get the permeability tensor
			FEPoroElasticMaterialPoint& pt = *(mp.ExtractData<FEPoroElasticMaterialPoint>());
			pt.m_p = 0;
			pt.m_w = vec3d(0,0,0);
					
			double K[3][3];
			if (poro) poro->Permeability(K, mp);
			else biph->Permeability(K, mp);

			eps = (K[0][0] + K[1][1] + K[2][2])/3;
		}
	}

	return eps;
}

//-----------------------------------------------------------------------------
void FESlidingInterface2::ProjectSurface(FESlidingSurface2& ss, FESlidingSurface2& ms, bool bupseg)
{
	FEMesh& mesh = m_pfem->m_mesh;
	FESurfaceElement* pme;
	vec3d r, nu;
	double rs[2];
	double Ln;

	double ps[4], p1;

	double R = m_srad*mesh.GetBoundingBox().radius();

	// loop over all integration points
	int n = 0;
	for (int i=0; i<ss.Elements(); ++i)
	{
		FESurfaceElement& el = ss.Element(i);
		bool sporo = PoroStatus(mesh, el);

		int ne = el.Nodes();
		int nint = el.GaussPoints();

		// get the nodal pressures
		if (sporo)
		{
			for (int j=0; j<ne; ++j) ps[j] = mesh.Node(el.m_node[j]).m_pt;
		}

		for (int j=0; j<nint; ++j, ++n)
		{
			// calculate the global position of the integration point
			r = ss.Local2Global(el, j);

			// get the pressure at the integration point
			if (sporo) p1 = el.eval(ps, j);

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
				else
				{
					pme = 0;
				}
			}

			// find the intersection point with the master surface
			if (pme == 0 && bupseg) pme = ms.FindIntersection(r, nu, rs, m_stol);

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

				double eps = m_epsn*ss.m_epsn[n];

				Ln = ss.m_Lmd[n] + eps*g;

				ss.m_gap[n] = (g <= R? g : 0);

				if ((Ln >= 0) && (g <= R))
				{

					// calculate the pressure gap function
					bool mporo = PoroStatus(mesh, *pme);
					if (sporo && mporo) {
						double pm[4];
						for (int k=0; k<pme->Nodes(); ++k) pm[k] = mesh.Node(pme->m_node[k]).m_pt;
						double p2 = pme->eval(pm, rs[0], rs[1]);
						ss.m_pg[n] = p1 - p2;
					}
				}
				else
				{
//					ss.m_Lmd[n] = 0;
//					ss.m_gap[n] = 0;
					ss.m_pme[n] = 0;
					if (sporo) {
//						ss.m_Lmp[n] = 0;
//						ss.m_pg[n] = 0;
					}
				}
			}
			else
			{
				// the node is not in contact
				ss.m_Lmd[n] = 0;
				ss.m_gap[n] = 0;
				if (sporo) {
					ss.m_Lmp[n] = 0;
					ss.m_pg[n] = 0;
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------

void FESlidingInterface2::Update()
{	
	int i, n, id, np;
	double rs[2];

	double R = m_srad*m_pfem->m_mesh.GetBoundingBox().radius();

	static int naug = 0;
	static int biter = 0;

	FEM& fem = dynamic_cast<FEM&>(*m_pfem);
	
	// get the iteration number
	// we need this number to see if we can do segment updates or not
	// also reset number of iterations after each augmentation
	if (fem.m_pStep->m_psolver->m_niter == 0) {
		biter = 0;
		naug = fem.m_pStep->m_psolver->m_naug;
	} else if (fem.m_pStep->m_psolver->m_naug > naug) {
		biter = fem.m_pStep->m_psolver->m_niter;
		naug = fem.m_pStep->m_psolver->m_naug;
	}
	int niter = fem.m_pStep->m_psolver->m_niter - biter;
	bool bupseg = ((m_nsegup == 0)? true : (niter <= m_nsegup));
	// get the logfile
//	Logfile& log = GetLogfile();
//	log.printf("seg_up iteration # %d\n", niter+1);
	
	// project the surfaces onto each other
	// this will update the gap functions as well
	ProjectSurface(m_ss, m_ms, bupseg);
	if (m_btwo_pass) ProjectSurface(m_ms, m_ss, bupseg);

	// Update the net contact pressures
	UpdateContactPressures();

	// set poro flag
	bool bporo = (m_ss.m_bporo || m_ms.m_bporo);

	// only continue if we are doing a poro-elastic simulation
	if (bporo == false) return;

	// update node normals
	m_ss.UpdateNodeNormals();
	if (bporo) m_ms.UpdateNodeNormals();

	// Now that the nodes have been projected, we need to figure out
	// if we need to modify the constraints on the pressure dofs.
	// If the nodes are not in contact, they must be free
	// draining. Since all nodes have been previously marked to be
	// free-draining in MarkFreeDraining(), we just need to reverse
	// this setting here, for nodes that are in contact.

	// Next, we loop over each surface, visiting the nodes
	// and finding out if that node is in contact or not
	int npass = (m_btwo_pass?2:1);
	for (np=0; np<npass; ++np)
	{
		FESlidingSurface2& ss = (np == 0? m_ss : m_ms);
		FESlidingSurface2& ms = (np == 0? m_ms : m_ss);

		// keep a running counter of integration points.
		int ni = 0;

		// loop over all elements of the primary surface
		for (n=0; n<ss.Elements(); ++n)
		{
			FESurfaceElement& el = ss.Element(n);
			int nint = el.GaussPoints();
			int neln = el.Nodes();

			// get the normal tractions at the integration points
			double ti[4], gap, eps;
			for (i=0; i<nint; ++i, ++ni) 
			{
				gap = ss.m_gap[ni];
				eps = m_epsn*ss.m_epsn[ni];
				ti[i] = MBRACKET(ss.m_Lmd[ni] + eps*gap);
			}

			// project the data to the nodes
			double tn[4];
			el.project_to_nodes(ti, tn);

			// see if we need to make changes to the BC's
			for (i=0; i<neln; ++i)
			{
				FENode& node = ss.Node(el.m_lnode[i]);
				id = node.m_ID[DOF_P];
				if ((id < -1) && (tn[i] > 0))
				{
					// mark node as non-free-draining (= pos ID)
					node.m_ID[DOF_P] = -id-2;
				}
			}
		}

		// loop over all nodes of the secondary surface
		// the secondary surface is trickier since we need
		// to look at the primary's surface projection
		for (n=0; n<ms.Nodes(); ++n)
		{
			// get the node
			FENode& node = ms.Node(n);
			
			// project it onto the primary surface
			int nei;
			FESurfaceElement* pse = ss.FindIntersection(node.m_rt, ms.m_nn[n], rs, m_stol, &nei);

			if (pse)
			{
				// we found an element, so let's see if it's even remotely close to contact
				// find the global location of the intersection point
				vec3d q = ss.Local2Global(*pse, rs[0], rs[1]);

				// calculate the gap function
				double g = ms.m_nn[n]*(node.m_rt - q);

				if (fabs(g) <= R)
				{
					// we found an element so let's calculate the nodal traction values for this element
					// get the normal tractions at the integration points
					double ti[4], gap, eps;
					int nint = pse->GaussPoints();
					int noff = ss.m_nei[nei];
					for (i=0; i<nint; ++i) 
					{
						gap = ss.m_gap[noff + i];
						eps = m_epsn*ss.m_epsn[noff + i];
						ti[i] = MBRACKET(ss.m_Lmd[noff + i] + eps*gap);
					}

					// project the data to the nodes
					double tn[4];
					pse->project_to_nodes(ti, tn);

					// now evaluate the traction at the intersection point
					double tp = pse->eval(tn, rs[0], rs[1]);

					// if tp > 0, mark node as non-free-draining. (= pos ID)
					id = node.m_ID[DOF_P];
					if ((id < -1) && (tp > 0))
					{
						// mark as non free-draining
						node.m_ID[DOF_P] = -id-2;
					}
				}
			}
		}
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
	double N[32];

	FEM& fem = dynamic_cast<FEM&>(*m_pfem);

	// get the mesh
	FEMesh* pm = m_ss.GetMesh();

	// get the solver
	FESolidSolver* psolver = dynamic_cast<FESolidSolver*>(fem.m_pStep->m_psolver);

	// if we're using the symmetric formulation
	// we need to multiply with the timestep
	double dt = fem.m_pStep->m_dt;

	// loop over the nr of passes
	int npass = (m_btwo_pass?2:1);
	for (int np=0; np<npass; ++np)
	{
		// get slave and master surface
		FESlidingSurface2& ss = (np == 0? m_ss : m_ms);
		FESlidingSurface2& ms = (np == 0? m_ms : m_ss);

		// keep a running counter of integration points
		int ni = 0;

		// loop over all slave elements
		for (i=0; i<ss.Elements(); ++i)
		{
			// get the surface element
			FESurfaceElement& se = ss.Element(i);

			bool sporo = PoroStatus(*pm, se);

			// get the nr of nodes and integration points
			int nseln = se.Nodes();
			int nint = se.GaussPoints();

			// get the element's LM vector
			ss.UnpackLM(se, sLM);

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

					bool mporo = PoroStatus(*pm, me);

					// get the nr of master element nodes
					int nmeln = me.Nodes();

					// get the element's LM vector
					ms.UnpackLM(me, mLM);

					// calculate degrees of freedom
					int ndof = 3*(nseln + nmeln);

					// build the LM vector
					LM.resize(ndof);
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
					en.resize(nseln+nmeln);
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
					double Lm = ss.m_Lmd[ni];

					// penalty 
					double eps = m_epsn*ss.m_epsn[ni];

					// contact traction
					double tn = Lm + eps*g;
					tn = MBRACKET(tn);

					// calculate the force vector
					fe.resize(ndof);
					zero(fe);

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

					for (k=0; k<ndof; ++k) fe[k] += tn*N[k]*detJ[j]*w[j];

					// assemble the global residual
					psolver->AssembleResidual(en, LM, fe, F);

					// do the biphasic stuff
					// TODO: I should only do this when the node is actually in contact
					//       in other words, when g >= 0
					if (sporo && mporo && (tn > 0))
					{
						// calculate nr of pressure dofs
						int ndof = nseln + nmeln;

						// calculate the flow rate
						double epsp = m_epsp*ss.m_epsp[ni];

						double wn = ss.m_Lmp[ni] + epsp*ss.m_pg[ni];

						// fill the LM
						LM.resize(ndof);
						for (k=0; k<nseln; ++k) LM[k        ] = sLM[3*nseln+k];
						for (k=0; k<nmeln; ++k) LM[k + nseln] = mLM[3*nmeln+k];

						// fill the force array
						fe.resize(ndof);
						zero(fe);
						for (k=0; k<nseln; ++k) N[k      ] =  Hs[k];
						for (k=0; k<nmeln; ++k) N[k+nseln] = -Hm[k];

						for (k=0; k<ndof; ++k) fe[k] += dt*wn*N[k]*detJ[j]*w[j];

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
	double N[32];
	matrix ke;

	FEM& fem = dynamic_cast<FEM&>(*m_pfem);

	// get the mesh
	FEMesh* pm = m_ss.GetMesh();

	// get the solver
	FESolidSolver* psolver = dynamic_cast<FESolidSolver*>(fem.m_pStep->m_psolver);

	// see how many reformations we've had to do so far
	int nref = psolver->m_nref;

	// set higher order stiffness mutliplier
	// NOTE: this algrotihm doesn't really need this
	// but I've added this functionality to compare with the other contact 
	// algorithms and to see the effect of the different stiffness contributions
	double knmult = m_knmult;
	if (m_knmult < 0)
	{
		int ni = int(-m_knmult);
		if (nref >= ni)
		{
			knmult = 1; 
			clog.printf("Higher order stiffness terms included.\n");
		}
		else knmult = 0;
	}

	// do single- or two-pass
	int npass = (m_btwo_pass?2:1);
	for (int np=0; np < npass; ++np)
	{
		// get the slave and master surface
		FESlidingSurface2& ss = (np == 0? m_ss : m_ms);
		FESlidingSurface2& ms = (np == 0? m_ms : m_ss);

		// keep a running counter of the integration points
		int ni = 0;

		// loop over all slave elements
		for (i=0; i<ss.Elements(); ++i)
		{
			// get ths slave element
			FESurfaceElement& se = ss.Element(i);

			bool sporo = PoroStatus(*pm, se);

			// get nr of nodes and integration points
			int nseln = se.Nodes();
			int nint = se.GaussPoints();

			// nodal pressures
			double pn[4];
			for (j=0; j<nseln; ++j) pn[j] = ss.GetMesh()->Node(se.m_node[j]).m_pt;

			// get the element's LM vector
			ss.UnpackLM(se, sLM);

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
				if (sporo)
				{
					pt[j] = se.eval(pn, j);
					dpr[j] = se.eval_deriv1(pn, j);
					dps[j] = se.eval_deriv2(pn, j);
				}
			}

			// loop over all integration points
			for (j=0; j<nint; ++j, ++ni)
			{
				// get the master element
				FESurfaceElement* pme = ss.m_pme[ni];
				if (pme)
				{
					FESurfaceElement& me = *pme;

					bool mporo = PoroStatus(*pm, me);

					// get the nr of master nodes
					int nmeln = me.Nodes();

					// nodal pressure
					double pm[4];
					for (k=0; k<nmeln; ++k) pm[k] = ms.GetMesh()->Node(me.m_node[k]).m_pt;

					// get the element's LM vector
					ms.UnpackLM(me, mLM);
					
					int ndpn;	// number of dofs per node
					int ndof;	// number of dofs in stiffness matrix

					if (sporo && mporo) {
						// calculate degrees of freedom for biphasic-on-biphasic contact
						ndpn = 4;
						ndof = ndpn*(nseln+nmeln);
						
						// build the LM vector
						LM.resize(ndof);
						
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
					}
					
					else {
						// calculate degrees of freedom for biphasic-on-elastic or elastic-on-elastic contact
						ndpn = 3;
						ndof = ndpn*(nseln + nmeln);
						
						// build the LM vector
						LM.resize(ndof);
						
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
					}
					
					// build the en vector
					en.resize(nseln+nmeln);
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
					double Lm = ss.m_Lmd[ni];

					// penalty 
					double eps = m_epsn*ss.m_epsn[ni];

					// contact traction
					double tn = Lm + eps*g;
					tn = MBRACKET(tn);

//					double dtn = m_eps*HEAVYSIDE(Lm + eps*g);
					double dtn = (tn > 0.? eps :0.);
					
					// create the stiffness matrix
					ke.Create(ndof, ndof); ke.zero();
					
					// --- S O L I D - S O L I D   C O N T A C T ---
					
					// a. NxN-term
					//------------------------------------
					
					// calculate the N-vector
					for (k=0; k<nseln; ++k)
					{
						N[ndpn*k  ] = Hs[k]*nu.x;
						N[ndpn*k+1] = Hs[k]*nu.y;
						N[ndpn*k+2] = Hs[k]*nu.z;
					}
					
					for (k=0; k<nmeln; ++k)
					{
						N[ndpn*(k+nseln)  ] = -Hm[k]*nu.x;
						N[ndpn*(k+nseln)+1] = -Hm[k]*nu.y;
						N[ndpn*(k+nseln)+2] = -Hm[k]*nu.z;
					}
					
					if (ndpn == 4) {
						for (k=0; k<nseln; ++k)
							N[ndpn*k+3] = 0;
						for (k=0; k<nmeln; ++k)
							N[ndpn*(k+nseln)+3] = 0;
					}
					
					for (k=0; k<ndof; ++k)
						for (l=0; l<ndof; ++l) ke[k][l] += dtn*N[k]*N[l]*detJ[j]*w[j];
					
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
					mat3d As[4];
					for (l=0; l<nseln; ++l)
						As[l] = S2*Gr[l] - S1*Gs[l];
					
					if (!m_bsymm)
					{	// non-symmetric
						for (l=0; l<nseln; ++l)
						{
							for (k=0; k<nseln+nmeln; ++k)
							{
								ke[k*ndpn  ][l*ndpn  ] -= knmult*tn*w[j]*N[k]*As[l][0][0];
								ke[k*ndpn  ][l*ndpn+1] -= knmult*tn*w[j]*N[k]*As[l][0][1];
								ke[k*ndpn  ][l*ndpn+2] -= knmult*tn*w[j]*N[k]*As[l][0][2];
								
								ke[k*ndpn+1][l*ndpn  ] -= knmult*tn*w[j]*N[k]*As[l][1][0];
								ke[k*ndpn+1][l*ndpn+1] -= knmult*tn*w[j]*N[k]*As[l][1][1];
								ke[k*ndpn+1][l*ndpn+2] -= knmult*tn*w[j]*N[k]*As[l][1][2];
								
								ke[k*ndpn+2][l*ndpn  ] -= knmult*tn*w[j]*N[k]*As[l][2][0];
								ke[k*ndpn+2][l*ndpn+1] -= knmult*tn*w[j]*N[k]*As[l][2][1];
								ke[k*ndpn+2][l*ndpn+2] -= knmult*tn*w[j]*N[k]*As[l][2][2];
							}
						}
					} 
					else 
					{	// symmetric
						for (l=0; l<nseln; ++l)
						{
							for (k=0; k<nseln+nmeln; ++k)
							{
								ke[k*ndpn  ][l*ndpn  ] -= 0.5*knmult*tn*w[j]*N[k]*As[l][0][0];
								ke[k*ndpn  ][l*ndpn+1] -= 0.5*knmult*tn*w[j]*N[k]*As[l][0][1];
								ke[k*ndpn  ][l*ndpn+2] -= 0.5*knmult*tn*w[j]*N[k]*As[l][0][2];
								
								ke[k*ndpn+1][l*ndpn  ] -= 0.5*knmult*tn*w[j]*N[k]*As[l][1][0];
								ke[k*ndpn+1][l*ndpn+1] -= 0.5*knmult*tn*w[j]*N[k]*As[l][1][1];
								ke[k*ndpn+1][l*ndpn+2] -= 0.5*knmult*tn*w[j]*N[k]*As[l][1][2];
								
								ke[k*ndpn+2][l*ndpn  ] -= 0.5*knmult*tn*w[j]*N[k]*As[l][2][0];
								ke[k*ndpn+2][l*ndpn+1] -= 0.5*knmult*tn*w[j]*N[k]*As[l][2][1];
								ke[k*ndpn+2][l*ndpn+2] -= 0.5*knmult*tn*w[j]*N[k]*As[l][2][2];
								
								ke[l*ndpn  ][k*ndpn  ] -= 0.5*knmult*tn*w[j]*N[k]*As[l][0][0];
								ke[l*ndpn+1][k*ndpn  ] -= 0.5*knmult*tn*w[j]*N[k]*As[l][0][1];
								ke[l*ndpn+2][k*ndpn  ] -= 0.5*knmult*tn*w[j]*N[k]*As[l][0][2];
								
								ke[l*ndpn  ][k*ndpn+1] -= 0.5*knmult*tn*w[j]*N[k]*As[l][1][0];
								ke[l*ndpn+1][k*ndpn+1] -= 0.5*knmult*tn*w[j]*N[k]*As[l][1][1];
								ke[l*ndpn+2][k*ndpn+1] -= 0.5*knmult*tn*w[j]*N[k]*As[l][1][2];
								
								ke[l*ndpn  ][k*ndpn+2] -= 0.5*knmult*tn*w[j]*N[k]*As[l][2][0];
								ke[l*ndpn+1][k*ndpn+2] -= 0.5*knmult*tn*w[j]*N[k]*As[l][2][1];
								ke[l*ndpn+2][k*ndpn+2] -= 0.5*knmult*tn*w[j]*N[k]*As[l][2][2];
							}
						}
					}
					
					// c. M-term
					//---------------------------------------
					
					vec3d Gm[2];
					ms.ContraBaseVectors(me, r, s, Gm);
					
					// evaluate master surface normal
					vec3d mnu = Gm[0] ^ Gm[1];
					mnu.unit();
					
					double Hmr[4], Hms[4];
					me.shape_deriv(Hmr, Hms, r, s);
					vec3d mm[4];
					for (k=0; k<nmeln; ++k) 
						mm[k] = Gm[0]*Hmr[k] + Gm[1]*Hms[k];
					
					if (!m_bsymm)
					{	// non-symmetric
						for (k=0; k<nmeln; ++k) 
						{
							for (l=0; l<nseln+nmeln; ++l)
							{
								ke[(k+nseln)*ndpn  ][l*ndpn  ] += tn*knmult*detJ[j]*w[j]*mnu.x*mm[k].x*N[l];
								ke[(k+nseln)*ndpn  ][l*ndpn+1] += tn*knmult*detJ[j]*w[j]*mnu.x*mm[k].y*N[l];
								ke[(k+nseln)*ndpn  ][l*ndpn+2] += tn*knmult*detJ[j]*w[j]*mnu.x*mm[k].z*N[l];
								
								ke[(k+nseln)*ndpn+1][l*ndpn  ] += tn*knmult*detJ[j]*w[j]*mnu.y*mm[k].x*N[l];
								ke[(k+nseln)*ndpn+1][l*ndpn+1] += tn*knmult*detJ[j]*w[j]*mnu.y*mm[k].y*N[l];
								ke[(k+nseln)*ndpn+1][l*ndpn+2] += tn*knmult*detJ[j]*w[j]*mnu.y*mm[k].z*N[l];
								
								ke[(k+nseln)*ndpn+2][l*ndpn  ] += tn*knmult*detJ[j]*w[j]*mnu.z*mm[k].x*N[l];
								ke[(k+nseln)*ndpn+2][l*ndpn+1] += tn*knmult*detJ[j]*w[j]*mnu.z*mm[k].y*N[l];
								ke[(k+nseln)*ndpn+2][l*ndpn+2] += tn*knmult*detJ[j]*w[j]*mnu.z*mm[k].z*N[l];
							}
						}
					}
					else
					{	// symmetric
						for (k=0; k<nmeln; ++k) 
						{
							for (l=0; l<nseln+nmeln; ++l)
							{
								ke[(k+nseln)*ndpn  ][l*ndpn  ] += 0.5*knmult*tn*detJ[j]*w[j]*mnu.x*mm[k].x*N[l];
								ke[(k+nseln)*ndpn  ][l*ndpn+1] += 0.5*knmult*tn*detJ[j]*w[j]*mnu.x*mm[k].y*N[l];
								ke[(k+nseln)*ndpn  ][l*ndpn+2] += 0.5*knmult*tn*detJ[j]*w[j]*mnu.x*mm[k].z*N[l];
								
								ke[(k+nseln)*ndpn+1][l*ndpn  ] += 0.5*knmult*tn*detJ[j]*w[j]*mnu.y*mm[k].x*N[l];
								ke[(k+nseln)*ndpn+1][l*ndpn+1] += 0.5*knmult*tn*detJ[j]*w[j]*mnu.y*mm[k].y*N[l];
								ke[(k+nseln)*ndpn+1][l*ndpn+2] += 0.5*knmult*tn*detJ[j]*w[j]*mnu.y*mm[k].z*N[l];
								
								ke[(k+nseln)*ndpn+2][l*ndpn  ] += 0.5*knmult*tn*detJ[j]*w[j]*mnu.z*mm[k].x*N[l];
								ke[(k+nseln)*ndpn+2][l*ndpn+1] += 0.5*knmult*tn*detJ[j]*w[j]*mnu.z*mm[k].y*N[l];
								ke[(k+nseln)*ndpn+2][l*ndpn+2] += 0.5*knmult*tn*detJ[j]*w[j]*mnu.z*mm[k].z*N[l];
								
								ke[l*ndpn  ][(k+nseln)*ndpn  ] += 0.5*knmult*tn*detJ[j]*w[j]*mnu.x*mm[k].x*N[l];
								ke[l*ndpn+1][(k+nseln)*ndpn  ] += 0.5*knmult*tn*detJ[j]*w[j]*mnu.y*mm[k].x*N[l];
								ke[l*ndpn+2][(k+nseln)*ndpn  ] += 0.5*knmult*tn*detJ[j]*w[j]*mnu.z*mm[k].x*N[l];
								
								ke[l*ndpn  ][(k+nseln)*ndpn+1] += 0.5*knmult*tn*detJ[j]*w[j]*mnu.x*mm[k].y*N[l];
								ke[l*ndpn+1][(k+nseln)*ndpn+1] += 0.5*knmult*tn*detJ[j]*w[j]*mnu.y*mm[k].y*N[l];
								ke[l*ndpn+2][(k+nseln)*ndpn+1] += 0.5*knmult*tn*detJ[j]*w[j]*mnu.z*mm[k].y*N[l];
								
								ke[l*ndpn  ][(k+nseln)*ndpn+2] += 0.5*knmult*tn*detJ[j]*w[j]*mnu.x*mm[k].z*N[l];
								ke[l*ndpn+1][(k+nseln)*ndpn+2] += 0.5*knmult*tn*detJ[j]*w[j]*mnu.y*mm[k].z*N[l];
								ke[l*ndpn+2][(k+nseln)*ndpn+2] += 0.5*knmult*tn*detJ[j]*w[j]*mnu.z*mm[k].z*N[l];
							}
						}
					}
					
					// --- B I P H A S I C   S T I F F N E S S ---
					if (sporo && mporo)
					{
						// the variable dt is either the timestep or one
						// depending on whether we are using the symmetric
						// poro version or not.
						double dt = fem.m_pStep->m_dt;
						
						double epsp = (tn > 0) ? m_epsp*ss.m_epsp[ni] : 0.;
						
						// --- S O L I D - P R E S S U R E   C O N T A C T ---
						
						if (!m_bsymm)
						{
							
							// a. q-term
							//-------------------------------------
							
							double dpmr, dpms;
							dpmr = me.eval_deriv1(pm, r, s);
							dpms = me.eval_deriv2(pm, r, s);
							
							for (k=0; k<nseln+nmeln; ++k)
								for (l=0; l<nseln+nmeln; ++l)
								{
									ke[4*k + 3][4*l  ] += dt*w[j]*detJ[j]*epsp*N[k]*N[l]*(dpmr*Gm[0].x + dpms*Gm[1].x);
									ke[4*k + 3][4*l+1] += dt*w[j]*detJ[j]*epsp*N[k]*N[l]*(dpmr*Gm[0].y + dpms*Gm[1].y);
									ke[4*k + 3][4*l+2] += dt*w[j]*detJ[j]*epsp*N[k]*N[l]*(dpmr*Gm[0].z + dpms*Gm[1].z);
								}
							
							double wn = ss.m_Lmp[ni] + epsp*ss.m_pg[ni];
							
							// b. A-term
							//-------------------------------------
							
							for (l=0; l<nseln; ++l)
								for (k=0; k<nseln+nmeln; ++k)
								{
									ke[4*k + 3][4*l  ] -= dt*w[j]*wn*N[k]*(As[l][0][0]*nu.x + As[l][0][1]*nu.y + As[l][0][2]*nu.z);
									ke[4*k + 3][4*l+1] -= dt*w[j]*wn*N[k]*(As[l][1][0]*nu.x + As[l][1][1]*nu.y + As[l][1][2]*nu.z);
									ke[4*k + 3][4*l+2] -= dt*w[j]*wn*N[k]*(As[l][2][0]*nu.x + As[l][2][1]*nu.y + As[l][2][2]*nu.z);
								}
							
							// c. m-term
							//---------------------------------------
							
							for (k=0; k<nmeln; ++k)
								for (l=0; l<nseln+nmeln; ++l)
								{
									ke[4*(k+nseln) + 3][4*l  ] += dt*w[j]*detJ[j]*wn*N[l]*mm[k].x;
									ke[4*(k+nseln) + 3][4*l+1] += dt*w[j]*detJ[j]*wn*N[l]*mm[k].y;
									ke[4*(k+nseln) + 3][4*l+2] += dt*w[j]*detJ[j]*wn*N[l]*mm[k].z;
								}
						}
						
						
						// --- P R E S S U R E - P R E S S U R E   C O N T A C T ---
						
						// calculate the N-vector
						for (k=0; k<nseln; ++k)
						{
							N[ndpn*k  ] = 0;
							N[ndpn*k+1] = 0;
							N[ndpn*k+2] = 0;
							N[ndpn*k+3] = Hs[k];
						}
						
						for (k=0; k<nmeln; ++k)
						{
							N[ndpn*(k+nseln)  ] = 0;
							N[ndpn*(k+nseln)+1] = 0;
							N[ndpn*(k+nseln)+2] = 0;
							N[ndpn*(k+nseln)+3] = -Hm[k];
						}
						
						for (k=0; k<ndof; ++k)
							for (l=0; l<ndof; ++l) ke[k][l] -= dt*epsp*w[j]*detJ[j]*N[k]*N[l];
						
					}
					
					// assemble the global stiffness
					psolver->AssembleStiffness(en, LM, ke);
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------
void FESlidingInterface2::UpdateContactPressures()
{
	int np, n, i, j, k;
	int npass = (m_btwo_pass?2:1);
	for (np=0; np<npass; ++np)
	{
		FESlidingSurface2& ss = (np == 0? m_ss : m_ms);
		FESlidingSurface2& ms = (np == 0? m_ms : m_ss);
		
		// keep a running counter of integration points.
		int ni = 0;
		
		// loop over all elements of the primary surface
		for (n=0; n<ss.Elements(); ++n)
		{
			FESurfaceElement& el = ss.Element(n);
			int nint = el.GaussPoints();
			
			// get the normal tractions at the integration points
			double gap, eps;
			for (i=0; i<nint; ++i, ++ni) 
			{
				gap = ss.m_gap[ni];
				eps = m_epsn*ss.m_epsn[ni];
				ss.m_Ln[ni] = MBRACKET(ss.m_Lmd[ni] + eps*gap);
				FESurfaceElement* pme = ss.m_pme[ni];
				if (m_btwo_pass && pme)
				{
					int mint = pme->GaussPoints();
					int noff = ms.m_nei[pme->m_lid];
					double ti[4];
					for (j=0; j<mint; ++j) {
						k = noff+j;
						gap = ms.m_gap[k];
						eps = m_epsn*ms.m_epsn[k];
						ti[j] = MBRACKET(ms.m_Lmd[k] + m_epsn*ms.m_epsn[k]*ms.m_gap[k]);
					}
					// project the data to the nodes
					double tn[4];
					pme->project_to_nodes(ti, tn);
					// now evaluate the traction at the intersection point
					double Ln = pme->eval(tn, ss.m_rs[ni][0], ss.m_rs[ni][1]);
					ss.m_Ln[ni] += MBRACKET(Ln);
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------
bool FESlidingInterface2::Augment(int naug)
{
	// make sure we need to augment
	if (!m_blaugon) return true;

	int i;
	double Ln, Lp;
	bool bconv = true;

	bool bporo = (m_ss.m_bporo && m_ms.m_bporo);
	int NS = m_ss.m_Lmd.size();
	int NM = m_ms.m_Lmd.size();

	// --- c a l c u l a t e   i n i t i a l   n o r m s ---
	// a. normal component
	double normL0 = 0, normP = 0, normDP = 0;
	for (i=0; i<NS; ++i) normL0 += m_ss.m_Lmd[i]*m_ss.m_Lmd[i];
	for (i=0; i<NM; ++i) normL0 += m_ms.m_Lmd[i]*m_ms.m_Lmd[i];

	// b. gap component
	// (is calculated during update)
	double maxgap = 0;
	double maxpg = 0;

	// update Lagrange multipliers
	double normL1 = 0, eps, epsp;
	for (i=0; i<NS; ++i)
	{
		// update Lagrange multipliers on slave surface
		eps = m_epsn*m_ss.m_epsn[i];
		Ln = m_ss.m_Lmd[i] + eps*m_ss.m_gap[i];
		m_ss.m_Lmd[i] = MBRACKET(Ln);

		normL1 += m_ss.m_Lmd[i]*m_ss.m_Lmd[i];
		
		if (m_ss.m_bporo) {
			Lp = 0;
			if (Ln > 0) {
				epsp = m_epsp*m_ss.m_epsp[i];
				Lp = m_ss.m_Lmp[i] + epsp*m_ss.m_pg[i];
				maxpg = max(maxpg,fabs(m_ss.m_pg[i]));
				normDP += m_ss.m_pg[i]*m_ss.m_pg[i];
			}
			m_ss.m_Lmp[i] = Lp;
		}
		
		if (Ln > 0) maxgap = max(maxgap,fabs(m_ss.m_gap[i]));
	}	
	
	for (i=0; i<NM; ++i)
	{
		// update Lagrange multipliers on master surface
		eps = m_epsn*m_ms.m_epsn[i];
		Ln = m_ms.m_Lmd[i] + eps*m_ms.m_gap[i];
		m_ms.m_Lmd[i] = MBRACKET(Ln);

		normL1 += m_ms.m_Lmd[i]*m_ms.m_Lmd[i];
		
		if (m_ms.m_bporo) {
			Lp = 0;
			if (Ln > 0) {
				epsp = m_epsp*m_ms.m_epsp[i];
				Lp = m_ms.m_Lmp[i] + epsp*m_ms.m_pg[i];
				maxpg = max(maxpg,fabs(m_ms.m_pg[i]));
				normDP += m_ms.m_pg[i]*m_ms.m_pg[i];
			}
			m_ms.m_Lmp[i] = Lp;
		}
		
		if (Ln > 0) maxgap = max(maxgap,fabs(m_ms.m_gap[i]));
	}

	// Ideally normP should be evaluated from the fluid pressure at the
	// contact interface (not easily accessible).  The next best thing
	// is to use the contact traction.
	normP = normL1;
	
	// calculate relative norms
	double lnorm = (normL1 != 0 ? fabs((normL1 - normL0) / normL1) : fabs(normL1 - normL0)); 
	double pnorm = (normP != 0 ? (normDP/normP) : normDP); 

	// check convergence
	if ((m_gtol > 0) && (maxgap > m_gtol)) bconv = false;
	if ((m_ptol > 0) && (bporo && maxpg > m_ptol)) bconv = false;

	if ((m_atol > 0) && (lnorm > m_atol)) bconv = false;
	if ((m_atol > 0) && (pnorm > m_atol)) bconv = false;

	clog.printf(" sliding interface # %d\n", m_nID);
	clog.printf("                        CURRENT        REQUIRED\n");
	clog.printf("    D multiplier : %15le", lnorm); if (m_atol > 0) clog.printf("%15le\n", m_atol); else clog.printf("       ***\n");
	if (bporo) { clog.printf("    P gap       : %15le", pnorm); if (m_atol > 0) clog.printf("%15le\n", m_atol); else clog.printf("       ***\n"); }

	clog.printf("    maximum gap  : %15le", maxgap);
	if (m_gtol > 0) clog.printf("%15le\n", m_gtol); else clog.printf("       ***\n");
	if (bporo) {
		clog.printf("    maximum pgap : %15le", maxpg);
		if (m_ptol > 0) clog.printf("%15le\n", m_ptol); else clog.printf("       ***\n");
	}
	
	return bconv;
}

//-----------------------------------------------------------------------------
void FESlidingInterface2::Serialize(DumpFile &ar)
{
	FEContactInterface::Serialize(ar);
	if (ar.IsSaving())
	{
		ar << m_knmult;
		ar << m_btwo_pass;
		ar << m_atol;
		ar << m_gtol;
		ar << m_ptol;
		ar << m_stol;
		ar << m_bsymm;
		ar << m_srad;
		ar << m_naugmax;
		ar << m_naugmin;
		ar << m_nsegup;
		ar << m_epsn;
		ar << m_bautopen;
		ar << m_epsp;

		m_ms.Serialize(ar);
		m_ss.Serialize(ar);
	}
	else
	{
		ar >> m_knmult;
		ar >> m_btwo_pass;
		ar >> m_atol;
		ar >> m_gtol;
		ar >> m_ptol;
		ar >> m_stol;
		ar >> m_bsymm;
		ar >> m_srad;
		ar >> m_naugmax;
		ar >> m_naugmin;
		ar >> m_nsegup;
		ar >> m_epsn;
		ar >> m_bautopen;
		ar >> m_epsp;

		m_ms.Serialize(ar);
		m_ss.Serialize(ar);
	}
}

//-----------------------------------------------------------------------------

bool FESlidingInterface2::PoroStatus(FEMesh& m, FESurfaceElement& el)
{
	bool status = false;

	// get the solid element this surface element belongs to
	FESolidElement* pe = dynamic_cast<FESolidElement*>(m.FindElementFromID(el.m_nelem));
	if (pe)
	{
		// get the material
		FEMaterial* pm = dynamic_cast<FEMaterial*>(m_pfem->GetMaterial(pe->GetMatID()));
		
		// see if this is a poro-elastic element
		FEPoroElastic* poro = dynamic_cast<FEPoroElastic*>(pm);
		FEBiphasic* biph = dynamic_cast<FEBiphasic*> (pm);
		if (poro || biph) status = true;
	}
	
	return status;
}

//-----------------------------------------------------------------------------

void FESlidingInterface2::MarkFreeDraining()
{	
	int i, id, np;

	// Mark all nodes as free-draining.  This needs to be done for ALL
	// contact interfaces prior to executing Update(), where nodes that are
	// in contact are subsequently marked as non free-draining.  This ensures
	// that for surfaces involved in more than one contact interface, nodes
	// that have been marked as non free-draining are not reset to 
	// free-draining.
	for (np=0; np<2; ++np)
	{
		FESlidingSurface2& s = (np == 0? m_ss : m_ms);
		
		if (s.m_bporo) {
			// first, mark all nodes as free-draining (= neg. ID)
			// this is done by setting the dof's equation number
			// to a negative number
			for (i=0; i<s.Nodes(); ++i) 
			{
				id = s.Node(i).m_ID[DOF_P];
				if (id >= 0) 
				{
					FENode& node = s.Node(i);
					// mark node as free-draining
					node.m_ID[DOF_P] = -id-2;
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------

void FESlidingInterface2::SetFreeDraining()
{	
	int i, np;
	
	// Set the pressure to zero for the free-draining nodes
	for (np=0; np<2; ++np)
	{
		FESlidingSurface2& s = (np == 0? m_ss : m_ms);
		
		if (s.m_bporo) {
			// loop over all nodes
			for (i=0; i<s.Nodes(); ++i) 
			{
				if (s.Node(i).m_ID[DOF_P] < -1)
				{
					FENode& node = s.Node(i);
					// set the fluid pressure to zero
					node.m_pt = 0;
				}
			}
		}
	}
}
