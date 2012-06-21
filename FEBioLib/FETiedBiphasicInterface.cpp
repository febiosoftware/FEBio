#include "stdafx.h"
#include "FETiedBiphasicInterface.h"
#include "FESolidSolver.h"
#include "FEBiphasic.h"
#include "log.h"

//-----------------------------------------------------------------------------
// Define sliding interface parameters
BEGIN_PARAMETER_LIST(FETiedBiphasicInterface, FEContactInterface)
ADD_PARAMETER(m_blaugon  , FE_PARAM_BOOL  , "laugon"             );
ADD_PARAMETER(m_atol     , FE_PARAM_DOUBLE, "tolerance"          );
ADD_PARAMETER(m_gtol     , FE_PARAM_DOUBLE, "gaptol"             );
ADD_PARAMETER(m_ptol     , FE_PARAM_DOUBLE, "ptol"               );
ADD_PARAMETER(m_epsn     , FE_PARAM_DOUBLE, "penalty"            );
ADD_PARAMETER(m_bautopen , FE_PARAM_BOOL  , "auto_penalty"       );
ADD_PARAMETER(m_btwo_pass, FE_PARAM_BOOL  , "two_pass"           );
ADD_PARAMETER(m_knmult   , FE_PARAM_DOUBLE, "knmult"             );
ADD_PARAMETER(m_stol     , FE_PARAM_DOUBLE, "search_tol"         );
ADD_PARAMETER(m_epsp     , FE_PARAM_DOUBLE, "pressure_penalty"   );
ADD_PARAMETER(m_bsymm    , FE_PARAM_BOOL  , "symmetric_stiffness");
ADD_PARAMETER(m_srad     , FE_PARAM_DOUBLE, "search_radius"      );
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
// FETiedBiphasicSurface
//-----------------------------------------------------------------------------

FETiedBiphasicSurface::FETiedBiphasicSurface(FEModel* pfem) : FEContactSurface(&pfem->m_mesh)
{ 
	m_bporo = false;
	m_pfem = pfem; 
}

//-----------------------------------------------------------------------------
bool FETiedBiphasicSurface::Init()
{
	// initialize surface data first
	if (FEContactSurface::Init() == false) return false;
	
	// count how many integration points we have
	int nint = 0, i;
	for (i=0; i<Elements(); ++i)
	{
		FESurfaceElement& el = Element(i);
		nint += el.GaussPoints();
	}
	
	// allocate data structures
	m_Gap.assign(nint, vec3d(0,0,0));
	m_dg.assign(nint, vec3d(0,0,0));
	m_nu.resize(nint);
	m_rs.resize(nint);
	m_Lmd.assign(nint, vec3d(0,0,0));
	m_pme.assign(nint, static_cast<FESurfaceElement*>(0));
	m_epsn.assign(nint, 1.0);
	
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
			FEBiphasic* biph = dynamic_cast<FEBiphasic*> (pm);
			if (biph) m_bporo = true;
		}
		++i;
	}
	
	// allocate biphasic stuff
	if (m_bporo)
	{
		m_pg.assign(nint, 0);
		m_Lmp.assign(nint, 0.0);
		m_epsp.assign(nint, 1.0);
	}
	
	return true;
}

//-----------------------------------------------------------------------------
void FETiedBiphasicSurface::ShallowCopy(FETiedBiphasicSurface &s)
{
	m_Lmd = s.m_Lmd;
	m_Gap = s.m_Gap;
	m_dg = s.m_dg;
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

void FETiedBiphasicSurface::UpdateNodeNormals()
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
void FETiedBiphasicSurface::Serialize(DumpFile& ar)
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
		ar << m_Gap;
		ar << m_dg;
		ar << m_nu;
		ar << m_rs;
		ar << m_Lmd;
		ar << m_Lmp;
		ar << m_nei;
		ar << m_epsn;
		ar << m_epsp;
		ar << m_nn;
		ar << m_pg;
	}
	else
	{
		ar >> m_Gap;
		ar >> m_dg;
		ar >> m_nu;
		ar >> m_rs;
		ar >> m_Lmd;
		ar >> m_Lmp;
		ar >> m_nei;
		ar >> m_epsn;
		ar >> m_epsp;
		ar >> m_nn;
		ar >> m_pg;
	}
}

//-----------------------------------------------------------------------------
// FETiedBiphasicInterface
//-----------------------------------------------------------------------------

FETiedBiphasicInterface::FETiedBiphasicInterface(FEModel* pfem) : FEContactInterface(pfem), m_ss(pfem), m_ms(pfem)
{
	m_ntype = FE_CONTACT_TIED_BIPHASIC;
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
	m_srad = 1.0;
	m_gtol = -1;	// we use augmentation tolerance by default
	m_ptol = -1;	// we use augmentation tolerance by default
	m_bautopen = false;
	
	m_naugmin = 0;
	m_naugmax = 10;
	
	m_ss.SetSibling(&m_ms);
	m_ms.SetSibling(&m_ss);
}

//-----------------------------------------------------------------------------

FETiedBiphasicInterface::~FETiedBiphasicInterface()
{
}

//-----------------------------------------------------------------------------
bool FETiedBiphasicInterface::Init()
{
	// initialize surface data
	if (m_ss.Init() == false) return false;
	if (m_ms.Init() == false) return false;
	
	bool bporo = (m_ss.m_bporo || m_ms.m_bporo);
	
	FENLSolver* psolver = m_pfem->GetCurrentStep()->m_psolver;
	
	// this contact implementation requires a non-symmetric stiffness matrix
	// so inform the FEM class
	if (!m_bsymm) 
	{
		// request a non-symmetric stiffness matrix
		psolver->m_bsymm = false;
	}
	
	// calculate the penalty
	if (m_bautopen) 
	{
		CalcAutoPenalty(m_ss);
		if (m_ss.m_bporo) CalcAutoPressurePenalty(m_ss);
		if (m_btwo_pass) {
			CalcAutoPenalty(m_ms);
			if (m_ms.m_bporo) CalcAutoPressurePenalty(m_ms);
		}
	}
	
	// project the surfaces onto each other
	// this will evaluate the gap functions in the reference configuration
	InitialProjection(m_ss, m_ms);
	if (m_btwo_pass) InitialProjection(m_ms, m_ss);
	
	return true;
}

//-----------------------------------------------------------------------------
void FETiedBiphasicInterface::CalcAutoPenalty(FETiedBiphasicSurface& s)
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
		double E = BulkModulus(el, s);
		
		// calculate penalty
		double eps = E*A/V;
		
		// assign to integation points of surface element
		int nint = el.GaussPoints();
		for (int j=0; j<nint; ++j, ++ni) s.m_epsn[ni] = eps;
	}
}

//-----------------------------------------------------------------------------
void FETiedBiphasicInterface::CalcAutoPressurePenalty(FETiedBiphasicSurface& s)
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

double FETiedBiphasicInterface::AutoPressurePenalty(FESurfaceElement& el, FETiedBiphasicSurface& s)
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
		FEBiphasic* biph = dynamic_cast<FEBiphasic*> (pm);
		if (biph)
		{
			// get a material point
			FEMaterialPoint& mp = *pe->m_State[0];
			FEElasticMaterialPoint& ept = *(mp.ExtractData<FEElasticMaterialPoint>());
			
			// setup the material point
			ept.F = mat3dd(1.0);
			ept.J = 1;
			ept.s.zero();
			
			// if this is a poroelastic element, then get the permeability tensor
			FEBiphasicMaterialPoint& pt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
			pt.m_p = 0;
			pt.m_w = vec3d(0,0,0);
			
			double K[3][3];
			biph->Permeability(K, mp);
			
			eps = (K[0][0] + K[1][1] + K[2][2])/3;
		}
	}
	
	return eps;
}

//-----------------------------------------------------------------------------
// Perform initial projection between tied surfaces in reference configuration
void FETiedBiphasicInterface::InitialProjection(FETiedBiphasicSurface& ss, FETiedBiphasicSurface& ms)
{
	bool bfirst = true;
	
	FESurfaceElement* pme;
	vec3d r, nu;
	double rs[2];
	
	// loop over all integration points
	int n = 0;
	for (int i=0; i<ss.Elements(); ++i)
	{
		FESurfaceElement& el = ss.Element(i);
		
		int nint = el.GaussPoints();
		
		for (int j=0; j<nint; ++j, ++n)
		{
			// calculate the global position of the integration point
			r = ss.Local2Global(el, j);
			
			// calculate the normal at this integration point
			nu = ss.SurfaceNormal(el, j);
			
			// find the intersection point with the master surface
			pme = ms.FindIntersection(r, nu, rs, bfirst, m_stol, m_srad);
			
			ss.m_pme[n] = pme;
			ss.m_rs[n][0] = rs[0];
			ss.m_rs[n][1] = rs[1];
			if (pme)
			{
				// the node could potentially be in contact
				// find the global location of the intersection point
				vec3d q = ms.Local2Global(*pme, rs[0], rs[1]);
				
				// calculate the gap function
				ss.m_Gap[n] = q - r;
			}
			else
			{
				// the node is not in contact
				ss.m_Gap[n] = 0;
			}
		}
	}
}

//-----------------------------------------------------------------------------
// Evaluate gap functions for position and fluid pressure
void FETiedBiphasicInterface::ProjectSurface(FETiedBiphasicSurface& ss, FETiedBiphasicSurface& ms)
{
	FEMesh& mesh = m_pfem->m_mesh;
	FESurfaceElement* pme;
	vec3d r;
	
	double ps[4], p1;
	
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
			ss.m_nu[n] = ss.SurfaceNormal(el, j);

			// if this node is tied, evaluate gap functions
			pme = ss.m_pme[n];
			if (pme)
			{
				// find the global location of the intersection point
				vec3d q = ms.Local2Global(*pme, ss.m_rs[n][0], ss.m_rs[n][1]);
				
				// calculate the gap function
				vec3d g = q - r;
				ss.m_dg[n] = g - ss.m_Gap[n];
				
				// calculate the pressure gap function
				bool mporo = PoroStatus(mesh, *pme);
				if (sporo && mporo) {
					double pm[4];
					for (int k=0; k<pme->Nodes(); ++k) pm[k] = mesh.Node(pme->m_node[k]).m_pt;
					double p2 = pme->eval(pm, ss.m_rs[n][0], ss.m_rs[n][1]);
					ss.m_pg[n] = p1 - p2;
				}
			}
			else
			{
				// the node is not tied
				ss.m_dg[n] = vec3d(0,0,0);
				if (sporo) ss.m_pg[n] = 0;
			}
		}
	}
}

//-----------------------------------------------------------------------------

void FETiedBiphasicInterface::Update(int niter)
{	

	// project the surfaces onto each other
	// this will update the gap functions as well
	ProjectSurface(m_ss, m_ms);
	if (m_btwo_pass) ProjectSurface(m_ms, m_ss);
	
}

//-----------------------------------------------------------------------------
void FETiedBiphasicInterface::ShallowCopy(FEContactInterface &ci)
{
	FETiedBiphasicInterface& si = dynamic_cast<FETiedBiphasicInterface&>(ci);
	m_ss.ShallowCopy(si.m_ss);
	m_ms.ShallowCopy(si.m_ms);
}

//-----------------------------------------------------------------------------
void FETiedBiphasicInterface::ContactForces(vector<double> &F, FENLSolver* psolver)
{
	int i, j, k;
	vector<int> sLM, mLM, LM, en;
	vector<double> fe;
	double detJ[4], w[4], *Hs, Hm[4];
	double N[32];

	// get time step
	double dt = m_pfem->GetCurrentStep()->m_dt;
	
	// get the mesh
	FEMesh* pm = m_ss.GetMesh();
	
	// if we're using the symmetric formulation
	// we need to multiply with the timestep
//	double dt = fem.m_pStep->m_dt;
	
	// loop over the nr of passes
	int npass = (m_btwo_pass?2:1);
	for (int np=0; np<npass; ++np)
	{
		// get slave and master surface
		FETiedBiphasicSurface& ss = (np == 0? m_ss : m_ms);
		FETiedBiphasicSurface& ms = (np == 0? m_ms : m_ss);
		
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
			
			// copy the LM vector; we'll need it later
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
					
					// copy LM vector
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
					vec3d dg = ss.m_dg[ni];
					
					// lagrange multiplier
					vec3d Lm = ss.m_Lmd[ni];
					
					// penalty 
					double eps = m_epsn*ss.m_epsn[ni];
					
					// contact traction
					vec3d t = Lm + dg*eps;
					
					// calculate the force vector
					fe.resize(ndof);
					zero(fe);
					
					for (k=0; k<nseln; ++k)
					{
						N[3*k  ] = Hs[k]*t.x;
						N[3*k+1] = Hs[k]*t.y;
						N[3*k+2] = Hs[k]*t.z;
					}
					
					for (k=0; k<nmeln; ++k)
					{
						N[3*(k+nseln)  ] = -Hm[k]*t.x;
						N[3*(k+nseln)+1] = -Hm[k]*t.y;
						N[3*(k+nseln)+2] = -Hm[k]*t.z;
					}
					
					for (k=0; k<ndof; ++k) fe[k] += N[k]*detJ[j]*w[j];
					
					// assemble the global residual
					psolver->AssembleResidual(en, LM, fe, F);
					
					// do the biphasic stuff
					// TODO: I should only do this when the node is actually in contact
					if (sporo && mporo && ss.m_pme[ni])
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
void FETiedBiphasicInterface::ContactStiffness(FENLSolver* psolver)
{
	int i, j, k, l;
	vector<int> sLM, mLM, LM, en;
	double detJ[4], w[4], *Hs, Hm[4], pt[4], dpr[4], dps[4];
	double N[32];
	matrix ke;

	// get time step
	double dt = m_pfem->GetCurrentStep()->m_dt;
	
	// get the mesh
	FEMesh* pm = m_ss.GetMesh();
	
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
		FETiedBiphasicSurface& ss = (np == 0? m_ss : m_ms);
		FETiedBiphasicSurface& ms = (np == 0? m_ms : m_ss);
		
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
			
			// copy the LM vector
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
					if (mporo) for (k=0; k<nmeln; ++k) pm[k] = ms.GetMesh()->Node(me.m_node[k]).m_pt;
					
					// copy the LM vector
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
					vec3d dg = ss.m_dg[ni];
					
					// lagrange multiplier
					vec3d Lm = ss.m_Lmd[ni];
					
					// penalty 
					double eps = m_epsn*ss.m_epsn[ni];
					
					// contact traction
					vec3d t = Lm + dg*eps;
					
					// create the stiffness matrix
					ke.resize(ndof, ndof); ke.zero();
					
					// --- S O L I D - S O L I D   C O N T A C T ---
					
					// a. I-term
					//------------------------------------
					
					for (k=0; k<nseln; ++k) {
						for (l=0; l<nseln; ++l)
						{
							ke[ndpn*k    ][ndpn*l    ] += eps*Hs[k]*Hs[l]*detJ[j]*w[j];
							ke[ndpn*k + 1][ndpn*l + 1] += eps*Hs[k]*Hs[l]*detJ[j]*w[j];
							ke[ndpn*k + 2][ndpn*l + 2] += eps*Hs[k]*Hs[l]*detJ[j]*w[j];
						}
						for (l=0; l<nmeln; ++l)
						{
							ke[ndpn*k    ][ndpn*(nseln+l)    ] += -eps*Hs[k]*Hm[l]*detJ[j]*w[j];
							ke[ndpn*k + 1][ndpn*(nseln+l) + 1] += -eps*Hs[k]*Hm[l]*detJ[j]*w[j];
							ke[ndpn*k + 2][ndpn*(nseln+l) + 2] += -eps*Hs[k]*Hm[l]*detJ[j]*w[j];
						}
					}
					
					for (k=0; k<nmeln; ++k) {
						for (l=0; l<nseln; ++l)
						{
							ke[ndpn*(nseln+k)    ][ndpn*l    ] += -eps*Hm[k]*Hs[l]*detJ[j]*w[j];
							ke[ndpn*(nseln+k) + 1][ndpn*l + 1] += -eps*Hm[k]*Hs[l]*detJ[j]*w[j];
							ke[ndpn*(nseln+k) + 2][ndpn*l + 2] += -eps*Hm[k]*Hs[l]*detJ[j]*w[j];
						}
						for (l=0; l<nmeln; ++l)
						{
							ke[ndpn*(nseln+k)    ][ndpn*(nseln+l)    ] += eps*Hm[k]*Hm[l]*detJ[j]*w[j];
							ke[ndpn*(nseln+k) + 1][ndpn*(nseln+l) + 1] += eps*Hm[k]*Hm[l]*detJ[j]*w[j];
							ke[ndpn*(nseln+k) + 2][ndpn*(nseln+l) + 2] += eps*Hm[k]*Hm[l]*detJ[j]*w[j];
						}
					}
					
					// b. A-term
					//-------------------------------------
					
					double* Gr = se.Gr(j);
					double* Gs = se.Gs(j);
					vec3d gs[2];
					ss.CoBaseVectors(se, j, gs);
					
					vec3d as[4];
					mat3d As[4];
					for (l=0; l<nseln; ++l) {
						as[l] = nu ^ (gs[1]*Gr[l] - gs[0]*Gs[l]);
						As[l] = t & as[l];
					}
					
					if (!m_bsymm)
					{
						// non-symmetric
						for (k=0; k<nseln; ++k) {
							for (l=0; l<nseln; ++l)
							{
								ke[ndpn*k    ][ndpn*l    ] += Hs[k]*As[l](0,0)*w[j];
								ke[ndpn*k    ][ndpn*l + 1] += Hs[k]*As[l](0,1)*w[j];
								ke[ndpn*k    ][ndpn*l + 2] += Hs[k]*As[l](0,2)*w[j];

								ke[ndpn*k + 1][ndpn*l    ] += Hs[k]*As[l](1,0)*w[j];
								ke[ndpn*k + 1][ndpn*l + 1] += Hs[k]*As[l](1,1)*w[j];
								ke[ndpn*k + 1][ndpn*l + 2] += Hs[k]*As[l](1,2)*w[j];

								ke[ndpn*k + 2][ndpn*l    ] += Hs[k]*As[l](2,0)*w[j];
								ke[ndpn*k + 2][ndpn*l + 1] += Hs[k]*As[l](2,1)*w[j];
								ke[ndpn*k + 2][ndpn*l + 2] += Hs[k]*As[l](2,2)*w[j];
							}
						}
						
						for (k=0; k<nmeln; ++k) {
							for (l=0; l<nseln; ++l)
							{
								ke[ndpn*(nseln+k)    ][ndpn*l    ] += -Hm[k]*As[l](0,0)*w[j];
								ke[ndpn*(nseln+k)    ][ndpn*l + 1] += -Hm[k]*As[l](0,1)*w[j];
								ke[ndpn*(nseln+k)    ][ndpn*l + 2] += -Hm[k]*As[l](0,2)*w[j];

								ke[ndpn*(nseln+k) + 1][ndpn*l    ] += -Hm[k]*As[l](1,0)*w[j];
								ke[ndpn*(nseln+k) + 1][ndpn*l + 1] += -Hm[k]*As[l](1,1)*w[j];
								ke[ndpn*(nseln+k) + 1][ndpn*l + 2] += -Hm[k]*As[l](1,2)*w[j];

								ke[ndpn*(nseln+k) + 2][ndpn*l    ] += -Hm[k]*As[l](2,0)*w[j];
								ke[ndpn*(nseln+k) + 2][ndpn*l + 1] += -Hm[k]*As[l](2,1)*w[j];
								ke[ndpn*(nseln+k) + 2][ndpn*l + 2] += -Hm[k]*As[l](2,2)*w[j];
							}
						}
						
					}
					else 
					{
						// symmetric
						for (k=0; k<nseln; ++k) {
							for (l=0; l<nseln; ++l)
							{
								ke[ndpn*k    ][ndpn*l    ] += 0.5*(Hs[k]*As[l](0,0)+Hs[l]*As[k](0,0))*w[j];
								ke[ndpn*k    ][ndpn*l + 1] += 0.5*(Hs[k]*As[l](0,1)+Hs[l]*As[k](1,0))*w[j];
								ke[ndpn*k    ][ndpn*l + 2] += 0.5*(Hs[k]*As[l](0,2)+Hs[l]*As[k](2,0))*w[j];
								
								ke[ndpn*k + 1][ndpn*l    ] += 0.5*(Hs[k]*As[l](1,0)+Hs[l]*As[k](0,1))*w[j];
								ke[ndpn*k + 1][ndpn*l + 1] += 0.5*(Hs[k]*As[l](1,1)+Hs[l]*As[k](1,1))*w[j];
								ke[ndpn*k + 1][ndpn*l + 2] += 0.5*(Hs[k]*As[l](1,2)+Hs[l]*As[k](2,1))*w[j];
								
								ke[ndpn*k + 2][ndpn*l    ] += 0.5*(Hs[k]*As[l](2,0)+Hs[l]*As[k](0,2))*w[j];
								ke[ndpn*k + 2][ndpn*l + 1] += 0.5*(Hs[k]*As[l](2,1)+Hs[l]*As[k](1,2))*w[j];
								ke[ndpn*k + 2][ndpn*l + 2] += 0.5*(Hs[k]*As[l](2,2)+Hs[l]*As[k](2,2))*w[j];
							}
						}
						
						for (k=0; k<nmeln; ++k) {
							for (l=0; l<nseln; ++l)
							{
								ke[ndpn*(nseln+k)    ][ndpn*l    ] += -0.5*Hm[k]*As[l](0,0)*w[j];
								ke[ndpn*(nseln+k)    ][ndpn*l + 1] += -0.5*Hm[k]*As[l](0,1)*w[j];
								ke[ndpn*(nseln+k)    ][ndpn*l + 2] += -0.5*Hm[k]*As[l](0,2)*w[j];
								
								ke[ndpn*(nseln+k) + 1][ndpn*l    ] += -0.5*Hm[k]*As[l](1,0)*w[j];
								ke[ndpn*(nseln+k) + 1][ndpn*l + 1] += -0.5*Hm[k]*As[l](1,1)*w[j];
								ke[ndpn*(nseln+k) + 1][ndpn*l + 2] += -0.5*Hm[k]*As[l](1,2)*w[j];
								
								ke[ndpn*(nseln+k) + 2][ndpn*l    ] += -0.5*Hm[k]*As[l](2,0)*w[j];
								ke[ndpn*(nseln+k) + 2][ndpn*l + 1] += -0.5*Hm[k]*As[l](2,1)*w[j];
								ke[ndpn*(nseln+k) + 2][ndpn*l + 2] += -0.5*Hm[k]*As[l](2,2)*w[j];
							}
						}
						
						for (k=0; k<nseln; ++k) {
							for (l=0; l<nmeln; ++l)
							{
								ke[ndpn*k    ][ndpn*(nseln+l)    ] += -0.5*Hm[l]*As[k](0,0)*w[j];
								ke[ndpn*k    ][ndpn*(nseln+l) + 1] += -0.5*Hm[l]*As[k](1,0)*w[j];
								ke[ndpn*k    ][ndpn*(nseln+l) + 2] += -0.5*Hm[l]*As[k](2,0)*w[j];
								
								ke[ndpn*k + 1][ndpn*(nseln+l)    ] += -0.5*Hm[l]*As[k](0,1)*w[j];
								ke[ndpn*k + 1][ndpn*(nseln+l) + 1] += -0.5*Hm[l]*As[k](1,1)*w[j];
								ke[ndpn*k + 1][ndpn*(nseln+l) + 2] += -0.5*Hm[l]*As[k](2,1)*w[j];
								
								ke[ndpn*k + 2][ndpn*(nseln+l)    ] += -0.5*Hm[l]*As[k](0,2)*w[j];
								ke[ndpn*k + 2][ndpn*(nseln+l) + 1] += -0.5*Hm[l]*As[k](1,2)*w[j];
								ke[ndpn*k + 2][ndpn*(nseln+l) + 2] += -0.5*Hm[l]*As[k](2,2)*w[j];
							}
						}
					}

					
					// --- B I P H A S I C   S T I F F N E S S ---
					if (sporo && mporo)
					{
						double epsp = (ss.m_pme[ni]) ? m_epsp*ss.m_epsp[ni] : 0.;
						
						// --- S O L I D - P R E S S U R E   C O N T A C T ---
						
						// b. A-term
						//-------------------------------------

						double wn = ss.m_Lmp[ni] + epsp*ss.m_pg[ni];
						
						if (!m_bsymm)
						{
							// non-symmetric
							for (k=0; k<nseln; ++k)
								for (l=0; l<nseln; ++l) {
								{
									ke[4*k + 3][4*l  ] += dt*w[j]*wn*Hs[k]*as[l].x;
									ke[4*k + 3][4*l+1] += dt*w[j]*wn*Hs[k]*as[l].y;
									ke[4*k + 3][4*l+2] += dt*w[j]*wn*Hs[k]*as[l].z;
								}
							}
							for (k=0; k<nmeln; ++k)
								for (l=0; l<nseln; ++l) {
									{
										ke[4*(k+nseln) + 3][4*l  ] += -dt*w[j]*wn*Hm[k]*as[l].x;
										ke[4*(k+nseln) + 3][4*l+1] += -dt*w[j]*wn*Hm[k]*as[l].y;
										ke[4*(k+nseln) + 3][4*l+2] += -dt*w[j]*wn*Hm[k]*as[l].z;
									}
								}
						}
						else 
						{
							// symmetric
							for (k=0; k<nseln; ++k)
								for (l=0; l<nseln; ++l) {
									{
										ke[4*k + 3][4*l  ] += dt*w[j]*wn*0.5*(Hs[k]*as[l].x+Hs[l]*as[k].x);
										ke[4*k + 3][4*l+1] += dt*w[j]*wn*0.5*(Hs[k]*as[l].y+Hs[l]*as[k].y);
										ke[4*k + 3][4*l+2] += dt*w[j]*wn*0.5*(Hs[k]*as[l].z+Hs[l]*as[k].z);
									}
								}
							for (k=0; k<nmeln; ++k)
								for (l=0; l<nseln; ++l) {
									{
										ke[4*(k+nseln) + 3][4*l  ] += -dt*w[j]*wn*0.5*Hm[k]*as[l].x;
										ke[4*(k+nseln) + 3][4*l+1] += -dt*w[j]*wn*0.5*Hm[k]*as[l].y;
										ke[4*(k+nseln) + 3][4*l+2] += -dt*w[j]*wn*0.5*Hm[k]*as[l].z;
									}
								}
							for (k=0; k<nseln; ++k)
								for (l=0; l<nmeln; ++l) {
									{
										ke[4*k + 3][4*(nseln+l)  ] += -dt*w[j]*wn*0.5*Hm[l]*as[k].x;
										ke[4*k + 3][4*(nseln+l)+1] += -dt*w[j]*wn*0.5*Hm[l]*as[k].y;
										ke[4*k + 3][4*(nseln+l)+2] += -dt*w[j]*wn*0.5*Hm[l]*as[k].z;
									}
								}
						}

						
						// --- P R E S S U R E - P R E S S U R E   C O N T A C T ---
						
						for (k=0; k<nseln; ++k) {
							for (l=0; l<nseln; ++l)
								ke[4*k + 3][4*l+3] += -dt*epsp*w[j]*detJ[j]*Hs[k]*Hs[l];
							for (l=0; l<nmeln; ++l)
								ke[4*k + 3][4*(nseln+l)+3] += dt*epsp*w[j]*detJ[j]*Hs[k]*Hm[l];
						}
						
						for (k=0; k<nmeln; ++k) {
							for (l=0; l<nseln; ++l)
								ke[4*(nseln+k)+3][4*l + 3] += dt*epsp*w[j]*detJ[j]*Hm[k]*Hs[l];
							for (l=0; l<nmeln; ++l)
								ke[4*(nseln+k)+3][4*(nseln+l) + 3] += -dt*epsp*w[j]*detJ[j]*Hm[k]*Hm[l];
						}
						
					}
					
					// assemble the global stiffness
					psolver->AssembleStiffness(en, LM, ke);
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------
bool FETiedBiphasicInterface::Augment(int naug)
{
	// make sure we need to augment
	if (!m_blaugon) return true;
	
	int i;
	vec3d Ln;
	double Lp;
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
		m_ss.m_Lmd[i] = m_ss.m_Lmd[i] + m_ss.m_dg[i]*eps;
		
		normL1 += m_ss.m_Lmd[i]*m_ss.m_Lmd[i];
		
		if (m_ss.m_bporo) {
			Lp = 0;
			if (m_ss.m_pme[i]) {
				epsp = m_epsp*m_ss.m_epsp[i];
				Lp = m_ss.m_Lmp[i] + epsp*m_ss.m_pg[i];
				maxpg = max(maxpg,fabs(m_ss.m_pg[i]));
				normDP += m_ss.m_pg[i]*m_ss.m_pg[i];
			}
			m_ss.m_Lmp[i] = Lp;
		}
		
		maxgap = max(maxgap,sqrt(m_ss.m_dg[i]*m_ss.m_dg[i]));
	}	
	
	for (i=0; i<NM; ++i)
	{
		// update Lagrange multipliers on master surface
		eps = m_epsn*m_ms.m_epsn[i];
		m_ms.m_Lmd[i] = m_ms.m_Lmd[i] + m_ms.m_dg[i]*eps;
		
		normL1 += m_ms.m_Lmd[i]*m_ms.m_Lmd[i];
		
		if (m_ms.m_bporo) {
			Lp = 0;
			if (m_ms.m_pme[i]) {
				epsp = m_epsp*m_ms.m_epsp[i];
				Lp = m_ms.m_Lmp[i] + epsp*m_ms.m_pg[i];
				maxpg = max(maxpg,fabs(m_ms.m_pg[i]));
				normDP += m_ms.m_pg[i]*m_ms.m_pg[i];
			}
			m_ms.m_Lmp[i] = Lp;
		}
		
		maxgap = max(maxgap,sqrt(m_ms.m_dg[i]*m_ms.m_dg[i]));
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
void FETiedBiphasicInterface::Serialize(DumpFile &ar)
{
	// store contact data
	FEContactInterface::Serialize(ar);

	// store contact surface data
	m_ms.Serialize(ar);
	m_ss.Serialize(ar);
}

//-----------------------------------------------------------------------------

bool FETiedBiphasicInterface::PoroStatus(FEMesh& m, FESurfaceElement& el)
{
	bool status = false;
	
	// get the solid element this surface element belongs to
	FESolidElement* pe = dynamic_cast<FESolidElement*>(m.FindElementFromID(el.m_nelem));
	if (pe)
	{
		// get the material
		FEMaterial* pm = dynamic_cast<FEMaterial*>(m_pfem->GetMaterial(pe->GetMatID()));
		
		// see if this is a poro-elastic element
		FEBiphasic* biph = dynamic_cast<FEBiphasic*> (pm);
		if (biph) status = true;
	}
	
	return status;
}
