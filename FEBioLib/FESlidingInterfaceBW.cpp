#include "stdafx.h"
#include "FESlidingInterfaceBW.h"
#include "FESolidSolver.h"
#include "log.h"

//-----------------------------------------------------------------------------
// Define sliding interface parameters
BEGIN_PARAMETER_LIST(FESlidingInterfaceBW, FEContactInterface)
	ADD_PARAMETER(m_blaugon  , FE_PARAM_BOOL  , "laugon"             );
	ADD_PARAMETER(m_atol     , FE_PARAM_DOUBLE, "tolerance"          );
	ADD_PARAMETER(m_gtol     , FE_PARAM_DOUBLE, "gaptol"             );
	ADD_PARAMETER(m_epsn     , FE_PARAM_DOUBLE, "penalty"            );
	ADD_PARAMETER(m_bautopen , FE_PARAM_BOOL  , "auto_penalty"       );
	ADD_PARAMETER(m_btwo_pass, FE_PARAM_BOOL  , "two_pass"           );
	ADD_PARAMETER(m_knmult   , FE_PARAM_INT   , "knmult"             );
	ADD_PARAMETER(m_stol     , FE_PARAM_DOUBLE, "search_tol"         );
	ADD_PARAMETER(m_bsymm    , FE_PARAM_BOOL  , "symmetric_stiffness");
	ADD_PARAMETER(m_srad     , FE_PARAM_DOUBLE, "search_radius"      );
	ADD_PARAMETER(m_nsegup   , FE_PARAM_INT   , "seg_up"             );
	ADD_PARAMETER(m_btension , FE_PARAM_BOOL  , "tension"            );
	ADD_PARAMETER(m_naugmin  , FE_PARAM_INT   , "minaug"             );
	ADD_PARAMETER(m_naugmax  , FE_PARAM_INT   , "maxaug"             );
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
// FESlidingSurfaceBW
//-----------------------------------------------------------------------------

FESlidingSurfaceBW::FESlidingSurfaceBW(FEModel* pfem) : FEContactSurface(&pfem->GetMesh())
{ 
	m_pfem = pfem; 
}

//-----------------------------------------------------------------------------
bool FESlidingSurfaceBW::Init()
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
	m_gap.assign(nint, 0.0);
	m_nu.resize(nint);
	m_rs.resize(nint);
	m_Lmd.assign(nint, 0.0);
	m_pme.assign(nint, static_cast<FESurfaceElement*>(0));
	m_epsn.assign(nint, 1.0);
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
	
	return true;
}

//-----------------------------------------------------------------------------
void FESlidingSurfaceBW::ShallowCopy(FESlidingSurfaceBW &s)
{
	m_Lmd = s.m_Lmd;
	m_gap = s.m_gap;
	m_Ln = s.m_Ln;
	zero(m_pme);
}

//-----------------------------------------------------------------------------
//! This function calculates the node normal. Due to the piecewise continuity
//! of the surface elements this normal is not uniquely defined so in order to
//! obtain a unique normal the normal is averaged for each node over all the 
//! element normals at the node

void FESlidingSurfaceBW::UpdateNodeNormals()
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
void FESlidingSurfaceBW::Serialize(DumpFile& ar)
{
	// We can serialize the base-class data
	FEContactSurface::Serialize(ar);
	
	// And finally, we serialize the surface data
	if (ar.IsSaving())
	{
		ar << m_gap;
		ar << m_nu;
		ar << m_rs;
		ar << m_Lmd;
		ar << m_nei;
		ar << m_epsn;
		ar << m_nn;
		ar << m_Ln;
	}
	else
	{
		ar >> m_gap;
		ar >> m_nu;
		ar >> m_rs;
		ar >> m_Lmd;
		ar >> m_nei;
		ar >> m_epsn;
		ar >> m_nn;
		ar >> m_Ln;
	}
}

//-----------------------------------------------------------------------------
// FESlidingInterfaceBW
//-----------------------------------------------------------------------------

FESlidingInterfaceBW::FESlidingInterfaceBW(FEModel* pfem) : FEContactInterface(pfem), m_ss(pfem), m_ms(pfem)
{
	m_ntype = FE_CONTACT_SLIDINGBW;
	static int count = 1;
	m_nID = count++;
	
	// initial values
	m_knmult = 1;
	m_atol = 0.1;
	m_epsn = 1;
	m_btwo_pass = false;
	m_stol = 0.01;
	m_bsymm = true;
	m_srad = 1.0;
	m_gtol = -1;	// we use augmentation tolerance by default
	m_bautopen = false;
	m_btension = false;
	
	m_naugmin = 0;
	m_naugmax = 10;
	
	m_ss.SetSibling(&m_ms);
	m_ms.SetSibling(&m_ss);
}

//-----------------------------------------------------------------------------

FESlidingInterfaceBW::~FESlidingInterfaceBW()
{
}

//-----------------------------------------------------------------------------
bool FESlidingInterfaceBW::Init()
{
	// initialize surface data
	if (m_ss.Init() == false) return false;
	if (m_ms.Init() == false) return false;

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
		CalcAutoPenalty(m_ms);
	}
	
	// update sliding interface data
	Update(0);
	
	return true;
}

//-----------------------------------------------------------------------------
void FESlidingInterfaceBW::CalcAutoPenalty(FESlidingSurfaceBW& s)
{
	// get the mesh
	FEMesh& m = m_pfem->GetMesh();
	
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
void FESlidingInterfaceBW::ProjectSurface(FESlidingSurfaceBW& ss, FESlidingSurfaceBW& ms, bool bupseg)
{
	bool bfirst = true;
	
	FEMesh& mesh = m_pfem->GetMesh();
	FESurfaceElement* pme;
	vec3d r, nu;
	double rs[2];
	double Ln;
	
	double R = m_srad*mesh.GetBoundingBox().radius();
	
	// loop over all integration points
	int n = 0;
	for (int i=0; i<ss.Elements(); ++i)
	{
		FESurfaceElement& el = ss.Element(i);
		
		int ne = el.Nodes();
		int nint = el.GaussPoints();
		
		for (int j=0; j<nint; ++j, ++n)
		{
			// calculate the global position of the integration point
			r = ss.Local2Global(el, j);
			
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
			if (pme == 0 && bupseg) pme = ms.FindIntersection(r, nu, rs, bfirst, m_stol, m_srad);
			
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
				
				if ((g > R) || ((!m_btension) && (Ln < 0)) )
//				if ((g > R) || (Ln < 0) )
					ss.m_pme[n] = 0;
			}
			else
			{
				// the node is not in contact
				ss.m_Lmd[n] = 0;
				ss.m_gap[n] = 0;
			}
		}
	}
}

//-----------------------------------------------------------------------------

void FESlidingInterfaceBW::Update(int nsolve_iter)
{	
	static int naug = 0;
	static int biter = 0;
	
	FEModel& fem = *m_pfem;
	
	// get the iteration number
	// we need this number to see if we can do segment updates or not
	// also reset number of iterations after each augmentation
	FEAnalysis* pstep = fem.GetCurrentStep();
	if (pstep->m_psolver->m_niter == 0) {
		biter = 0;
		naug = pstep->m_psolver->m_naug;
	} else if (pstep->m_psolver->m_naug > naug) {
		biter = pstep->m_psolver->m_niter;
		naug = pstep->m_psolver->m_naug;
	}
	int niter = pstep->m_psolver->m_niter - biter;
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

	return;
}

//-----------------------------------------------------------------------------
void FESlidingInterfaceBW::ShallowCopy(FEContactInterface &ci)
{
	FESlidingInterfaceBW& si = dynamic_cast<FESlidingInterfaceBW&>(ci);
	m_ss.ShallowCopy(si.m_ss);
	m_ms.ShallowCopy(si.m_ms);
}

//-----------------------------------------------------------------------------
void FESlidingInterfaceBW::ContactForces(FEGlobalVector& R)
{
	int i, j, k;
	vector<int> sLM, mLM, LM, en;
	vector<double> fe;
	double detJ[4], w[4], *Hs, Hm[4];
	double N[24];

	FEModel& fem = *m_pfem;
	
	// get the mesh
	FEMesh* pm = m_ss.GetMesh();
	
	// loop over the nr of passes
	int npass = (m_btwo_pass?2:1);
	for (int np=0; np<npass; ++np)
	{
		// get slave and master surface
		FESlidingSurfaceBW& ss = (np == 0? m_ss : m_ms);
		FESlidingSurfaceBW& ms = (np == 0? m_ms : m_ss);
		
		// keep a running counter of integration points
		int ni = 0;
		
		// loop over all slave elements
		for (i=0; i<ss.Elements(); ++i)
		{
			// get the surface element
			FESurfaceElement& se = ss.Element(i);
			
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
					double g = ss.m_gap[ni];
					
					// lagrange multiplier
					double Lm = ss.m_Lmd[ni];
					
					// penalty 
					double eps = m_epsn*ss.m_epsn[ni];
					
					// contact traction
					double tn = Lm + eps*g;
					tn = m_btension ? tn : MBRACKET(tn);
//					tn = MBRACKET(tn);
					
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
					R.Assemble(en, LM, fe);
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------
void FESlidingInterfaceBW::ContactStiffness(FENLSolver* psolver)
{
	int i, j, k, l;
	vector<int> sLM, mLM, LM, en;
	double detJ[4], w[4], *Hs, Hm[4];
	double N[24];
	matrix ke;
	
	FEModel& fem = *m_pfem;
	
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
		FESlidingSurfaceBW& ss = (np == 0? m_ss : m_ms);
		FESlidingSurfaceBW& ms = (np == 0? m_ms : m_ss);
		
		// keep a running counter of the integration points
		int ni = 0;
		
		// loop over all slave elements
		for (i=0; i<ss.Elements(); ++i)
		{
			// get ths slave element
			FESurfaceElement& se = ss.Element(i);
			
			// get nr of nodes and integration points
			int nseln = se.Nodes();
			int nint = se.GaussPoints();
			
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
				
			}
			
			// loop over all integration points
			for (j=0; j<nint; ++j, ++ni)
			{
				// get the master element
				FESurfaceElement* pme = ss.m_pme[ni];
				if (pme)
				{
					FESurfaceElement& me = *pme;
					
					// get the nr of master nodes
					int nmeln = me.Nodes();
					
					// copy the LM vector
					ms.UnpackLM(me, mLM);
					
					// calculate degrees of freedom
					int ndpn = 3;
					int ndof = ndpn*(nseln + nmeln);
					
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
					tn = m_btension ? tn : MBRACKET(tn);
//					tn = MBRACKET(tn);
					
					//					double dtn = m_eps*HEAVYSIDE(Lm + eps*g);
					double dtn = m_btension ? eps : (tn > 0.? eps :0.);
//					double dtn = (tn > 0.? eps :0.);
					
					// create the stiffness matrix
					ke.resize(ndof, ndof); ke.zero();
					
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
					
					// assemble the global stiffness
					psolver->AssembleStiffness(en, LM, ke);
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------
void FESlidingInterfaceBW::UpdateContactPressures()
{
	int np, n, i, j, k;
	int npass = (m_btwo_pass?2:1);
	for (np=0; np<npass; ++np)
	{
		FESlidingSurfaceBW& ss = (np == 0? m_ss : m_ms);
		FESlidingSurfaceBW& ms = (np == 0? m_ms : m_ss);
		
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
				ss.m_Ln[ni] = m_btension ? (ss.m_Lmd[ni] + eps*gap) : MBRACKET(ss.m_Lmd[ni] + eps*gap);
//				ss.m_Ln[ni] = MBRACKET(ss.m_Lmd[ni] + eps*gap);
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
						ti[j] = m_btension ? (ms.m_Lmd[k] + m_epsn*ms.m_epsn[k]*ms.m_gap[k]) : 
						MBRACKET(ms.m_Lmd[k] + m_epsn*ms.m_epsn[k]*ms.m_gap[k]);
//						ti[j] = MBRACKET(ms.m_Lmd[k] + m_epsn*ms.m_epsn[k]*ms.m_gap[k]);
					}
					// project the data to the nodes
					double tn[4];
					pme->project_to_nodes(ti, tn);
					// now evaluate the traction at the intersection point
					double Ln = pme->eval(tn, ss.m_rs[ni][0], ss.m_rs[ni][1]);
					ss.m_Ln[ni] += (m_btension ? Ln : MBRACKET(Ln));
//					ss.m_Ln[ni] += MBRACKET(Ln);
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------
bool FESlidingInterfaceBW::Augment(int naug)
{
	// make sure we need to augment
	if (!m_blaugon) return true;
	
	int i;
	double Ln;
	bool bconv = true;
	
	int NS = m_ss.m_Lmd.size();
	int NM = m_ms.m_Lmd.size();
	
	// --- c a l c u l a t e   i n i t i a l   n o r m s ---
	// a. normal component
	double normL0 = 0;
	for (i=0; i<NS; ++i) normL0 += m_ss.m_Lmd[i]*m_ss.m_Lmd[i];
	for (i=0; i<NM; ++i) normL0 += m_ms.m_Lmd[i]*m_ms.m_Lmd[i];
	
	// b. gap component
	// (is calculated during update)
	double maxgap = 0;
	
	// update Lagrange multipliers
	double normL1 = 0, eps;
	for (i=0; i<NS; ++i)
	{
		// update Lagrange multipliers on slave surface
		eps = m_epsn*m_ss.m_epsn[i];
		Ln = m_ss.m_Lmd[i] + eps*m_ss.m_gap[i];
		m_ss.m_Lmd[i] = m_btension ? Ln : MBRACKET(Ln);
//		m_ss.m_Lmd[i] = MBRACKET(Ln);
		
		normL1 += m_ss.m_Lmd[i]*m_ss.m_Lmd[i];
		
		if (m_btension)
			maxgap = max(maxgap,fabs(m_ss.m_gap[i]));
		else if (Ln > 0) maxgap = max(maxgap,fabs(m_ss.m_gap[i]));
//		if (Ln > 0) maxgap = max(maxgap,fabs(m_ss.m_gap[i]));
	}	
	
	for (i=0; i<NM; ++i)
	{
		// update Lagrange multipliers on master surface
		eps = m_epsn*m_ms.m_epsn[i];
		Ln = m_ms.m_Lmd[i] + eps*m_ms.m_gap[i];
		m_ms.m_Lmd[i] = m_btension ? Ln : MBRACKET(Ln);
//		m_ms.m_Lmd[i] = MBRACKET(Ln);
		
		normL1 += m_ms.m_Lmd[i]*m_ms.m_Lmd[i];
		
		if (m_btension)
			maxgap = max(maxgap,fabs(m_ms.m_gap[i]));
		else if (Ln > 0) maxgap = max(maxgap,fabs(m_ms.m_gap[i]));
//		if (Ln > 0) maxgap = max(maxgap,fabs(m_ms.m_gap[i]));
	}
	
	// calculate relative norms
	double lnorm = (normL1 != 0 ? fabs((normL1 - normL0) / normL1) : fabs(normL1 - normL0)); 
	
	// check convergence
	if ((m_gtol > 0) && (maxgap > m_gtol)) bconv = false;
	if ((m_atol > 0) && (lnorm > m_atol)) bconv = false;
	if (naug < m_naugmin ) bconv = false;
	if (naug >= m_naugmax) bconv = true;
	
	clog.printf(" sliding interface # %d\n", m_nID);
	clog.printf("                        CURRENT        REQUIRED\n");
	clog.printf("    D multiplier : %15le", lnorm); if (m_atol > 0) clog.printf("%15le\n", m_atol); else clog.printf("       ***\n");
	
	clog.printf("    maximum gap  : %15le", maxgap);
	if (m_gtol > 0) clog.printf("%15le\n", m_gtol); else clog.printf("       ***\n");
	
	return bconv;
}

//-----------------------------------------------------------------------------
void FESlidingInterfaceBW::Serialize(DumpFile &ar)
{
	// serialize contact data
	FEContactInterface::Serialize(ar);
	
	// serialize contact surface data
	m_ms.Serialize(ar);
	m_ss.Serialize(ar);
}
