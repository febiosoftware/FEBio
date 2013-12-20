#include "stdafx.h"
#include "FESlidingInterfaceBW.h"
#include "FESolidSolver.h"
#include "FECore/log.h"

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
FESlidingSurfaceBW::Data::Data()
{
	m_gap = 0.0;
	m_Lmd = 0.0;
	m_Ln = 0.0;
	m_epsn = 1.0;
	m_nu = vec3d(0,0,0);
	m_rs = vec2d(0,0);
	m_pme = (FESurfaceElement*)0;
}

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
	
	// allocate data structure
	const int NE = Elements();
	m_Data.resize(NE);
	for (int i=0; i<Elements(); ++i)
	{
		FESurfaceElement& el = Element(i);
		int nint = el.GaussPoints();
		m_Data[i].resize(nint);
	}
	return true;
}

//-----------------------------------------------------------------------------
void FESlidingSurfaceBW::InitProjection()
{
	m_OT.Init();
}

//-----------------------------------------------------------------------------
//! \todo Originally, we only copied Lmd, gap, Ln and reset pme to zero.
//!       Need to check if this achieves the same
void FESlidingSurfaceBW::ShallowCopy(DumpStream& dmp, bool bsave)
{
	// And finally, we serialize the surface data
	if (bsave)
	{
		for (int i=0; i<(int) m_Data.size(); ++i)
		{
			vector<Data>& di = m_Data[i];
			int nint = (int) di.size();
			for (int j=0; j<nint; ++j)
			{
				Data& d = di[j];
				dmp << d.m_gap;
				dmp << d.m_nu;
				dmp << d.m_rs;
				dmp << d.m_Lmd;
				dmp << d.m_epsn;
				dmp << d.m_Ln;
			}
		}
	}
	else
	{
		for (int i=0; i<(int) m_Data.size(); ++i)
		{
			vector<Data>& di = m_Data[i];
			int nint = (int) di.size();
			for (int j=0; j<nint; ++j)
			{
				Data& d = di[j];
				dmp >> d.m_gap;
				dmp >> d.m_nu;
				dmp >> d.m_rs;
				dmp >> d.m_Lmd;
				dmp >> d.m_epsn;
				dmp >> d.m_Ln;
			}
		}
	}
}

//-----------------------------------------------------------------------------
void FESlidingSurfaceBW::Serialize(DumpFile& ar)
{
	// We can serialize the base-class data
	FEContactSurface::Serialize(ar);
	
	// And finally, we serialize the surface data
	if (ar.IsSaving())
	{
		for (int i=0; i<(int) m_Data.size(); ++i)
		{
			vector<Data>& di = m_Data[i];
			int nint = (int) di.size();
			for (int j=0; j<nint; ++j)
			{
				Data& d = di[j];
				ar << d.m_gap;
				ar << d.m_nu;
				ar << d.m_rs;
				ar << d.m_Lmd;
				ar << d.m_epsn;
				ar << d.m_Ln;
			}
		}
	}
	else
	{
		for (int i=0; i<(int) m_Data.size(); ++i)
		{
			vector<Data>& di = m_Data[i];
			int nint = (int) di.size();
			for (int j=0; j<nint; ++j)
			{
				Data& d = di[j];
				ar >> d.m_gap;
				ar >> d.m_nu;
				ar >> d.m_rs;
				ar >> d.m_Lmd;
				ar >> d.m_epsn;
				ar >> d.m_Ln;
			}
		}
	}
}

//-----------------------------------------------------------------------------
vec3d FESlidingSurfaceBW::GetContactForce()
{
	// initialize contact force
	vec3d f(0,0,0);
	
	// loop over all elements of the primary surface
	for (int n=0; n<Elements(); ++n)
	{
		FESurfaceElement& el = Element(n);
		int nint = el.GaussPoints();
		
		// evaluate the contact force for that element
		for (int i=0; i<nint; ++i) 
		{
			// get data for this integration point
			Data& data = m_Data[n][i];

			// get the base vectors
			vec3d g[2];
			CoBaseVectors(el, i, g);

			// normal (magnitude = area)
			vec3d n = g[0] ^ g[1];

			// gauss weight
			double w = el.GaussWeights()[i];

			// contact force
			f += n*(w*data.m_Ln);
		}
	}
	
	return f;
}

//-----------------------------------------------------------------------------
void FESlidingSurfaceBW::GetNodalContactGap(int nface, double* pg)
{
	FESurfaceElement& el = Element(nface);
	int ne = el.Nodes();
	int ni = el.GaussPoints();
	double gi[FEElement::MAX_INTPOINTS];
	for (int k=0; k<ni; ++k) gi[k] = m_Data[nface][k].m_gap;
	el.project_to_nodes(gi, pg);
}

//-----------------------------------------------------------------------------
void FESlidingSurfaceBW::GetNodalContactPressure(int nface, double* pg)
{
	FESurfaceElement& el = Element(nface);
	int ne = el.Nodes();
	int ni = el.GaussPoints();
	double ti[FEElement::MAX_INTPOINTS];
	for (int k=0; k<ni; ++k) ti[k] = m_Data[nface][k].m_Ln;
	el.project_to_nodes(ti, pg);
}

//-----------------------------------------------------------------------------
void FESlidingSurfaceBW::GetNodalContactTraction(int nface, vec3d* pt)
{
	FESurfaceElement& el = Element(nface);
	int ne = el.Nodes();
	int ni = el.GaussPoints();

	vec3d t;
	const int MFI = FEElement::MAX_INTPOINTS;
	double tix[MFI], tiy[MFI], tiz[MFI];
	for (int k=0; k<ni; ++k)
	{
		Data& pt = m_Data[nface][k];
		double gi = pt.m_gap;
		double Li = pt.m_Ln;
		vec3d  ti = pt.m_nu;
		if (gi > 0) t = ti*(Li); else t = vec3d(0,0,0);
		tix[k] = t.x; tiy[k] = t.y; tiz[k] = t.z;
	}

	// project traction to nodes
	const int MFN = FEElement::MAX_NODES;
	double tnx[MFN], tny[MFN], tnz[MFN];
	el.project_to_nodes(tix, tnx);
	el.project_to_nodes(tiy, tny);
	el.project_to_nodes(tiz, tnz);

	// store data
	for (int k=0; k<ne; ++k)
	{
		pt[k].x = tnx[k];
		pt[k].y = tny[k];
		pt[k].z = tnz[k];
	}
}

//-----------------------------------------------------------------------------
// FESlidingInterfaceBW
//-----------------------------------------------------------------------------

FESlidingInterfaceBW::FESlidingInterfaceBW(FEModel* pfem) : FEContactInterface(pfem), m_ss(pfem), m_ms(pfem)
{
	static int count = 1;
	m_nID = count++;
	
	// initial values
	m_knmult = 1;
	m_atol = 0.1;
	m_epsn = 1;
	m_btwo_pass = false;
	m_gtol = 0;
	m_stol = 0.01;
	m_bsymm = true;
	m_srad = 1.0;
	m_nsegup = 0;
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
	
	return true;
}

//-----------------------------------------------------------------------------
void FESlidingInterfaceBW::Activate()
{
	// don't forget to call the base class!
	FEContactInterface::Activate();

	// this contact implementation requires a non-symmetric stiffness matrix
	// so inform the FEM class
	if (!m_bsymm) 
	{
		// request a non-symmetric stiffness matrix
		FESolver* psolver = GetFEModel()->GetCurrentStep()->m_psolver;
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
}

//-----------------------------------------------------------------------------
//! build the matrix profile for use in the stiffness matrix
void FESlidingInterfaceBW::BuildMatrixProfile(FEStiffnessMatrix& K)
{
	FEMesh& mesh = GetFEModel()->GetMesh();

	vector<int> lm(6*FEElement::MAX_NODES*2);
					
	int npass = (m_btwo_pass?2:1);
	for (int np=0; np<npass; ++np)
	{
		FESlidingSurfaceBW& ss = (np == 0? m_ss : m_ms);
		FESlidingSurfaceBW& ms = (np == 0? m_ms : m_ss);
						
		int k, l;
		for (int j=0; j<ss.Elements(); ++j)
		{
			FESurfaceElement& se = ss.Element(j);
			int nint = se.GaussPoints();
			int* sn = &se.m_node[0];
			for (k=0; k<nint; ++k)
			{
				FESurfaceElement* pe = ss.m_Data[j][k].m_pme;
				if (pe != 0)
				{
					FESurfaceElement& me = dynamic_cast<FESurfaceElement&> (*pe);
					int* mn = &me.m_node[0];
									
					assign(lm, -1);
									
					int nseln = se.Nodes();
					int nmeln = me.Nodes();
									
					for (l=0; l<nseln; ++l)
					{
						int* id = mesh.Node(sn[l]).m_ID;
						lm[6*l  ] = id[DOF_X];
						lm[6*l+1] = id[DOF_Y];
						lm[6*l+2] = id[DOF_Z];
						lm[6*l+3] = id[DOF_RU];
						lm[6*l+4] = id[DOF_RV];
						lm[6*l+5] = id[DOF_RW];
					}
									
					for (l=0; l<nmeln; ++l)
					{
						int* id = mesh.Node(mn[l]).m_ID;
						lm[6*(l+nseln)  ] = id[DOF_X];
						lm[6*(l+nseln)+1] = id[DOF_Y];
						lm[6*(l+nseln)+2] = id[DOF_Z];
						lm[6*(l+nseln)+3] = id[DOF_RU];
						lm[6*(l+nseln)+4] = id[DOF_RV];
						lm[6*(l+nseln)+5] = id[DOF_RW];
					}
									
					K.build_add(lm);
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------
void FESlidingInterfaceBW::CalcAutoPenalty(FESlidingSurfaceBW& s)
{
	// get the mesh
	FEMesh& m = GetFEModel()->GetMesh();
	
	// loop over all surface elements
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
		double E = AutoPenalty(el, s);
		
		// calculate penalty
		double eps = E*A/V;
		
		// assign to integation points of surface element
		int nint = el.GaussPoints();
		for (int j=0; j<nint; ++j)
		{
			FESlidingSurfaceBW::Data& data = s.m_Data[i][j];
			data.m_epsn = eps;
		}
	}
}

//-----------------------------------------------------------------------------
void FESlidingInterfaceBW::ProjectSurface(FESlidingSurfaceBW& ss, FESlidingSurfaceBW& ms, bool bupseg)
{
	FEMesh& mesh = GetFEModel()->GetMesh();
	double R = m_srad*mesh.GetBoundingBox().radius();

	// initialize projection data
	ms.InitProjection();
	bool binit = false;
	
	// loop over all integration points
	#pragma omp parallel for shared(mesh, R, binit, ss, ms, bupseg)
	for (int i=0; i<ss.Elements(); ++i)
	{
		FESurfaceElement& el = ss.Element(i);
		
		int ne = el.Nodes();
		int nint = el.GaussPoints();
		
		for (int j=0; j<nint; ++j)
		{
			// get the integration point data
			FESlidingSurfaceBW::Data& data = ss.m_Data[i][j];

			// calculate the global position of the integration point
			vec3d r = ss.Local2Global(el, j);
			
			// calculate the normal at this integration point
			vec3d nu = ss.SurfaceNormal(el, j);
			
			// first see if the old intersected face is still good enough
			FESurfaceElement* pme = data.m_pme;
			double rs[2] = {0,0};
			if (pme)
			{
				double g;
				
				// see if the ray intersects this element
				if (ms.Intersect(*pme, r, nu, rs, g, m_stol))
				{
					data.m_rs[0] = rs[0];
					data.m_rs[1] = rs[1];
				}
				else
				{
					pme = 0;
				}
			}
			
			// find the intersection point with the master surface
			if (pme == 0 && bupseg) pme = ms.FindIntersection(r, nu, rs, binit, m_stol, R);
			
			data.m_pme = pme;
			data.m_nu = nu;
			data.m_rs[0] = rs[0];
			data.m_rs[1] = rs[1];
			if (pme)
			{
				// the node could potentially be in contact
				// find the global location of the intersection point
				vec3d q = ms.Local2Global(*pme, rs[0], rs[1]);
				
				// calculate the gap function
				// NOTE: this has the opposite sign compared
				// to Gerard's notes.
				double g = nu*(r - q);
				
				double eps = m_epsn*data.m_epsn;
				
				double Ln = data.m_Lmd + eps*g;
				
				data.m_gap = (g <= R? g : 0);
				
				if ((g > R) || ((!m_btension) && (Ln < 0)) ) {
//				if ((g > R) || (Ln < 0) )
					data.m_pme = 0;
					data.m_gap = 0;
				}
			}
			else
			{
				// the node is not in contact
				data.m_Lmd = 0;
				data.m_gap = 0;
			}
		}
	}
}

//-----------------------------------------------------------------------------

void FESlidingInterfaceBW::Update(int nsolve_iter)
{	
	static int naug = 0;
	static int biter = 0;
	
	FEModel& fem = *GetFEModel();
	
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
void FESlidingInterfaceBW::ShallowCopy(DumpStream& dmp, bool bsave)
{
	m_ss.ShallowCopy(dmp, bsave);
	m_ms.ShallowCopy(dmp, bsave);
}

//-----------------------------------------------------------------------------
void FESlidingInterfaceBW::ContactForces(FEGlobalVector& R)
{
	int i, j, k;
	vector<int> sLM, mLM, LM, en;
	vector<double> fe;
	const int MN = FEElement::MAX_NODES;
	double detJ[MN], w[MN], *Hs, Hm[MN];
	double N[MN*6];

	// loop over the nr of passes
	int npass = (m_btwo_pass?2:1);
	for (int np=0; np<npass; ++np)
	{
		// get slave and master surface
		FESlidingSurfaceBW& ss = (np == 0? m_ss : m_ms);
		FESlidingSurfaceBW& ms = (np == 0? m_ms : m_ss);
		
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
			for (j=0; j<nint; ++j)
			{
				// get integration point data
				FESlidingSurfaceBW::Data& data = ss.m_Data[i][j];

				// get the master element
				FESurfaceElement* pme = data.m_pme;
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
					double r = data.m_rs[0];
					double s = data.m_rs[1];
					me.shape_fnc(Hm, r, s);
					
					// get normal vector
					vec3d nu = data.m_nu;
					
					// gap function
					double g = data.m_gap;
					
					// lagrange multiplier
					double Lm = data.m_Lmd;
					
					// penalty 
					double eps = m_epsn*data.m_epsn;
					
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
void FESlidingInterfaceBW::ContactStiffness(FESolver* psolver)
{
	int i, j, k, l;
	vector<int> sLM, mLM, LM, en;
	const int MN = FEElement::MAX_NODES;
	double detJ[MN], w[MN], *Hs, Hm[MN];
	double N[MN*6];
	matrix ke;
	
	FEModel& fem = *GetFEModel();
	
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
			felog.printf("Higher order stiffness terms included.\n");
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
			for (j=0; j<nint; ++j)
			{
				// get integration point data
				FESlidingSurfaceBW::Data& data = ss.m_Data[i][j];

				// get the master element
				FESurfaceElement* pme = data.m_pme;
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
					double r = data.m_rs[0];
					double s = data.m_rs[1];
					me.shape_fnc(Hm, r, s);
					
					// get slave normal vector
					vec3d nu = data.m_nu;
					
					// gap function
					double g = data.m_gap;
					
					// lagrange multiplier
					double Lm = data.m_Lmd;
					
					// penalty 
					double eps = m_epsn*data.m_epsn;
					
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
					mat3d As[MN];
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
					
					double Hmr[MN], Hms[MN];
					me.shape_deriv(Hmr, Hms, r, s);
					vec3d mm[MN];
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
	int npass = (m_btwo_pass?2:1);
	const int MN = FEElement::MAX_NODES;
	for (int np=0; np<npass; ++np)
	{
		FESlidingSurfaceBW& ss = (np == 0? m_ss : m_ms);
		FESlidingSurfaceBW& ms = (np == 0? m_ms : m_ss);
		
		// loop over all elements of the primary surface
		for (int n=0; n<ss.Elements(); ++n)
		{
			FESurfaceElement& el = ss.Element(n);
			int nint = el.GaussPoints();
			
			// get the normal tractions at the integration points
			double gap, eps;
			for (int i=0; i<nint; ++i) 
			{
				// get integration point data
				FESlidingSurfaceBW::Data& data = ss.m_Data[n][i];

				gap = data.m_gap;
				eps = m_epsn*data.m_epsn;
				data.m_Ln = m_btension ? (data.m_Lmd + eps*gap) : MBRACKET(data.m_Lmd + eps*gap);
//				data.m_Ln = MBRACKET(data.m_Lmd + eps*gap);
				FESurfaceElement* pme = data.m_pme;
				if (m_btwo_pass && pme)
				{
					// get master element data
					vector<FESlidingSurfaceBW::Data>& mdv = ms.m_Data[pme->m_lid];
					int mint = pme->GaussPoints();
					double ti[MN];
					for (int j=0; j<mint; ++j) 
					{
						FESlidingSurfaceBW::Data& md = mdv[j];
						gap = md.m_gap;
						eps = m_epsn*md.m_epsn;
						ti[j] = m_btension ? (md.m_Lmd + m_epsn*md.m_epsn*md.m_gap) : 
						MBRACKET(md.m_Lmd + m_epsn*md.m_epsn*md.m_gap);
//						ti[j] = MBRACKET(md.m_Lmd + m_epsn*md.m_epsn*md.m_gap);
					}
					// project the data to the nodes
					double tn[MN];
					pme->project_to_nodes(ti, tn);
					// now evaluate the traction at the intersection point
					double Ln = pme->eval(tn, data.m_rs[0], data.m_rs[1]);
					data.m_Ln += (m_btension ? Ln : MBRACKET(Ln));
//					data.m_Ln += MBRACKET(Ln);
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
	
	double Ln;
	bool bconv = true;

	int NS = (int) m_ss.m_Data.size();
	int NM = (int) m_ms.m_Data.size();

	// --- c a l c u l a t e   i n i t i a l   n o r m s ---
	// a. normal component
	double normL0 = 0;
	for (int i=0; i<NS; ++i)
	{
		vector<FESlidingSurfaceBW::Data>& sd = m_ss.m_Data[i];
		for (int j=0; j<(int)sd.size(); ++j)
		{
			FESlidingSurfaceBW::Data& ds = sd[j];
			normL0 += ds.m_Lmd*ds.m_Lmd;
		}
	}
	for (int i=0; i<NM; ++i)
	{
		vector<FESlidingSurfaceBW::Data>& md = m_ms.m_Data[i];
		for (int j=0; j<(int)md.size(); ++j)
		{
			FESlidingSurfaceBW::Data& dm = md[j];
			normL0 += dm.m_Lmd*dm.m_Lmd;
		}
	}
	
	// b. gap component
	// (is calculated during update)
	double maxgap = 0;
	
	// update Lagrange multipliers
	double normL1 = 0, eps;
	for (int i=0; i<NS; ++i)
	{
		vector<FESlidingSurfaceBW::Data>& sd = m_ss.m_Data[i];
		for (int j=0; j<(int)sd.size(); ++j)
		{
			FESlidingSurfaceBW::Data& ds = sd[j];

			// update Lagrange multipliers on slave surface
			eps = m_epsn*ds.m_epsn;
			Ln = ds.m_Lmd + eps*ds.m_gap;
			ds.m_Lmd = m_btension ? Ln : MBRACKET(Ln);
//			ds.m_Lmd = MBRACKET(Ln);
		
			normL1 += ds.m_Lmd*ds.m_Lmd;
		
			if (m_btension)
				maxgap = max(maxgap,fabs(ds.m_gap));
			else if (Ln > 0) maxgap = max(maxgap,fabs(ds.m_gap));
//			if (Ln > 0) maxgap = max(maxgap,fabs(ds.m_gap));
		}
	}	
	
	for (int i=0; i<NM; ++i)
	{
		vector<FESlidingSurfaceBW::Data>& md = m_ms.m_Data[i];
		for (int j=0; j<(int)md.size(); ++j)
		{
			FESlidingSurfaceBW::Data& dm = md[j];

			// update Lagrange multipliers on master surface
			eps = m_epsn*dm.m_epsn;
			Ln = dm.m_Lmd + eps*dm.m_gap;
			dm.m_Lmd = m_btension ? Ln : MBRACKET(Ln);
//			dm.m_Lmd = MBRACKET(Ln);
		
			normL1 += dm.m_Lmd*dm.m_Lmd;
		
			if (m_btension)
				maxgap = max(maxgap,fabs(dm.m_gap));
			else if (Ln > 0) maxgap = max(maxgap,fabs(dm.m_gap));
//			if (Ln > 0) maxgap = max(maxgap,fabs(dm.m_gap));
		}
	}
	
	// calculate relative norms
	double lnorm = (normL1 != 0 ? fabs((normL1 - normL0) / normL1) : fabs(normL1 - normL0)); 
	
	// check convergence
	if ((m_gtol > 0) && (maxgap > m_gtol)) bconv = false;
	if ((m_atol > 0) && (lnorm > m_atol)) bconv = false;
	if (naug < m_naugmin ) bconv = false;
	if (naug >= m_naugmax) bconv = true;
	
	felog.printf(" sliding interface # %d\n", m_nID);
	felog.printf("                        CURRENT        REQUIRED\n");
	felog.printf("    D multiplier : %15le", lnorm); if (m_atol > 0) felog.printf("%15le\n", m_atol); else felog.printf("       ***\n");
	
	felog.printf("    maximum gap  : %15le", maxgap);
	if (m_gtol > 0) felog.printf("%15le\n", m_gtol); else felog.printf("       ***\n");
	
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
