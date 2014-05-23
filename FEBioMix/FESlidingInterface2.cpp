#include "FESlidingInterface2.h"
#include "FEBiphasic.h"
#include "FEBioMech/FEStiffnessMatrix.h"
#include "FECore/FEModel.h"
#include "FECore/FENormalProjection.h"
#include "FECore/log.h"

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
	ADD_PARAMETER(m_naugmin  , FE_PARAM_INT   , "minaug"             );
	ADD_PARAMETER(m_naugmax  , FE_PARAM_INT   , "maxaug"             );
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FESlidingSurface2::Data::Data()
{
	m_gap  = 0.0;
	m_Lmd  = 0.0;
	m_Lmp  = 0.0;
	m_epsn = 1.0;
	m_epsp = 1.0;
	m_Ln   = 0.0;
	m_pg   = 0.0;
	m_nu   = vec3d(0,0,0);
	m_rs   = vec2d(0,0);
	m_pme  = (FESurfaceElement*)0;
}

//-----------------------------------------------------------------------------
// FESlidingSurface2
//-----------------------------------------------------------------------------

FESlidingSurface2::FESlidingSurface2(FEModel* pfem) : FEBiphasicContactSurface(&pfem->GetMesh())
{ 
	m_bporo = false;
	m_pfem = pfem; 
}

//-----------------------------------------------------------------------------
bool FESlidingSurface2::Init()
{
	// initialize surface data first
	if (FEContactSurface::Init() == false) return false;

	// allocate data structures
	int NE = Elements();
	m_Data.resize(NE);
	for (int i=0; i<NE; ++i)
	{
		FESurfaceElement& el = Element(i);
		int nint = el.GaussPoints();
		m_Data[i].resize(nint);
	}

	// allocate node normals
	m_nn.assign(Nodes(), 0);

	// determine biphasic status
	m_poro.resize(Elements(),false);
	for (int i=0; i<Elements(); ++i)
	{
		// get the surface element
		FESurfaceElement& se = Element(i);
		
		// get the solid element this surface element belongs to
		FESolidElement* pe = dynamic_cast<FESolidElement*>(m_pMesh->FindElementFromID(se.m_nelem));
		if (pe)
		{
			// get the material
			FEMaterial* pm = m_pfem->GetMaterial(pe->GetMatID());
			
			// see if this is a poro-elastic element
			FEBiphasic* biph = dynamic_cast<FEBiphasic*> (pm);
			if (biph) {
				m_poro[i] = true;
				m_bporo = true;
			}
		}
	}

	return true;
}

//-----------------------------------------------------------------------------
//! This function calculates the node normal. Due to the piecewise continuity
//! of the surface elements this normal is not uniquely defined so in order to
//! obtain a unique normal the normal is averaged for each node over all the 
//! element normals at the node

void FESlidingSurface2::UpdateNodeNormals()
{
	const int MN = FEElement::MAX_NODES;
	vec3d y[MN];

	// zero nodal normals
	zero(m_nn);

	// loop over all elements
	for (int i=0; i<Elements(); ++i)
	{
		FESurfaceElement& el = Element(i);
		int ne = el.Nodes();

		// get the nodal coordinates
		for (int j=0; j<ne; ++j) y[j] = Node(el.m_lnode[j]).m_rt;

		// calculate the normals
		for (int j=0; j<ne; ++j)
		{
			int jp1 = (j+1)%ne;
			int jm1 = (j+ne-1)%ne;
			vec3d n = (y[jp1] - y[j]) ^ (y[jm1] - y[j]);
			m_nn[el.m_lnode[j]] += n;
		}
	}

	// normalize all vectors
	const int N = Nodes();
	for (int i=0; i<N; ++i) m_nn[i].unit();
}

//-----------------------------------------------------------------------------
vec3d FESlidingSurface2::GetContactForce()
{
	int n, i;
	
	// initialize contact force
	vec3d f(0,0,0);
	
	// loop over all elements of the surface
	for (n=0; n<Elements(); ++n)
	{
		FESurfaceElement& el = Element(n);
		int nint = el.GaussPoints();
		
		// evaluate the contact force for that element
		for (i=0; i<nint; ++i)
		{
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
double FESlidingSurface2::GetContactArea()
{
	// initialize contact area
	double a = 0;
	
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
            double s = (data.m_Ln > 0) ? 1 : 0;
            
			// get the base vectors
			vec3d g[2];
			CoBaseVectors(el, i, g);
            
			// normal (magnitude = area)
			vec3d n = g[0] ^ g[1];
            
			// gauss weight
			double w = el.GaussWeights()[i];
            
			// contact force
			a += n.norm()*(w*s);
		}
	}
	
	return a;
}

//-----------------------------------------------------------------------------
vec3d FESlidingSurface2::GetFluidForce()
{
	int n, i;
	const int MN = FEElement::MAX_NODES;
	double pn[MN];
	
	// initialize contact force
	vec3d f(0,0,0);
	
	// loop over all elements of the surface
	for (n=0; n<Elements(); ++n)
	{
		FESurfaceElement& el = Element(n);
		int nseln = el.Nodes();
        
		// nodal pressures
		for (i=0; i<nseln; ++i) pn[i] = GetMesh()->Node(el.m_node[i]).m_pt;
		
		int nint = el.GaussPoints();
		
		// evaluate the fluid force for that element
		for (i=0; i<nint; ++i)
		{
			// get the base vectors
			vec3d g[2];
			CoBaseVectors(el, i, g);
			// normal (magnitude = area)
			vec3d n = g[0] ^ g[1];
			// gauss weight
			double w = el.GaussWeights()[i];
			// fluid pressure
			double p = el.eval(pn, i);
			// contact force
			f += n*(w*p);
		}
	}
	
	return f;
}

//-----------------------------------------------------------------------------
void FESlidingSurface2::ShallowCopy(DumpStream& dmp, bool bsave)
{
	if (bsave)
	{
		dmp << m_bporo;

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
				dmp << d.m_Lmp;
				dmp << d.m_epsn;
				dmp << d.m_epsp;
				dmp << d.m_pg;
				dmp << d.m_Ln;
			}
		}
		dmp << m_poro;
		dmp << m_nn;
	}
	else
	{
		dmp >> m_bporo;

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
				dmp >> d.m_Lmp;
				dmp >> d.m_epsn;
				dmp >> d.m_epsp;
				dmp >> d.m_pg;
				dmp >> d.m_Ln;
			}
		}
		dmp >> m_poro;
		dmp >> m_nn;

		// reset element pointers
		for (int i=0; i<Elements(); ++i)
		{
			vector<Data>& di = m_Data[i];
			int n = (int) di.size();
			for (int j=0; j<n; ++j) di[j].m_pme = 0;
		}
	}
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
				ar << d.m_Lmp;
				ar << d.m_epsn;
				ar << d.m_epsp;
				ar << d.m_pg;
				ar << d.m_Ln;
			}
		}
		ar << m_poro;
		ar << m_nn;
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
				ar >> d.m_Lmp;
				ar >> d.m_epsn;
				ar >> d.m_epsp;
				ar >> d.m_pg;
				ar >> d.m_Ln;
			}
		}
		ar >> m_poro;
		ar >> m_nn;
	}
}

//-----------------------------------------------------------------------------
void FESlidingSurface2::GetNodalContactGap(int nface, double* pg)
{
	FESurfaceElement& el = Element(nface);
	int ni = el.GaussPoints();
	double gi[FEElement::MAX_INTPOINTS];
	for (int k=0; k<ni; ++k) gi[k] = m_Data[nface][k].m_gap;
	el.project_to_nodes(gi, pg);
}

//-----------------------------------------------------------------------------
void FESlidingSurface2::GetNodalContactPressure(int nface, double* pg)
{
	FESurfaceElement& el = Element(nface);
	int ni = el.GaussPoints();
	double ti[FEElement::MAX_INTPOINTS];
	for (int k=0; k<ni; ++k) ti[k] = m_Data[nface][k].m_Ln;
	el.project_to_nodes(ti, pg);
}

//-----------------------------------------------------------------------------
void FESlidingSurface2::GetNodalPressureGap(int nface, double* pg)
{
	FESurfaceElement& el = Element(nface);
	int ni = el.GaussPoints();
	double gi[FEElement::MAX_INTPOINTS];
	for (int k=0; k<ni; ++k) gi[k] = m_Data[nface][k].m_pg;
	el.project_to_nodes(gi, pg);
}

//-----------------------------------------------------------------------------
void FESlidingSurface2::GetNodalContactTraction(int nface, vec3d* pt)
{
	FESurfaceElement& el = Element(nface);
	int ne = el.Nodes();
	int ni = el.GaussPoints();

	vec3d t;
	const int MFI = FEElement::MAX_INTPOINTS;
	double tix[MFI], tiy[MFI], tiz[MFI];
	for (int k=0; k<ni; ++k)
	{
		FESlidingSurface2::Data& pt = m_Data[nface][k];
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
// FESlidingInterface2
//-----------------------------------------------------------------------------

FESlidingInterface2::FESlidingInterface2(FEModel* pfem) : FEContactInterface(pfem), m_ss(pfem), m_ms(pfem)
{
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
	m_gtol = 0;
	m_ptol = 0;
	m_nsegup = 0;
	m_bautopen = false;

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
bool FESlidingInterface2::Init()
{
	// initialize surface data
	if (m_ss.Init() == false) return false;
	if (m_ms.Init() == false) return false;
	
	// this contact implementation requires a non-symmetric stiffness matrix
	// so inform the FEM class
	if (!m_bsymm) 
	{
		// request a non-symmetric stiffness matrix
		FESolver* psolver = GetFEModel()->GetCurrentStep()->m_psolver;
		psolver->m_bsymm = false;
	}

	return true;
}

//-----------------------------------------------------------------------------
void FESlidingInterface2::BuildMatrixProfile(FEStiffnessMatrix& K)
{
	FEMesh& mesh = GetFEModel()->GetMesh();

	vector<int> lm(7*FEElement::MAX_NODES*2);

	int npass = (m_btwo_pass?2:1);
	for (int np=0; np<npass; ++np)
	{
		FESlidingSurface2& ss = (np == 0? m_ss : m_ms);
		FESlidingSurface2& ms = (np == 0? m_ms : m_ss);

		int k, l;
		for (int j=0; j<ss.Elements(); ++j)
		{
			FESurfaceElement& se = ss.Element(j);
			int nint = se.GaussPoints();
			int* sn = &se.m_node[0];
			for (k=0; k<nint; ++k)
			{
				FESlidingSurface2::Data& pt = ss.m_Data[j][k];
				FESurfaceElement* pe = pt.m_pme;
				if (pe != 0)
				{
					FESurfaceElement& me = *pe;
					int* mn = &me.m_node[0];

					assign(lm, -1);

					int nseln = se.Nodes();
					int nmeln = me.Nodes();

					for (l=0; l<nseln; ++l)
					{
						vector<int>& id = mesh.Node(sn[l]).m_ID;
						lm[7*l  ] = id[DOF_X];
						lm[7*l+1] = id[DOF_Y];
						lm[7*l+2] = id[DOF_Z];
						lm[7*l+3] = id[DOF_P];
						lm[7*l+4] = id[DOF_RU];
						lm[7*l+5] = id[DOF_RV];
						lm[7*l+6] = id[DOF_RW];
					}

					for (l=0; l<nmeln; ++l)
					{
						vector<int>& id = mesh.Node(mn[l]).m_ID;
						lm[7*(l+nseln)  ] = id[DOF_X];
						lm[7*(l+nseln)+1] = id[DOF_Y];
						lm[7*(l+nseln)+2] = id[DOF_Z];
						lm[7*(l+nseln)+3] = id[DOF_P];
						lm[7*(l+nseln)+4] = id[DOF_RU];
						lm[7*(l+nseln)+5] = id[DOF_RV];
						lm[7*(l+nseln)+6] = id[DOF_RW];
					}

					K.build_add(lm);
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! This function is called during the initialization
void FESlidingInterface2::Activate()
{
	// don't forget to call base member
	FEContactInterface::Activate();

	// calculate the penalty
	if (m_bautopen) 
	{
		CalcAutoPenalty(m_ss);
		CalcAutoPenalty(m_ms);
		CalcAutoPressurePenalty(m_ss);
		CalcAutoPressurePenalty(m_ms);
	}

	// update sliding interface data
	Update(0);
}

//-----------------------------------------------------------------------------
void FESlidingInterface2::CalcAutoPenalty(FESlidingSurface2& s)
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
		double K = AutoPenalty(el, s);

		// calculate penalty
		double eps = K*A/V;

		// assign to integation points of surface element
		int nint = el.GaussPoints();
		for (int j=0; j<nint; ++j)
		{
			FESlidingSurface2::Data& pt = s.m_Data[i][j];
			pt.m_epsn = eps;
		}
	}
}

//-----------------------------------------------------------------------------
void FESlidingInterface2::CalcAutoPressurePenalty(FESlidingSurface2& s)
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
		double k = AutoPressurePenalty(el, s);

		// calculate penalty
		double eps = k*A/V;

		// assign to integation points of surface element
		int nint = el.GaussPoints();
		for (int j=0; j<nint; ++j)
		{
			FESlidingSurface2::Data& pt = s.m_Data[i][j];
			pt.m_epsp = eps;
		}
	}
}

//-----------------------------------------------------------------------------

double FESlidingInterface2::AutoPressurePenalty(FESurfaceElement& el, FESlidingSurface2& s)
{
	// get the mesh
	FEMesh& m = GetFEModel()->GetMesh();

	// evaluate element surface normal at parametric center
	vec3d t[2];
	s.CoBaseVectors0(el, 0, 0, t);
	vec3d n = t[0] ^ t[1];
	n.unit();

	double eps = 0;

	// get the solid element this surface element belongs to
	FESolidElement* pe = dynamic_cast<FESolidElement*>(m.FindElementFromID(el.m_nelem));
	if (pe)
	{
		// get the material
		FEMaterial* pm = GetFEModel()->GetMaterial(pe->GetMatID());

		// see if this is a poro-elastic element
		FEBiphasic* biph = dynamic_cast<FEBiphasic*> (pm);
		if (biph)
		{
			// get a material point
			FEMaterialPoint& mp = *pe->GetMaterialPoint(0);
			FEElasticMaterialPoint& ept = *(mp.ExtractData<FEElasticMaterialPoint>());

			// setup the material point
			ept.m_F = mat3dd(1.0);
			ept.m_J = 1;
			ept.m_s.zero();

			// if this is a poroelastic element, then get the permeability tensor
			FEBiphasicMaterialPoint& pt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
			pt.m_p = 0;
			pt.m_w = vec3d(0,0,0);
					
			double K[3][3];
			biph->Permeability(K, mp);

			eps = n.x*(K[0][0]*n.x+K[0][1]*n.y+K[0][2]*n.z)
			+n.y*(K[1][0]*n.x+K[1][1]*n.y+K[1][2]*n.z)
			+n.z*(K[2][0]*n.x+K[2][1]*n.y+K[2][2]*n.z);
		}
	}

	return eps;
}

//-----------------------------------------------------------------------------
void FESlidingInterface2::ProjectSurface(FESlidingSurface2& ss, FESlidingSurface2& ms, bool bupseg)
{
	FEMesh& mesh = GetFEModel()->GetMesh();
	FESurfaceElement* pme;
	vec3d r, nu;
	double rs[2];
	double Ln;

	double ps[FEElement::MAX_NODES], p1;

	double R = m_srad*mesh.GetBoundingBox().radius();

	FENormalProjection np(ms);
	np.SetTolerance(m_stol);
	np.SetSearchRadius(R);
	np.Init();

	// loop over all integration points
 //   #pragma omp parallel for shared(R, bupseg)
	for (int i=0; i<ss.Elements(); ++i)
	{
		FESurfaceElement& el = ss.Element(i);
		bool sporo = ss.m_poro[i];

		int ne = el.Nodes();
		int nint = el.GaussPoints();

		// get the nodal pressures
		if (sporo)
		{
			for (int j=0; j<ne; ++j) ps[j] = mesh.Node(el.m_node[j]).m_pt;
		}

		for (int j=0; j<nint; ++j)
		{
			// get the integration point data
			FESlidingSurface2::Data& pt = ss.m_Data[i][j];

			// calculate the global position of the integration point
			r = ss.Local2Global(el, j);

			// get the pressure at the integration point
			if (sporo) p1 = el.eval(ps, j);

			// calculate the normal at this integration point
			nu = ss.SurfaceNormal(el, j);

			// first see if the old intersected face is still good enough
			pme = pt.m_pme;
			if (pme)
			{
				double g;

				// see if the ray intersects this element
				if (ms.Intersect(*pme, r, nu, rs, g, m_stol))
				{
					pt.m_rs[0] = rs[0];
					pt.m_rs[1] = rs[1];
				}
				else
				{
					pme = 0;
				}
			}

			// find the intersection point with the master surface
			if (pme == 0 && bupseg) pme = np.Project(r, nu, rs);

			pt.m_pme = pme;
			pt.m_nu = nu;
			pt.m_rs[0] = rs[0];
			pt.m_rs[1] = rs[1];
			if (pme)
			{
				// the node could potentially be in contact
				// find the global location of the intersection point
				vec3d q = ms.Local2Global(*pme, rs[0], rs[1]);

				// calculate the gap function
				// NOTE: this has the opposite sign compared
				// to Gerard's notes.
				double g = nu*(r - q);

				double eps = m_epsn*pt.m_epsn;

				Ln = pt.m_Lmd + eps*g;

				pt.m_gap = (g <= R? g : 0);

				if ((Ln >= 0) && (g <= R))
				{

					// calculate the pressure gap function
					bool mporo = ms.m_poro[pme->m_lid];
					if (sporo && mporo) {
						double pm[FEElement::MAX_NODES];
						for (int k=0; k<pme->Nodes(); ++k) pm[k] = mesh.Node(pme->m_node[k]).m_pt;
						double p2 = pme->eval(pm, rs[0], rs[1]);
						pt.m_pg = p1 - p2;
					}
				}
				else
				{
//					pt.m_Lmd = 0;
					pt.m_gap = 0;
					pt.m_pme = 0;
					if (sporo) {
//						pt.m_Lmp = 0;
						pt.m_pg = 0;
					}
				}
			}
			else
			{
				// the node is not in contact
				pt.m_Lmd = 0;
				pt.m_gap = 0;
				if (sporo) {
					pt.m_Lmp = 0;
					pt.m_pg = 0;
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------

void FESlidingInterface2::Update(int niter)
{	
	int i, n, id, np;
	double rs[2];

	double R = m_srad*GetFEModel()->GetMesh().GetBoundingBox().radius();

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
	niter = pstep->m_psolver->m_niter - biter;
	bool bupseg = ((m_nsegup == 0)? true : (niter <= m_nsegup));
	// get the logfile
//	Logfile& log = GetLogfile();
//	log.printf("seg_up iteration # %d\n", niter+1);
	
	// project the surfaces onto each other
	// this will update the gap functions as well
	ProjectSurface(m_ss, m_ms, bupseg);
	if (m_btwo_pass || m_ms.m_bporo) ProjectSurface(m_ms, m_ss, bupseg);

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

		// loop over all elements of the primary surface
		for (n=0; n<ss.Elements(); ++n)
		{
			FESurfaceElement& el = ss.Element(n);
			int nint = el.GaussPoints();
			int neln = el.Nodes();

			// get the normal tractions at the integration points
			double ti[FEElement::MAX_NODES], gap, eps;
			for (i=0; i<nint; ++i) 
			{
				FESlidingSurface2::Data& pt = ss.m_Data[n][i];
				gap = pt.m_gap;
				eps = m_epsn*pt.m_epsn;
				ti[i] = MBRACKET(pt.m_Lmd + eps*gap);
			}

			// project the data to the nodes
			double tn[FEElement::MAX_NODES];
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
		if (ms.m_bporo) {
			FENormalProjection np(ss);
			np.SetTolerance(m_stol);
			np.SetSearchRadius(R);
			np.Init();

			for (n=0; n<ms.Nodes(); ++n)
			{
				// get the node
				FENode& node = ms.Node(n);
				
				// project it onto the primary surface
				FESurfaceElement* pse = np.Project(node.m_rt, ms.m_nn[n], rs);
				
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
						double ti[FEElement::MAX_NODES], gap, eps;
						int nint = pse->GaussPoints();
						for (i=0; i<nint; ++i) 
						{
							FESlidingSurface2::Data& pt = ss.m_Data[pse->m_lid][i];
							gap = pt.m_gap;
							eps = m_epsn*pt.m_epsn;
							ti[i] = MBRACKET(pt.m_Lmd + eps*gap);
						}
						
						// project the data to the nodes
						double tn[FEElement::MAX_NODES];
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
}

//-----------------------------------------------------------------------------
void FESlidingInterface2::ShallowCopy(DumpStream& dmp, bool bsave)
{
	m_ss.ShallowCopy(dmp, bsave);
	m_ms.ShallowCopy(dmp, bsave);
}

//-----------------------------------------------------------------------------
void FESlidingInterface2::ContactForces(FEGlobalVector& R)
{
	int i, j, k;
	vector<int> sLM, mLM, LM, en;
	vector<double> fe;
	const int MN = FEElement::MAX_NODES;
	double detJ[MN], w[MN], *Hs, Hm[MN];
	double N[4*MN*2]; // TODO: is the size correct?

	FEModel& fem = *GetFEModel();

	// get the mesh
	FEMesh* pm = m_ss.GetMesh();

	// if we're using the symmetric formulation
	// we need to multiply with the timestep
	double dt = fem.GetCurrentStep()->m_dt;

	// loop over the nr of passes
	int npass = (m_btwo_pass?2:1);
	for (int np=0; np<npass; ++np)
	{
		// get slave and master surface
		FESlidingSurface2& ss = (np == 0? m_ss : m_ms);
		FESlidingSurface2& ms = (np == 0? m_ms : m_ss);

		// loop over all slave elements
		for (i=0; i<ss.Elements(); ++i)
		{
			// get the surface element
			FESurfaceElement& se = ss.Element(i);

			bool sporo = ss.m_poro[i];

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
				// get the integration point data
				FESlidingSurface2::Data& pt = ss.m_Data[i][j];

				// get the master element
				FESurfaceElement* pme = pt.m_pme;
				if (pme)
				{
					// get the master element
					FESurfaceElement& me = *pme;

					bool mporo = ms.m_poro[pme->m_lid];

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
					double r = pt.m_rs[0];
					double s = pt.m_rs[1];
					me.shape_fnc(Hm, r, s);

					// get normal vector
					vec3d nu = pt.m_nu;

					// gap function
					double g = pt.m_gap;
					
					// lagrange multiplier
					double Lm = pt.m_Lmd;

					// penalty 
					double eps = m_epsn*pt.m_epsn;

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
					R.Assemble(en, LM, fe);

					// do the biphasic stuff
					// TODO: I should only do this when the node is actually in contact
					//       in other words, when g >= 0
					if (sporo && mporo && (tn > 0))
					{
						// calculate nr of pressure dofs
						int ndof = nseln + nmeln;

						// calculate the flow rate
						double epsp = m_epsp*pt.m_epsp;

						double wn = pt.m_Lmp + epsp*pt.m_pg;

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
						R.Assemble(en, LM, fe);
					}
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------
void FESlidingInterface2::ContactStiffness(FESolver* psolver)
{
	int i, j, k, l;
	vector<int> sLM, mLM, LM, en;
	const int MN = FEElement::MAX_NODES;
	double detJ[MN], w[MN], *Hs, Hm[MN], pt[MN], dpr[MN], dps[MN];
	double N[4*MN*2];
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
		FESlidingSurface2& ss = (np == 0? m_ss : m_ms);
		FESlidingSurface2& ms = (np == 0? m_ms : m_ss);

		// loop over all slave elements
		for (i=0; i<ss.Elements(); ++i)
		{
			// get ths slave element
			FESurfaceElement& se = ss.Element(i);

			bool sporo = ss.m_poro[i];

			// get nr of nodes and integration points
			int nseln = se.Nodes();
			int nint = se.GaussPoints();

			// nodal pressures
			double pn[MN];
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
			for (j=0; j<nint; ++j)
			{
				// get integration point data
				FESlidingSurface2::Data& pt = ss.m_Data[i][j];

				// get the master element
				FESurfaceElement* pme = pt.m_pme;
				if (pme)
				{
					FESurfaceElement& me = *pme;

					bool mporo = ms.m_poro[pme->m_lid];

					// get the nr of master nodes
					int nmeln = me.Nodes();

					// nodal pressure
					double pm[MN];
					for (k=0; k<nmeln; ++k) pm[k] = ms.GetMesh()->Node(me.m_node[k]).m_pt;

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
					double r = pt.m_rs[0];
					double s = pt.m_rs[1];
					me.shape_fnc(Hm, r, s);

					// get slave normal vector
					vec3d nu = pt.m_nu;

					// gap function
					double g = pt.m_gap;
					
					// lagrange multiplier
					double Lm = pt.m_Lmd;

					// penalty 
					double eps = m_epsn*pt.m_epsn;

					// contact traction
					double tn = Lm + eps*g;
					tn = MBRACKET(tn);

//					double dtn = m_eps*HEAVYSIDE(Lm + eps*g);
					double dtn = (tn > 0.? eps :0.);
					
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
					mat3d As[FEElement::MAX_NODES];
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
					
					double Hmr[FEElement::MAX_NODES], Hms[FEElement::MAX_NODES];
					me.shape_deriv(Hmr, Hms, r, s);
					vec3d mm[FEElement::MAX_NODES];
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
						double dt = fem.GetCurrentStep()->m_dt;
						
						double epsp = (tn > 0) ? m_epsp*pt.m_epsp : 0.;
						
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
							
							double wn = pt.m_Lmp + epsp*pt.m_pg;
							
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
	int npass = (m_btwo_pass?2:1);
	for (int np=0; np<npass; ++np)
	{
		FESlidingSurface2& ss = (np == 0? m_ss : m_ms);
		FESlidingSurface2& ms = (np == 0? m_ms : m_ss);
		
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
				FESlidingSurface2::Data& pt = ss.m_Data[n][i];

				gap = pt.m_gap;
				eps = m_epsn*pt.m_epsn;
				pt.m_Ln = MBRACKET(pt.m_Lmd + eps*gap);
				FESurfaceElement* pme = pt.m_pme;
				if (m_btwo_pass && pme)
				{
					int mint = pme->GaussPoints();
					double ti[FEElement::MAX_NODES];
					for (int j=0; j<mint; ++j) {
						FESlidingSurface2::Data& mp = ms.m_Data[pme->m_lid][j];
						gap = mp.m_gap;
						eps = m_epsn*mp.m_epsn;
						ti[j] = MBRACKET(mp.m_Lmd + m_epsn*mp.m_epsn*mp.m_gap);
					}
					// project the data to the nodes
					double tn[FEElement::MAX_NODES];
					pme->project_to_nodes(ti, tn);
					// now evaluate the traction at the intersection point
					double Ln = pme->eval(tn, pt.m_rs[0], pt.m_rs[1]);
					pt.m_Ln += MBRACKET(Ln);
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
	int NS = m_ss.m_Data.size();
	int NM = m_ms.m_Data.size();

	// --- c a l c u l a t e   i n i t i a l   n o r m s ---
	// a. normal component
	double normL0 = 0, normP = 0, normDP = 0;
	for (int i=0; i<NS; ++i)
	{
		vector<FESlidingSurface2::Data>& sd = m_ss.m_Data[i];
		for (int j=0; j<(int)sd.size(); ++j)
		{
			FESlidingSurface2::Data& ds = sd[j];
			normL0 += ds.m_Lmd*ds.m_Lmd;
		}
	}
	for (int i=0; i<NM; ++i)
	{
		vector<FESlidingSurface2::Data>& md = m_ms.m_Data[i];
		for (int j=0; j<(int)md.size(); ++j)
		{
			FESlidingSurface2::Data& dm = md[j];
			normL0 += dm.m_Lmd*dm.m_Lmd;
		}
	}

	// b. gap component
	// (is calculated during update)
	double maxgap = 0;
	double maxpg = 0;

	// update Lagrange multipliers
	double normL1 = 0, eps, epsp;
	for (i=0; i<NS; ++i)
	{
		vector<FESlidingSurface2::Data>& sd = m_ss.m_Data[i];
		for (int j=0; j<(int)sd.size(); ++j)
		{
			FESlidingSurface2::Data& ds = sd[j];

			// update Lagrange multipliers on slave surface
			eps = m_epsn*ds.m_epsn;
			Ln = ds.m_Lmd + eps*ds.m_gap;
			ds.m_Lmd = MBRACKET(Ln);

			normL1 += ds.m_Lmd*ds.m_Lmd;
		
			if (m_ss.m_bporo) {
				Lp = 0;
				if (Ln > 0) {
					epsp = m_epsp*ds.m_epsp;
					Lp = ds.m_Lmp + epsp*ds.m_pg;
					maxpg = max(maxpg,fabs(ds.m_pg));
					normDP += ds.m_pg*ds.m_pg;
				}
				ds.m_Lmp = Lp;
			}
		
			if (Ln > 0) maxgap = max(maxgap,fabs(ds.m_gap));
		}
	}	
	
	for (i=0; i<NM; ++i)
	{
		vector<FESlidingSurface2::Data>& md = m_ms.m_Data[i];
		for (int j=0; j<(int)md.size(); ++j)
		{
			FESlidingSurface2::Data& dm = md[j];

			// update Lagrange multipliers on master surface
			eps = m_epsn*dm.m_epsn;
			Ln = dm.m_Lmd + eps*dm.m_gap;
			dm.m_Lmd = MBRACKET(Ln);

			normL1 += dm.m_Lmd*dm.m_Lmd;
		
			if (m_ms.m_bporo) {
				Lp = 0;
				if (Ln > 0) {
					epsp = m_epsp*dm.m_epsp;
					Lp = dm.m_Lmp + epsp*dm.m_pg;
					maxpg = max(maxpg,fabs(dm.m_pg));
					normDP += dm.m_pg*dm.m_pg;
				}
				dm.m_Lmp = Lp;
			}
		
			if (Ln > 0) maxgap = max(maxgap,fabs(dm.m_gap));
		}
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

	if (naug < m_naugmin ) bconv = false;
	if (naug >= m_naugmax) bconv = true;

	felog.printf(" sliding interface # %d\n", m_nID);
	felog.printf("                        CURRENT        REQUIRED\n");
	felog.printf("    D multiplier : %15le", lnorm); if (m_atol > 0) felog.printf("%15le\n", m_atol); else felog.printf("       ***\n");
	if (bporo) { felog.printf("    P gap        : %15le", pnorm); if (m_atol > 0) felog.printf("%15le\n", m_atol); else felog.printf("       ***\n"); }

	felog.printf("    maximum gap  : %15le", maxgap);
	if (m_gtol > 0) felog.printf("%15le\n", m_gtol); else felog.printf("       ***\n");
	if (bporo) {
		felog.printf("    maximum pgap : %15le", maxpg);
		if (m_ptol > 0) felog.printf("%15le\n", m_ptol); else felog.printf("       ***\n");
	}
	
	return bconv;
}

//-----------------------------------------------------------------------------
void FESlidingInterface2::Serialize(DumpFile &ar)
{
	// serialize contact data
	FEContactInterface::Serialize(ar);

	// serialize contact surface data
	m_ms.Serialize(ar);
	m_ss.Serialize(ar);
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
