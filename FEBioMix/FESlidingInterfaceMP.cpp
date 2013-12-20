#include "stdafx.h"
#include "FESlidingInterfaceMP.h"
#include "FEBiphasic.h"
#include "FEBiphasicSolute.h"
#include "FEMultiphasic.h"
#include "FEBioMech/FEStiffnessMatrix.h"
#include "FECore/FEModel.h"
#include "FECore/log.h"

//-----------------------------------------------------------------------------
// Define sliding interface parameters
BEGIN_PARAMETER_LIST(FESlidingInterfaceMP, FEContactInterface)
ADD_PARAMETER(m_blaugon  , FE_PARAM_BOOL  , "laugon"               ); 
ADD_PARAMETER(m_atol     , FE_PARAM_DOUBLE, "tolerance"            );
ADD_PARAMETER(m_gtol     , FE_PARAM_DOUBLE, "gaptol"               );
ADD_PARAMETER(m_ptol     , FE_PARAM_DOUBLE, "ptol"                 );
ADD_PARAMETER(m_ctol     , FE_PARAM_DOUBLE, "ctol"                 );
ADD_PARAMETER(m_epsn     , FE_PARAM_DOUBLE, "penalty"              );
ADD_PARAMETER(m_bautopen , FE_PARAM_BOOL  , "auto_penalty"         );
ADD_PARAMETER(m_btwo_pass, FE_PARAM_BOOL  , "two_pass"             );
ADD_PARAMETER(m_knmult   , FE_PARAM_INT   , "knmult"               );
ADD_PARAMETER(m_stol     , FE_PARAM_DOUBLE, "search_tol"           );
ADD_PARAMETER(m_epsp     , FE_PARAM_DOUBLE, "pressure_penalty"     );
ADD_PARAMETER(m_epsc     , FE_PARAM_DOUBLE, "concentration_penalty");
ADD_PARAMETER(m_bsymm    , FE_PARAM_BOOL  , "symmetric_stiffness"  );
ADD_PARAMETER(m_srad     , FE_PARAM_DOUBLE, "search_radius"        );
ADD_PARAMETER(m_nsegup   , FE_PARAM_INT   , "seg_up"               );
ADD_PARAMETER(m_naugmin  , FE_PARAM_INT   , "minaug"             );
ADD_PARAMETER(m_naugmax  , FE_PARAM_INT   , "maxaug"             );
ADD_PARAMETER(m_ambp     , FE_PARAM_DOUBLE, "ambient_pressure"     );
ADD_PARAMETER(m_ambc[0]  , FE_PARAM_DOUBLE, "ambient_concentration_1");
ADD_PARAMETER(m_ambc[1]  , FE_PARAM_DOUBLE, "ambient_concentration_2");
ADD_PARAMETER(m_ambc[2]  , FE_PARAM_DOUBLE, "ambient_concentration_3");
ADD_PARAMETER(m_ambc[3]  , FE_PARAM_DOUBLE, "ambient_concentration_4");
ADD_PARAMETER(m_ambc[4]  , FE_PARAM_DOUBLE, "ambient_concentration_5");
ADD_PARAMETER(m_ambc[5]  , FE_PARAM_DOUBLE, "ambient_concentration_6");
ADD_PARAMETER(m_ambc[6]  , FE_PARAM_DOUBLE, "ambient_concentration_7");
ADD_PARAMETER(m_ambc[7]  , FE_PARAM_DOUBLE, "ambient_concentration_8");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FESlidingSurfaceMP::Data::Data()
{
	m_gap = 0.0;
	m_nu = vec3d(0,0,0);
	m_rs = vec2d(0,0);
	m_Lmd = 0.0;
	m_Lmp = 0.0;
	m_epsn = 1.0;
	m_epsp = 1.0;
	m_Ln = 0.0;
	m_pg = 0.0;
    
	m_pme = (FESurfaceElement*)(0);
}

//-----------------------------------------------------------------------------
// FESlidingSurfaceMP
//-----------------------------------------------------------------------------

FESlidingSurfaceMP::FESlidingSurfaceMP(FEModel* pfem) : FEContactSurface(&pfem->GetMesh())
{ 
	m_bporo = m_bsolu = false;
	m_pfem = pfem; 
}

//-----------------------------------------------------------------------------
bool FESlidingSurfaceMP::Init()
{
	// initialize surface data first
	if (FEContactSurface::Init() == false) return false;
	
	m_nn.assign(Nodes(), 0);
	
	// determine solutes for this surface using the first surface element
	// TODO: Check that all elements use the same set of solutes as the first element
	int nsol = 0;
	if (Elements()) {
		FESurfaceElement& se = Element(0);
		// get the solid element this surface element belongs to
		FESolidElement* pe = dynamic_cast<FESolidElement*>(m_pMesh->FindElementFromID(se.m_nelem));
		if (pe)
		{
			// get the material
			FEMaterial* pm = dynamic_cast<FEMaterial*>(m_pfem->GetMaterial(pe->GetMatID()));
			
			// check type of element
			FEBiphasic* pb = dynamic_cast<FEBiphasic*> (pm);
			FEBiphasicSolute* pbs = dynamic_cast<FEBiphasicSolute*> (pm);
			FEMultiphasic* pmp = dynamic_cast<FEMultiphasic*> (pm);
			if (pb) {
				m_bporo = true;
				nsol = 0;
			} else if (pbs) {
				m_bporo = m_bsolu = true;
				nsol = 1;
				m_sid.assign(nsol, pbs->GetSolute()->GetSoluteID());
			} else if (pmp) {
				m_bporo = m_bsolu = true;
				nsol = pmp->Solutes();
				m_sid.resize(nsol);
				for (int isol=0; isol<nsol; ++isol) {
					m_sid[isol] = pmp->GetSolute(isol)->GetSoluteID();
				}
			}
		}
	}
	
	// allocate data structures
	int NE = Elements();
	m_Data.resize(NE);
	for (int i=0; i<NE; ++i)
	{
		FESurfaceElement& el = Element(i);
		int nint = el.GaussPoints();
		m_Data[i].resize(nint);
        if (nsol) {
            for (int j=0; j<nint; ++j) {
                m_Data[i][j].m_Lmc.resize(nsol);
                m_Data[i][j].m_epsc.resize(nsol);
                m_Data[i][j].m_cg.resize(nsol);
            }
        }
	}
	
	return true;
}

//-----------------------------------------------------------------------------
void FESlidingSurfaceMP::ShallowCopy(DumpStream& dmp, bool bsave)
{
	// We need to store the m_bporo and m_bsolu flags first
	// since we need them before we initialize the surface data
	if (bsave)
	{
		dmp << m_bporo;
		dmp << m_bsolu;
	}
	else
	{
		dmp >> m_bporo;
		dmp >> m_bsolu;
	}
	
	// And finally, we serialize the surface data
	if (bsave)
	{
		int N = (int) m_Data.size();
		for (int i=0; i<N; ++i)
		{
			vector<Data>& di = m_Data[i];
			for (int j=0; j<(int) di.size(); ++j)
			{
				Data& d = di[j];
				dmp << d.m_gap;
				dmp << d.m_nu;
				dmp << d.m_rs;
				dmp << d.m_Lmd;
				dmp << d.m_Lmp;
				dmp << d.m_Lmc;
				dmp << d.m_epsn;
				dmp << d.m_epsp;
				dmp << d.m_epsc;
				dmp << d.m_pg;
				dmp << d.m_cg;
				dmp << d.m_Ln;
			}
		}
		dmp << m_nn;
        dmp << m_sid;
	}
	else
	{
		int N = (int) m_Data.size();
		for (int i=0; i<N; ++i)
		{
			vector<Data>& di = m_Data[i];
			for (int j=0; j<(int) di.size(); ++j)
			{
				Data& d = di[j];
				dmp >> d.m_gap;
				dmp >> d.m_nu;
				dmp >> d.m_rs;
				dmp >> d.m_Lmd;
				dmp >> d.m_Lmp;
				dmp >> d.m_Lmc;
				dmp >> d.m_epsn;
				dmp >> d.m_epsp;
				dmp >> d.m_epsc;
				dmp >> d.m_pg;
				dmp >> d.m_cg;
				dmp >> d.m_Ln;
			}
		}
		dmp >> m_nn;
        dmp >> m_sid;
        
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
//! This function calculates the node normal. Due to the piecewise continuity
//! of the surface elements this normal is not uniquely defined so in order to
//! obtain a unique normal the normal is averaged for each node over all the 
//! element normals at the node

void FESlidingSurfaceMP::UpdateNodeNormals()
{
	int N = Nodes(), i, j, ne, jp1, jm1;
	const int MN = FEElement::MAX_NODES;
	vec3d y[MN], n;
	
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
vec3d FESlidingSurfaceMP::GetContactForce()
{
	int n, i;
	
	// initialize contact force
	vec3d f(0,0,0);
	
	// loop over all elements of the primary surface
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
vec3d FESlidingSurfaceMP::GetFluidForce()
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
void FESlidingSurfaceMP::Serialize(DumpFile& ar)
{
	// We need to store the m_bporo and m_bsolu flags first
	// since we need them before we initialize the surface data
	if (ar.IsSaving())
	{
		ar << m_bporo;
		ar << m_bsolu;
	}
	else
	{
		ar >> m_bporo;
		ar >> m_bsolu;
	}
	
	// Next, we can serialize the base-class data
	FEContactSurface::Serialize(ar);
	
	// And finally, we serialize the surface data
	if (ar.IsSaving())
	{
		int N = (int) m_Data.size();
		for (int i=0; i<N; ++i)
		{
			vector<Data>& di = m_Data[i];
			for (int j=0; j<(int) di.size(); ++j)
			{
				Data& d = di[j];
				ar << d.m_gap;
				ar << d.m_nu;
				ar << d.m_rs;
				ar << d.m_Lmd;
				ar << d.m_Lmp;
				ar << d.m_Lmc;
				ar << d.m_epsn;
				ar << d.m_epsp;
				ar << d.m_epsc;
				ar << d.m_pg;
				ar << d.m_cg;
				ar << d.m_Ln;
			}
		}
		ar << m_nn;
        ar << m_sid;
	}
	else
	{
		int N = (int) m_Data.size();
		for (int i=0; i<N; ++i)
		{
			vector<Data>& di = m_Data[i];
			for (int j=0; j<(int) di.size(); ++j)
			{
				Data& d = di[j];
				ar >> d.m_gap;
				ar >> d.m_nu;
				ar >> d.m_rs;
				ar >> d.m_Lmd;
				ar >> d.m_Lmp;
				ar >> d.m_Lmc;
				ar >> d.m_epsn;
				ar >> d.m_epsp;
				ar >> d.m_epsc;
				ar >> d.m_pg;
				ar >> d.m_cg;
				ar >> d.m_Ln;
			}
		}
		ar >> m_nn;
        ar >> m_sid;
	}
}

//-----------------------------------------------------------------------------
void FESlidingSurfaceMP::GetNodalContactGap(int nface, double* gn)
{
	double gi[FEElement::MAX_INTPOINTS];
	FESurfaceElement& el = Element(nface);
	int ni = el.GaussPoints();
	for (int k=0; k<ni; ++k) gi[k] = m_Data[nface][k].m_gap;
	el.project_to_nodes(gi, gn);
}

//-----------------------------------------------------------------------------
void FESlidingSurfaceMP::GetNodalContactPressure(int nface, double* pn)
{
	FESurfaceElement& el = Element(nface);
	int ni = el.GaussPoints();
	double ti[FEElement::MAX_INTPOINTS];
	for (int k=0; k<ni; ++k) ti[k] = m_Data[nface][k].m_Ln;
	el.project_to_nodes(ti, pn);
}

//-----------------------------------------------------------------------------
void FESlidingSurfaceMP::GetNodalContactTraction(int nface, vec3d* tn)
{
	FESurfaceElement& el = Element(nface);
	int ne = el.Nodes();
	int ni = el.GaussPoints();
    
	vec3d t;
	const int MFI = FEElement::MAX_INTPOINTS;
	double tix[MFI], tiy[MFI], tiz[MFI];
	for (int k=0; k<ni; ++k)
	{
		FESlidingSurfaceMP::Data& pt = m_Data[nface][k];
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
		tn[k].x = tnx[k];
		tn[k].y = tny[k];
		tn[k].z = tnz[k];
	}
}

//-----------------------------------------------------------------------------
// FESlidingInterfaceMP
//-----------------------------------------------------------------------------

FESlidingInterfaceMP::FESlidingInterfaceMP(FEModel* pfem) : FEContactInterface(pfem), m_ss(pfem), m_ms(pfem)
{
	static int count = 1;
	m_nID = count++;
	
	// initial values
	m_knmult = 1;
	m_atol = 0.1;
	m_epsn = 1;
	m_epsp = 1;
	m_epsc = 1;
	m_btwo_pass = false;
	m_stol = 0.01;
	m_bsymm = true;
	m_srad = 1.0;
	m_gtol = 0;
	m_ptol = 0;
	m_ctol = 0;
	m_ambp = 0;
	for (int isol=0; isol<MAX_CDOFS; ++isol) m_ambc[isol] = 0;
	m_nsegup = 0;
	m_bautopen = false;
	
	m_naugmin = 0;
	m_naugmax = 10;
	
	m_ss.SetSibling(&m_ms);
	m_ms.SetSibling(&m_ss);
}

//-----------------------------------------------------------------------------

FESlidingInterfaceMP::~FESlidingInterfaceMP()
{
}

//-----------------------------------------------------------------------------
bool FESlidingInterfaceMP::Init()
{
	m_Rgas = GetFEModel()->GetGlobalConstant("R");
	m_Tabs = GetFEModel()->GetGlobalConstant("T");
	
	// initialize surface data
	if (m_ss.Init() == false) return false;
	if (m_ms.Init() == false) return false;
	
	// determine which solutes are common to both contact surfaces
	m_sid.clear(); m_ssl.clear(); m_msl.clear();
	for (int is=0; is<m_ss.m_sid.size(); ++is) {
		for (int im=0; im<m_ms.m_sid.size(); ++im) {
			if (m_ms.m_sid[im] == m_ss.m_sid[is]) {
				m_sid.push_back(m_ss.m_sid[is]);
				m_ssl.push_back(is);
				m_msl.push_back(im);
			}
		}
	}
	
	return true;
}

//-----------------------------------------------------------------------------
//! build the matrix profile for use in the stiffness matrix
void FESlidingInterfaceMP::BuildMatrixProfile(FEStiffnessMatrix& K)
{
	FEMesh& mesh = GetFEModel()->GetMesh();
    
	vector<int> lm(8*FEElement::MAX_NODES*2);
    
	int npass = (m_btwo_pass?2:1);
	for (int np=0; np<npass; ++np)
	{
		FESlidingSurfaceMP& ss = (np == 0? m_ss : m_ms);
		FESlidingSurfaceMP& ms = (np == 0? m_ms : m_ss);
        
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
						lm[8*l  ] = id[DOF_X];
						lm[8*l+1] = id[DOF_Y];
						lm[8*l+2] = id[DOF_Z];
						lm[8*l+3] = id[DOF_P];
						lm[8*l+4] = id[DOF_RU];
						lm[8*l+5] = id[DOF_RV];
						lm[8*l+6] = id[DOF_RW];
                        for (int m=0; m<ss.m_sid.size(); ++m) {
                            lm[8*l+7+m] = id[DOF_C + ss.m_sid[m]];
                        }
					}
                    
					for (l=0; l<nmeln; ++l)
					{
						int* id = mesh.Node(mn[l]).m_ID;
						lm[8*(l+nseln)  ] = id[DOF_X];
						lm[8*(l+nseln)+1] = id[DOF_Y];
						lm[8*(l+nseln)+2] = id[DOF_Z];
						lm[8*(l+nseln)+3] = id[DOF_P];
						lm[8*(l+nseln)+4] = id[DOF_RU];
						lm[8*(l+nseln)+5] = id[DOF_RV];
						lm[8*(l+nseln)+6] = id[DOF_RW];
                        for (int m=0; m<ms.m_sid.size(); ++m) {
                            lm[8*(l+nseln)+7+m] = id[DOF_C + ms.m_sid[m]];
                        }
					}
                    
					K.build_add(lm);
				}
			}
		}
	}
    
}

//-----------------------------------------------------------------------------
void FESlidingInterfaceMP::Activate()
{
	// don't forget to call the base members
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
		if (m_ss.m_bporo) CalcAutoPressurePenalty(m_ss);
		for (int is=0; is<m_ssl.size(); ++is)
			CalcAutoConcentrationPenalty(m_ss, m_ssl[is]);
		if (m_ms.m_bporo) CalcAutoPressurePenalty(m_ms);
		for (int im=0; im<m_msl.size(); ++im)
			CalcAutoConcentrationPenalty(m_ms, m_msl[im]);
	}
	
	// update sliding interface data
	Update(0);
}

//-----------------------------------------------------------------------------
void FESlidingInterfaceMP::CalcAutoPenalty(FESlidingSurfaceMP& s)
{
	// get the mesh
	FEMesh& m = GetFEModel()->GetMesh();
	
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
		double E = AutoPenalty(el, s);
		
		// calculate penalty
		double eps = E*A/V;
		
		// assign to integation points of surface element
		int nint = el.GaussPoints();
		for (int j=0; j<nint; ++j, ++ni)
        {
			FESlidingSurfaceMP::Data& pt = s.m_Data[i][j];
            pt.m_epsn = eps;
        }
	}
}

//-----------------------------------------------------------------------------
void FESlidingInterfaceMP::CalcAutoPressurePenalty(FESlidingSurfaceMP& s)
{
	// get the mesh
	FEMesh& m = GetFEModel()->GetMesh();
	
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
		for (int j=0; j<nint; ++j, ++ni)
        {
			FESlidingSurfaceMP::Data& pt = s.m_Data[i][j];
            pt.m_epsp = eps;
        }
	}
}

//-----------------------------------------------------------------------------

double FESlidingInterfaceMP::AutoPressurePenalty(FESurfaceElement& el, FESlidingSurfaceMP& s)
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
		FEMaterial* pm = dynamic_cast<FEMaterial*>(GetFEModel()->GetMaterial(pe->GetMatID()));
		
		// check type of element
		FEBiphasic* pb = dynamic_cast<FEBiphasic*> (pm);
		FEBiphasicSolute* pbs = dynamic_cast<FEBiphasicSolute*> (pm);
		FEMultiphasic* pmp = dynamic_cast<FEMultiphasic*> (pm);
		if (pb)
		{
			// get a material point
			FEMaterialPoint& mp = *pe->m_State[0];
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
			pb->Permeability(K, mp);
			
			eps = n.x*(K[0][0]*n.x+K[0][1]*n.y+K[0][2]*n.z)
			+n.y*(K[1][0]*n.x+K[1][1]*n.y+K[1][2]*n.z)
			+n.z*(K[2][0]*n.x+K[2][1]*n.y+K[2][2]*n.z);
		}
		else if (pbs)
		{
			// get a material point
			FEMaterialPoint& mp = *pe->m_State[0];
			FEElasticMaterialPoint& ept = *(mp.ExtractData<FEElasticMaterialPoint>());
			
			// setup the material point
			ept.m_F = mat3dd(1.0);
			ept.m_J = 1;
			ept.m_s.zero();
			
			// if this is a biphasic-solute element, then get the permeability tensor
			FEBiphasicMaterialPoint& ppt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
			FESoluteMaterialPoint& spt = *(mp.ExtractData<FESoluteMaterialPoint>());
			ppt.m_p = 0;
			ppt.m_w = vec3d(0,0,0);
			spt.m_c = 0;
			spt.m_j = vec3d(0,0,0);
			
			mat3ds K = pbs->GetPermeability()->Permeability(mp);
			
			eps = n*(K*n);
		}
		else if (pmp)
		{
			// get a material point
			FEMaterialPoint& mp = *pe->m_State[0];
			FEElasticMaterialPoint& ept = *(mp.ExtractData<FEElasticMaterialPoint>());
			
			// setup the material point
			ept.m_F = mat3dd(1.0);
			ept.m_J = 1;
			ept.m_s.zero();
			
			// get the permeability tensor
			FEBiphasicMaterialPoint& ppt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
			FESolutesMaterialPoint& spt = *(mp.ExtractData<FESolutesMaterialPoint>());
			ppt.m_p = 0;
			ppt.m_w = vec3d(0,0,0);
			for (int i=0; i<pmp->Solutes(); ++i) {
				spt.m_c[i] = 0;
				spt.m_j[i] = vec3d(0,0,0);
			}
			
			mat3ds K = pmp->GetPermeability()->Permeability(mp);
			
			eps = n*(K*n);
		}
	}
	
	return eps;
}

//-----------------------------------------------------------------------------
void FESlidingInterfaceMP::CalcAutoConcentrationPenalty(FESlidingSurfaceMP& s,
														const int isol)
{
	// get the mesh
	FEMesh& m = GetFEModel()->GetMesh();
	
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
		double d = AutoConcentrationPenalty(el, s, isol);
		
		// calculate penalty
		double eps = d*A/V;
		
		// assign to integation points of surface element
		int nint = el.GaussPoints();
		for (int j=0; j<nint; ++j, ++ni)
        {
			FESlidingSurfaceMP::Data& pt = s.m_Data[i][j];
            pt.m_epsc[isol] = eps;
        }
	}
}

//-----------------------------------------------------------------------------

double FESlidingInterfaceMP::AutoConcentrationPenalty(FESurfaceElement& el, 
													  FESlidingSurfaceMP& s,
													  const int isol)
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
		FEMaterial* pm = dynamic_cast<FEMaterial*>(GetFEModel()->GetMaterial(pe->GetMatID()));
		
		// see if this is a biphasic-solute or multiphasic element
		FEBiphasicSolute* pbs = dynamic_cast<FEBiphasicSolute*> (pm);
		FEMultiphasic* pmp = dynamic_cast<FEMultiphasic*> (pm);
		if (pbs)
		{
			// get a material point
			FEMaterialPoint& mp = *pe->m_State[0];
			FEElasticMaterialPoint& ept = *(mp.ExtractData<FEElasticMaterialPoint>());
			
			// setup the material point
			ept.m_F = mat3dd(1.0);
			ept.m_J = 1;
			ept.m_s.zero();
			
			// for a biphasic-solute element, get the diffusivity tensor
			FEBiphasicMaterialPoint& ppt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
			FESoluteMaterialPoint& spt = *(mp.ExtractData<FESoluteMaterialPoint>());
			ppt.m_p = 0;
			ppt.m_w = vec3d(0,0,0);
			spt.m_c = 0;
			spt.m_j = vec3d(0,0,0);
			
			mat3ds D = pbs->GetSolute()->m_pDiff->Diffusivity(mp)
			*(pbs->Porosity(mp)*pbs->GetSolute()->m_pSolub->Solubility(mp));
			
			// evaluate normal component of diffusivity
			eps = n*(D*n);
		}
		else if (pmp)
		{
			// get a material point
			FEMaterialPoint& mp = *pe->m_State[0];
			FEElasticMaterialPoint& ept = *(mp.ExtractData<FEElasticMaterialPoint>());
			
			// setup the material point
			ept.m_F = mat3dd(1.0);
			ept.m_J = 1;
			ept.m_s.zero();
			
			// for a multiphasic element, get the diffusivity tensor for that solute
			FEBiphasicMaterialPoint& ppt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
			FESolutesMaterialPoint& spt = *(mp.ExtractData<FESolutesMaterialPoint>());
			ppt.m_p = 0;
			ppt.m_w = vec3d(0,0,0);
			for (int i=0; i<pmp->Solutes(); ++i) {
				spt.m_c[i] = 0;
				spt.m_j[i] = vec3d(0,0,0);
			}
			
			mat3ds D = pmp->GetSolute(isol)->m_pDiff->Diffusivity(mp)
			*(pmp->Porosity(mp)*pmp->GetSolute(isol)->m_pSolub->Solubility(mp));
			
			// evaluate normal component of diffusivity
			eps = n*(D*n);
		}
	}
	
	return eps;
}

//-----------------------------------------------------------------------------
void FESlidingInterfaceMP::ProjectSurface(FESlidingSurfaceMP& ss, FESlidingSurfaceMP& ms, bool bupseg)
{
	bool bfirst = true;
	
	FEMesh& mesh = GetFEModel()->GetMesh();
	FESurfaceElement* pme;
	vec3d r, nu;
	double rs[2];
	double Ln;
	
	const int MN = FEElement::MAX_NODES;
	int nsol = (int)m_sid.size();
	double ps[MN], p1;
	vector< vector<double> > cs(nsol, vector<double>(MN));
	vector<double> c1(nsol);
	
	double R = m_srad*mesh.GetBoundingBox().radius();
	
	// loop over all integration points
	for (int i=0; i<ss.Elements(); ++i)
	{
		FESurfaceElement& el = ss.Element(i);
		
		bool sporo = ss.m_bporo;
		
		int ne = el.Nodes();
		int nint = el.GaussPoints();
		
		// get the nodal pressures
		if (sporo)
		{
			for (int j=0; j<ne; ++j) ps[j] = mesh.Node(el.m_node[j]).m_pt;
		}
		
		// get the nodal concentrations
		for (int isol=0; isol<nsol; ++isol) {
			for (int j=0; j<ne; ++j) cs[isol][j] = mesh.Node(el.m_node[j]).m_ct[m_sid[isol]];
		}
		
		for (int j=0; j<nint; ++j)
		{
			FESlidingSurfaceMP::Data& pt = ss.m_Data[i][j];
            
			// calculate the global position of the integration point
			r = ss.Local2Global(el, j);
			
			// get the pressure at the integration point
			if (sporo) p1 = el.eval(ps, j);
			
			// get the concentration at the integration point
			for (int isol=0; isol<nsol; ++isol) c1[isol] = el.eval(&cs[isol][0], j);
			
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
			if (pme == 0 && bupseg) pme = ms.FindIntersection(r, nu, rs, bfirst, m_stol, R);
			
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
					bool mporo = ms.m_bporo;
					if (sporo && mporo) {
						double pm[MN];
						for (int k=0; k<pme->Nodes(); ++k) pm[k] = mesh.Node(pme->m_node[k]).m_pt;
						double p2 = pme->eval(pm, rs[0], rs[1]);
						pt.m_pg = p1 - p2;
					}
					for (int isol=0; isol<nsol; ++isol) {
						int sid = m_sid[isol];
						double cm[MN];
						for (int k=0; k<pme->Nodes(); ++k) cm[k] = mesh.Node(pme->m_node[k]).m_ct[sid];
						double c2 = pme->eval(cm, rs[0], rs[1]);
						pt.m_cg[m_ssl[isol]] = c1[isol] - c2;
					}
				}
				else
				{
					pt.m_gap = 0;
					pt.m_pme = 0;
					if (sporo) {
						pt.m_pg = 0;
					}
					for (int isol=0; isol<nsol; ++isol) {
						pt.m_cg[m_ssl[isol]] = 0;
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
				for (int isol=0; isol<nsol; ++isol) {
					pt.m_Lmc[m_ssl[isol]] = 0;
					pt.m_cg[m_ssl[isol]] = 0;
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------

void FESlidingInterfaceMP::Update(int niter)
{
	int i, j, n, id, np;
	double rs[2];
	const int MN = FEElement::MAX_NODES;
	
	FEModel& fem = *GetFEModel();
	
	double R = m_srad*GetFEModel()->GetMesh().GetBoundingBox().radius();
	
	static int naug = 0;
	static int biter = 0;
	
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
	if (m_btwo_pass || m_ss.m_bporo) ProjectSurface(m_ms, m_ss, bupseg);
	
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
	// if we need to modify the constraints on the pressure and concentration dofs.
	// If the nodes are not in contact, they must match ambient
	// conditions. Since all nodes have been previously marked to be
	// under ambient conditions in MarkAmbient(), we just need to reverse
	// this setting here, for nodes that are in contact.
	
	// Next, we loop over each surface, visiting the nodes
	// and finding out if that node is in contact or not
	int npass = (m_btwo_pass?2:1);
	for (np=0; np<npass; ++np)
	{
		FESlidingSurfaceMP& ss = (np == 0? m_ss : m_ms);
		FESlidingSurfaceMP& ms = (np == 0? m_ms : m_ss);
		
		// loop over all elements of the primary surface
		for (n=0; n<ss.Elements(); ++n)
		{
			FESurfaceElement& el = ss.Element(n);
			int nint = el.GaussPoints();
			int neln = el.Nodes();
			
			// get the normal tractions at the integration points
			double ti[MN], gap, eps;
			for (i=0; i<nint; ++i)
			{
				FESlidingSurfaceMP::Data& pt = ss.m_Data[n][i];
				gap = pt.m_gap;
				eps = m_epsn*pt.m_epsn;
				ti[i] = MBRACKET(pt.m_Lmd + eps*gap);
			}
			
			// project the data to the nodes
			double tn[MN];
			el.project_to_nodes(ti, tn);
			
			// see if we need to make changes to the BC's
			for (i=0; i<neln; ++i)
			{
				if (tn[i] > 0)
				{
					FENode& node = ss.Node(el.m_lnode[i]);
					id = node.m_ID[DOF_P];
					if (id < -1)
						// mark node as non-ambient (= pos ID)
						node.m_ID[DOF_P] = -id-2;
					for (j=0; j<MAX_CDOFS; ++j) {
						id = node.m_ID[DOF_C+j];
						if (id < -1)
							// mark node as non-ambient (= pos ID)
							node.m_ID[DOF_C+j] = -id-2;
					}
				}
			}
		}
		
		// loop over all nodes of the secondary surface
		// the secondary surface is trickier since we need
		// to look at the primary's surface projection
		if (ms.m_bporo) {
			for (n=0; n<ms.Nodes(); ++n)
			{
				// get the node
				FENode& node = ms.Node(n);
				
				// project it onto the primary surface
				bool bfirst = false;
				FESurfaceElement* pse = ss.FindIntersection(node.m_rt, ms.m_nn[n], rs, bfirst, m_stol, R);
				
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
						double ti[MN], gap, eps;
						int nint = pse->GaussPoints();
						vector<FESlidingSurfaceMP::Data>& sd = ss.m_Data[pse->m_lid];
						for (i=0; i<nint; ++i)
						{
							gap = sd[i].m_gap;
							eps = m_epsn*sd[i].m_epsn;
							ti[i] = MBRACKET(sd[i].m_Lmd + eps*gap);
						}
						
						// project the data to the nodes
						double tn[MN];
						pse->project_to_nodes(ti, tn);
						
						// now evaluate the traction at the intersection point
						double tp = pse->eval(tn, rs[0], rs[1]);
						
						// if tp > 0, mark node as non-ambient. (= pos ID)
						if (tp > 0)  {
							id = node.m_ID[DOF_P];
							if (id < -1)
								// mark as non-ambient
								node.m_ID[DOF_P] = -id-2;
							for (j=0; j<MAX_CDOFS; ++j) {
								id = node.m_ID[DOF_C+j];
								if (id < -1)
									// mark as non-ambient
									node.m_ID[DOF_C+j] = -id-2;
							}
						}
					}
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------
void FESlidingInterfaceMP::ShallowCopy(DumpStream& dmp, bool bsave)
{
	m_ss.ShallowCopy(dmp, bsave);
	m_ms.ShallowCopy(dmp, bsave);
}

//-----------------------------------------------------------------------------
void FESlidingInterfaceMP::ContactForces(FEGlobalVector& R)
{
	int i, j, k;
	vector<int> sLM, mLM, LM, en;
	vector<double> fe;
	const int MN = FEElement::MAX_NODES;
	double detJ[MN], w[MN], *Hs, Hm[MN];
	double N[MN*10];
	int nsol = (int)m_sid.size();
	
	FEModel& fem = *GetFEModel();
	
	double dt = fem.GetCurrentStep()->m_dt;
	
	// loop over the nr of passes
	int npass = (m_btwo_pass?2:1);
	for (int np=0; np<npass; ++np)
	{
		// get slave and master surface
		FESlidingSurfaceMP& ss = (np == 0? m_ss : m_ms);
		FESlidingSurfaceMP& ms = (np == 0? m_ms : m_ss);
		vector<int>& sl = (np == 0? m_ssl : m_msl);
		
		// loop over all slave elements
		for (i=0; i<ss.Elements(); ++i)
		{
			// get the surface element
			FESurfaceElement& se = ss.Element(i);
			
			bool sporo = ss.m_bporo;
			
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
				FESlidingSurfaceMP::Data& pt = ss.m_Data[i][j];
				// get the master element
				FESurfaceElement* pme = pt.m_pme;
				if (pme)
				{
					// get the master element
					FESurfaceElement& me = *pme;
					
					bool mporo = ms.m_bporo;
					
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
					// only do this when the node is actually in contact
					// in other words, when tn > 0
					if (tn > 0)
					{
						if (sporo && mporo)
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
						for (int isol=0; isol<nsol; ++isol)
						{
							int sid = m_sid[isol];
							int l = sl[isol];
							
							// calculate nr of concentration dofs
							int ndof = nseln + nmeln;
							
							// calculate the flow rate
							double epsc = m_epsc*pt.m_epsc[l];
							
							double jn = pt.m_Lmc[l] + epsc*pt.m_cg[l];
							
							// fill the LM
							LM.resize(ndof);
							for (k=0; k<nseln; ++k) LM[k        ] = sLM[(DOF_C+sid)*nseln+k];
							for (k=0; k<nmeln; ++k) LM[k + nseln] = mLM[(DOF_C+sid)*nmeln+k];
							
							// fill the force array
							fe.resize(ndof);
							zero(fe);
							for (k=0; k<nseln; ++k) N[k      ] =  Hs[k];
							for (k=0; k<nmeln; ++k) N[k+nseln] = -Hm[k];
							
							for (k=0; k<ndof; ++k) fe[k] += dt*jn*N[k]*detJ[j]*w[j];
							
							// assemble residual
							R.Assemble(en, LM, fe);
						}
					}
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------
void FESlidingInterfaceMP::ContactStiffness(FESolver* psolver)
{
	int i, j, k, l;
	vector<int> sLM, mLM, LM, en;
	const int MN = FEElement::MAX_NODES;
	double detJ[MN], w[MN], *Hs, Hm[MN];
	double pt[MN], dpr[MN], dps[MN];
	double N[MN*10];
	matrix ke;
	int nsol = (int)m_sid.size();
	vector< vector<double> > ct(nsol,vector<double>(MN));
	vector< vector<double> > dcr(nsol,vector<double>(MN));
	vector< vector<double> > dcs(nsol,vector<double>(MN));
	
	FEModel& fem = *GetFEModel();
 	
	// see how many reformations we've had to do so far
	int nref = psolver->m_nref;
	
	// set higher order stiffness mutliplier
	// NOTE: this algorithm doesn't really need this
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
		FESlidingSurfaceMP& ss = (np == 0? m_ss : m_ms);
		FESlidingSurfaceMP& ms = (np == 0? m_ms : m_ss);
		vector<int>& sl = (np == 0? m_ssl : m_msl);
		
		// loop over all slave elements
		for (i=0; i<ss.Elements(); ++i)
		{
			// get ths slave element
			FESurfaceElement& se = ss.Element(i);
			
			bool sporo = ss.m_bporo;
			
			// get nr of nodes and integration points
			int nseln = se.Nodes();
			int nint = se.GaussPoints();
			
			double pn[MN];
			vector< vector<double> >cn(nsol,vector<double>(MN));
			for (j=0; j<nseln; ++j)
			{
				pn[j] = ss.GetMesh()->Node(se.m_node[j]).m_pt;
				for (int isol=0; isol<nsol; ++isol) {
					cn[isol][j] = ss.GetMesh()->Node(se.m_node[j]).m_ct[m_sid[isol]];
				}
			}
			
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
				// concentration
				for (int isol=0; isol<nsol; ++isol) {
					ct[isol][j] = se.eval(&cn[isol][0], j);
					dcr[isol][j] = se.eval_deriv1(&cn[isol][0], j);
					dcs[isol][j] = se.eval_deriv2(&cn[isol][0], j);
				}
			}
			
			// loop over all integration points
			for (j=0; j<nint; ++j)
			{
				FESlidingSurfaceMP::Data& pt = ss.m_Data[i][j];
				// get the master element
				FESurfaceElement* pme = pt.m_pme;
				if (pme)
				{
					FESurfaceElement& me = *pme;
					
					bool mporo = ms.m_bporo;
					
					// get the nr of master nodes
					int nmeln = me.Nodes();
					
					// nodal data
					double pm[MN];
					vector< vector<double> > cm(nsol,vector<double>(MN));
					for (k=0; k<nmeln; ++k) 
					{
						pm[k] = ms.GetMesh()->Node(me.m_node[k]).m_pt;
						for (int isol=0; isol<nsol; ++isol) {
							cm[isol][k] = ms.GetMesh()->Node(me.m_node[k]).m_ct[m_sid[isol]];
						}
					}
					
					// copy the LM vector
					ms.UnpackLM(me, mLM);
					
					int ndpn;	// number of dofs per node
					int ndof;	// number of dofs in stiffness matrix
					
					if (nsol) {
						// calculate dofs for biphasic-solute contact
						ndpn = 4+nsol;
						ndof = ndpn*(nseln+nmeln);
						
						// build the LM vector
						LM.resize(ndof);
						
						for (k=0; k<nseln; ++k)
						{
							LM[ndpn*k  ] = sLM[3*k  ];			// x-dof
							LM[ndpn*k+1] = sLM[3*k+1];			// y-dof
							LM[ndpn*k+2] = sLM[3*k+2];			// z-dof
							LM[ndpn*k+3] = sLM[3*nseln+k];		// p-dof
							for (int isol=0; isol<nsol; ++isol)
								LM[ndpn*k+4+isol] = sLM[(DOF_C+m_sid[isol])*nseln+k];		// c-dof
						}
						for (k=0; k<nmeln; ++k)
						{
							LM[ndpn*(k+nseln)  ] = mLM[3*k  ];			// x-dof
							LM[ndpn*(k+nseln)+1] = mLM[3*k+1];			// y-dof
							LM[ndpn*(k+nseln)+2] = mLM[3*k+2];			// z-dof
							LM[ndpn*(k+nseln)+3] = mLM[3*nmeln+k];		// p-dof
							for (int isol=0; isol<nsol; ++isol)
								LM[ndpn*(k+nseln)+4+isol] = mLM[(DOF_C+m_sid[isol])*nmeln+k];		// c-dof
						}
					}
					
					else if (sporo && mporo) {
						// calculate dofs for biphasic contact
						ndpn = 4;
						ndof = ndpn*(nseln+nmeln);
						
						// build the LM vector
						LM.resize(ndof);
						
						for (k=0; k<nseln; ++k)
						{
							LM[ndpn*k  ] = sLM[3*k  ];			// x-dof
							LM[ndpn*k+1] = sLM[3*k+1];			// y-dof
							LM[ndpn*k+2] = sLM[3*k+2];			// z-dof
							LM[ndpn*k+3] = sLM[3*nseln+k];		// p-dof
						}
						for (k=0; k<nmeln; ++k)
						{
							LM[ndpn*(k+nseln)  ] = mLM[3*k  ];			// x-dof
							LM[ndpn*(k+nseln)+1] = mLM[3*k+1];			// y-dof
							LM[ndpn*(k+nseln)+2] = mLM[3*k+2];			// z-dof
							LM[ndpn*(k+nseln)+3] = mLM[3*nmeln+k];		// p-dof
						}
					}
					
					else {
						// calculate dofs for elastic contact
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
					
					//	double dtn = m_eps*HEAVYSIDE(Lm + eps*g);
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
					
					if (ndpn > 3) {
						for (k=0; k<nseln; ++k)
						{
							N[ndpn*k+3] = 0;
							for (int isol=0; isol<nsol; ++isol)
								N[ndpn*k+4+isol] = 0;
						}
						for (k=0; k<nmeln; ++k)
						{
							N[ndpn*(k+nseln)+3] = 0;
							for (int isol=0; isol<nsol; ++isol)
								N[ndpn*(k+nseln)+4+isol] = 0;
						}
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
					
					// --- M U L T I P H A S I C   S T I F F N E S S ---
					if ((tn > 0) && nsol)
					{
						double dt = fem.GetCurrentStep()->m_dt;
						
						double epsp = m_epsp*pt.m_epsp;
												
						// --- S O L I D - P R E S S U R E / S O L U T E   C O N T A C T ---
						
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
									ke[ndpn*k + 3][ndpn*l  ] += dt*w[j]*detJ[j]*epsp*N[k]*N[l]*(dpmr*Gm[0].x + dpms*Gm[1].x);
									ke[ndpn*k + 3][ndpn*l+1] += dt*w[j]*detJ[j]*epsp*N[k]*N[l]*(dpmr*Gm[0].y + dpms*Gm[1].y);
									ke[ndpn*k + 3][ndpn*l+2] += dt*w[j]*detJ[j]*epsp*N[k]*N[l]*(dpmr*Gm[0].z + dpms*Gm[1].z);
									
									for (int isol=0; isol<nsol; ++isol) {
										double epsc = m_epsc*pt.m_epsc[sl[isol]];
										double dcmr = me.eval_deriv1(&cm[isol][0], r, s);
										double dcms = me.eval_deriv2(&cm[isol][0], r, s);
										
										ke[ndpn*k + 4 + isol][ndpn*l  ] += dt*w[j]*detJ[j]*epsc*N[k]*N[l]*(dcmr*Gm[0].x + dcms*Gm[1].x);
										ke[ndpn*k + 4 + isol][ndpn*l+1] += dt*w[j]*detJ[j]*epsc*N[k]*N[l]*(dcmr*Gm[0].y + dcms*Gm[1].y);
										ke[ndpn*k + 4 + isol][ndpn*l+2] += dt*w[j]*detJ[j]*epsc*N[k]*N[l]*(dcmr*Gm[0].z + dcms*Gm[1].z);
									}
								}
							
							double wn = pt.m_Lmp + epsp*pt.m_pg;
							
							// b. A-term
							//-------------------------------------
							
							for (l=0; l<nseln; ++l)
								for (k=0; k<nseln+nmeln; ++k)
								{
									ke[ndpn*k + 3][ndpn*l  ] -= dt*w[j]*wn*N[k]*(As[l][0][0]*nu.x + As[l][0][1]*nu.y + As[l][0][2]*nu.z);
									ke[ndpn*k + 3][ndpn*l+1] -= dt*w[j]*wn*N[k]*(As[l][1][0]*nu.x + As[l][1][1]*nu.y + As[l][1][2]*nu.z);
									ke[ndpn*k + 3][ndpn*l+2] -= dt*w[j]*wn*N[k]*(As[l][2][0]*nu.x + As[l][2][1]*nu.y + As[l][2][2]*nu.z);
									
									for (int isol=0; isol<nsol; ++isol) {
										double epsc = m_epsc*pt.m_epsc[sl[isol]];
										double jn = pt.m_Lmc[sl[isol]] + epsc*pt.m_cg[sl[isol]];
										ke[ndpn*k + 4 + isol][ndpn*l  ] -= dt*w[j]*jn*N[k]*(As[l][0][0]*nu.x + As[l][0][1]*nu.y + As[l][0][2]*nu.z);
										ke[ndpn*k + 4 + isol][ndpn*l+1] -= dt*w[j]*jn*N[k]*(As[l][1][0]*nu.x + As[l][1][1]*nu.y + As[l][1][2]*nu.z);
										ke[ndpn*k + 4 + isol][ndpn*l+2] -= dt*w[j]*jn*N[k]*(As[l][2][0]*nu.x + As[l][2][1]*nu.y + As[l][2][2]*nu.z);
									}
								}
							
							// c. m-term
							//---------------------------------------
							
							for (k=0; k<nmeln; ++k)
								for (l=0; l<nseln+nmeln; ++l)
								{
									ke[ndpn*(k+nseln) + 3][ndpn*l  ] += dt*w[j]*detJ[j]*wn*N[l]*mm[k].x;
									ke[ndpn*(k+nseln) + 3][ndpn*l+1] += dt*w[j]*detJ[j]*wn*N[l]*mm[k].y;
									ke[ndpn*(k+nseln) + 3][ndpn*l+2] += dt*w[j]*detJ[j]*wn*N[l]*mm[k].z;
									
									for (int isol=0; isol<nsol; ++isol) {
										double epsc = m_epsc*pt.m_epsc[sl[isol]];
										double jn = pt.m_Lmc[sl[isol]] + epsc*pt.m_cg[sl[isol]];
										ke[ndpn*(k+nseln) + 4 + isol][ndpn*l  ] += dt*w[j]*detJ[j]*jn*N[l]*mm[k].x;
										ke[ndpn*(k+nseln) + 4 + isol][ndpn*l+1] += dt*w[j]*detJ[j]*jn*N[l]*mm[k].y;
										ke[ndpn*(k+nseln) + 4 + isol][ndpn*l+2] += dt*w[j]*detJ[j]*jn*N[l]*mm[k].z;
									}
								}
						}
						
						
						// --- P R E S S U R E - P R E S S U R E   C O N T A C T ---
						
						// calculate the N-vector
						for (k=0; k<nseln; ++k)
						{
							N[ndpn*k+3] = Hs[k];
							for (int isol=0; isol<nsol; ++isol)
								N[ndpn*k+4+isol] = Hs[k];
						}
						
						for (k=0; k<nmeln; ++k)
						{
							N[ndpn*(k+nseln)+3] = -Hm[k];
							for (int isol=0; isol<nsol; ++isol)
								N[ndpn*(k+nseln)+4+isol] = -Hm[k];
						}
						
						for (k=3; k<ndof; k+=ndpn)
							for (l=3; l<ndof; l+=ndpn) ke[k][l] -= dt*epsp*w[j]*detJ[j]*N[k]*N[l];
						
						// --- C O N C E N T R A T I O N - C O N C E N T R A T I O N   C O N T A C T ---
						
						for (int isol=0; isol<nsol; ++isol) {
							double epsc = m_epsc*pt.m_epsc[sl[isol]];
							for (k=4+isol; k<ndof; k+=ndpn)
								for (l=4+isol; l<ndof; l+=ndpn) ke[k][l] -= dt*epsc*w[j]*detJ[j]*N[k]*N[l];
						}
						
					}
					
					// --- B I P H A S I C   S T I F F N E S S ---
					else if ((tn > 0) && (sporo && mporo))
					{
						// the variable dt is either the timestep or one
						// depending on whether we are using the symmetric
						// poro version or not.
						double dt = fem.GetCurrentStep()->m_dt;
						
						double epsp = m_epsp*pt.m_epsp;
						
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
									ke[ndpn*k + 3][ndpn*l  ] += dt*w[j]*detJ[j]*epsp*N[k]*N[l]*(dpmr*Gm[0].x + dpms*Gm[1].x);
									ke[ndpn*k + 3][ndpn*l+1] += dt*w[j]*detJ[j]*epsp*N[k]*N[l]*(dpmr*Gm[0].y + dpms*Gm[1].y);
									ke[ndpn*k + 3][ndpn*l+2] += dt*w[j]*detJ[j]*epsp*N[k]*N[l]*(dpmr*Gm[0].z + dpms*Gm[1].z);
								}
							
							double wn = pt.m_Lmp + epsp*pt.m_pg;
							
							// b. A-term
							//-------------------------------------
							
							for (l=0; l<nseln; ++l)
								for (k=0; k<nseln+nmeln; ++k)
								{
									ke[ndpn*k + 3][ndpn*l  ] -= dt*w[j]*wn*N[k]*(As[l][0][0]*nu.x + As[l][0][1]*nu.y + As[l][0][2]*nu.z);
									ke[ndpn*k + 3][ndpn*l+1] -= dt*w[j]*wn*N[k]*(As[l][1][0]*nu.x + As[l][1][1]*nu.y + As[l][1][2]*nu.z);
									ke[ndpn*k + 3][ndpn*l+2] -= dt*w[j]*wn*N[k]*(As[l][2][0]*nu.x + As[l][2][1]*nu.y + As[l][2][2]*nu.z);
								}
							
							// c. m-term
							//---------------------------------------
							
							for (k=0; k<nmeln; ++k)
								for (l=0; l<nseln+nmeln; ++l)
								{
									ke[ndpn*(k+nseln) + 3][ndpn*l  ] += dt*w[j]*detJ[j]*wn*N[l]*mm[k].x;
									ke[ndpn*(k+nseln) + 3][ndpn*l+1] += dt*w[j]*detJ[j]*wn*N[l]*mm[k].y;
									ke[ndpn*(k+nseln) + 3][ndpn*l+2] += dt*w[j]*detJ[j]*wn*N[l]*mm[k].z;
								}
						}
						
						
						// --- P R E S S U R E - P R E S S U R E   C O N T A C T ---
						
						// calculate the N-vector
						for (k=0; k<nseln; ++k)
						{
							N[ndpn*k+3] = Hs[k];
						}
						
						for (k=0; k<nmeln; ++k)
						{
							N[ndpn*(k+nseln)+3] = -Hm[k];
						}
						
						for (k=3; k<ndof; k+=ndpn)
							for (l=3; l<ndof; l+=ndpn) ke[k][l] -= dt*epsp*w[j]*detJ[j]*N[k]*N[l];
						
					}
					
					// assemble the global stiffness
					psolver->AssembleStiffness(en, LM, ke);
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------
void FESlidingInterfaceMP::UpdateContactPressures()
{
	int np, n, i, j;
	const int MN = FEElement::MAX_NODES;
	int npass = (m_btwo_pass?2:1);
	for (np=0; np<npass; ++np)
	{
		FESlidingSurfaceMP& ss = (np == 0? m_ss : m_ms);
		FESlidingSurfaceMP& ms = (np == 0? m_ms : m_ss);
		
		// loop over all elements of the primary surface
		for (n=0; n<ss.Elements(); ++n)
		{
			FESurfaceElement& el = ss.Element(n);
			int nint = el.GaussPoints();
			
			// get the normal tractions at the integration points
			double gap, eps;
			for (i=0; i<nint; ++i)
			{
				FESlidingSurfaceMP::Data& pt = ss.m_Data[n][i];
				gap = pt.m_gap;
				eps = m_epsn*pt.m_epsn;
				pt.m_Ln = MBRACKET(pt.m_Lmd + eps*gap);
				FESurfaceElement* pme = pt.m_pme;
				if (m_btwo_pass && pme)
				{
					int mint = pme->GaussPoints();
					vector<FESlidingSurfaceMP::Data>& md = ms.m_Data[pme->m_lid];
					double ti[MN];
					for (j=0; j<mint; ++j) {
						gap = md[j].m_gap;
						eps = m_epsn*md[j].m_epsn;
						ti[j] = MBRACKET(md[j].m_Lmd + m_epsn*md[j].m_epsn*md[j].m_gap);
					}
					// project the data to the nodes
					double tn[MN];
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
bool FESlidingInterfaceMP::Augment(int naug)
{
	// make sure we need to augment
	if (!m_blaugon) return true;
	
	double Ln, Lp;
	int nsol = (int)m_sid.size();
	vector<double>Lc(nsol);
	bool bconv = true;
	
	bool bporo = (m_ss.m_bporo && m_ms.m_bporo);
	bool bsolu = (m_ss.m_bsolu && m_ms.m_bsolu);
	int NS = (int)m_ss.m_Data.size();
	int NM = (int)m_ms.m_Data.size();
	
	// --- c a l c u l a t e   i n i t i a l   n o r m s ---
	// a. normal component
	double normL0 = 0, normP = 0, normDP = 0, normC = 0;
	vector<double>normDC(nsol,0);
	for (int i=0; i<NS; ++i)
	{
		vector<FESlidingSurfaceMP::Data>& sd = m_ss.m_Data[i];
		for (int j=0; j<(int)sd.size(); ++j)
		{
			FESlidingSurfaceMP::Data& ds = sd[j];
			normL0 += ds.m_Lmd*ds.m_Lmd;
		}
	}
	for (int i=0; i<NM; ++i)
	{
		vector<FESlidingSurfaceMP::Data>& md = m_ms.m_Data[i];
		for (int j=0; j<(int)md.size(); ++j)
		{
			FESlidingSurfaceMP::Data& dm = md[j];
			normL0 += dm.m_Lmd*dm.m_Lmd;
		}
	}
	
	// b. gap component
	// (is calculated during update)
	double maxgap = 0, maxpg = 0;
	vector<double> maxcg(nsol,0);
	
	// update Lagrange multipliers
	double normL1 = 0, eps, epsp, epsc;
	for (int i=0; i<NS; ++i)
	{
		vector<FESlidingSurfaceMP::Data>& sd = m_ss.m_Data[i];
		for (int j=0; j<(int)sd.size(); ++j)
		{
			FESlidingSurfaceMP::Data& ds = sd[j];
            
            // update Lagrange multipliers on slave surface
			eps = m_epsn*ds.m_epsn;
			Ln = ds.m_Lmd + eps*ds.m_gap;
			ds.m_Lmd = MBRACKET(Ln);
            
			normL1 += ds.m_Lmd*ds.m_Lmd;
            
            if (m_ss.m_bporo) {
                Lp = 0;
                Lc.assign(nsol, 0);
                if (Ln > 0) {
					epsp = m_epsp*ds.m_epsp;
					Lp = ds.m_Lmp + epsp*ds.m_pg;
					maxpg = max(maxpg,fabs(ds.m_pg));
					normDP += ds.m_pg*ds.m_pg;
                    for (int isol=0; isol<nsol; ++isol) {
                        int l = m_ssl[isol];
                        epsc = m_epsc*ds.m_epsc[l];
                        Lc[isol] = ds.m_Lmc[l] + epsc*ds.m_cg[l];
                        maxcg[isol] = max(maxcg[isol],fabs(ds.m_cg[l]));
                        normDC[isol] += ds.m_cg[l]*ds.m_cg[l];
                    }
                }
                ds.m_Lmp = Lp;
                for (int isol=0; isol<nsol; ++isol) ds.m_Lmc[m_ssl[isol]] = Lc[isol];
            }
            
            if (Ln > 0) maxgap = max(maxgap,fabs(ds.m_gap));
        }
	}
	
	for (int i=0; i<NM; ++i)
	{
		vector<FESlidingSurfaceMP::Data>& md = m_ms.m_Data[i];
		for (int j=0; j<(int)md.size(); ++j)
		{
			FESlidingSurfaceMP::Data& dm = md[j];
            
            // update Lagrange multipliers on master surface
			eps = m_epsn*dm.m_epsn;
			Ln = dm.m_Lmd + eps*dm.m_gap;
			dm.m_Lmd = MBRACKET(Ln);
            
			normL1 += dm.m_Lmd*dm.m_Lmd;
            
            if (m_ms.m_bporo) {
                Lp = 0;
                Lc.assign(nsol, 0);
                if (Ln > 0) {
					epsp = m_epsp*dm.m_epsp;
					Lp = dm.m_Lmp + epsp*dm.m_pg;
					maxpg = max(maxpg,fabs(dm.m_pg));
					normDP += dm.m_pg*dm.m_pg;
                    for (int isol=0; isol<nsol; ++isol) {
                        int l = m_msl[isol];
                        epsc = m_epsc*dm.m_epsc[l];
                        Lc[isol] = dm.m_Lmc[l] + epsc*dm.m_cg[l];
                        maxcg[isol] = max(maxcg[isol],fabs(dm.m_cg[l]));
                        normDC[isol] += dm.m_cg[l]*dm.m_cg[l];
                    }
                }
                dm.m_Lmp = Lp;
                for (int isol=0; isol<nsol; ++isol) dm.m_Lmc[m_msl[isol]] = Lc[isol];
            }
            
            if (Ln > 0) maxgap = max(maxgap,fabs(dm.m_gap));
        }
	}
	
	// normP should be a measure of the fluid pressure at the
	// contact interface.  However, since it could be zero,
	// use an average measure of the contact traction instead.
	normP = normL1;
	normC = normL1/(m_Rgas*m_Tabs);
	
	// calculate relative norms
	double lnorm = (normL1 != 0 ? fabs((normL1 - normL0) / normL1) : fabs(normL1 - normL0)); 
	double pnorm = (normP != 0 ? (normDP/normP) : normDP);
	vector<double> cnorm(nsol);
	for (int isol=0; isol<nsol; ++isol)
		cnorm[isol] = (normC != 0 ? (normDC[isol]/normC) : normDC[isol]);
	
	// check convergence
	if ((m_gtol > 0) && (maxgap > m_gtol)) bconv = false;
	if ((m_ptol > 0) && (bporo && maxpg > m_ptol)) bconv = false;
	for (int isol=0; isol<nsol; ++isol)
		if ((m_ctol > 0) && (bsolu && maxcg[isol] > m_ctol)) bconv = false;
	
	if ((m_atol > 0) && (lnorm > m_atol)) bconv = false;
	if ((m_atol > 0) && (pnorm > m_atol)) bconv = false;
	for (int isol=0; isol<nsol; ++isol)
		if ((m_atol > 0) && (cnorm[isol] > m_atol)) bconv = false;
	
	if (naug < m_naugmin ) bconv = false;
	if (naug >= m_naugmax) bconv = true;
	
	felog.printf(" sliding interface # %d\n", m_nID);
	felog.printf("                        CURRENT        REQUIRED\n");
	felog.printf("    D multiplier : %15le", lnorm);
	if (m_atol > 0) felog.printf("%15le\n", m_atol);
	else felog.printf("       ***\n");
	if (bporo) { felog.printf("    P gap       : %15le", pnorm);
		if (m_atol > 0) felog.printf("%15le\n", m_atol);
		else felog.printf("       ***\n");
	}
	for (int isol=0; isol<nsol; ++isol) {
		felog.printf("    C[%d] gap   : %15le", m_sid[isol], cnorm[isol]);
		if (m_atol > 0) felog.printf("%15le\n", m_atol);
		else felog.printf("       ***\n");
	}
	
	felog.printf("    maximum gap  : %15le", maxgap);
	if (m_gtol > 0) felog.printf("%15le\n", m_gtol); else felog.printf("       ***\n");
	if (bporo) {
		felog.printf("    maximum pgap : %15le", maxpg);
		if (m_ptol > 0) felog.printf("%15le\n", m_ptol); else felog.printf("       ***\n");
	}
	for (int isol=0; isol<nsol; ++isol) {
		felog.printf("    maximum cgap[%d] : %15le", m_sid[isol], maxcg[isol]);
		if (m_ctol > 0) felog.printf("%15le\n", m_ctol); else felog.printf("       ***\n");
	}
	
	return bconv;
}

//-----------------------------------------------------------------------------
void FESlidingInterfaceMP::Serialize(DumpFile &ar)
{
	// serialize contact data
	FEContactInterface::Serialize(ar);
	
	// serialize contact surface data
	m_ms.Serialize(ar);
	m_ss.Serialize(ar);
}

//-----------------------------------------------------------------------------

void FESlidingInterfaceMP::MarkAmbient()
{	
	int i, j, id, np;
	
	// Mark all nodes as free-draining.  This needs to be done for ALL
	// contact interfaces prior to executing Update(), where nodes that are
	// in contact are subsequently marked as non free-draining.  This ensures
	// that for surfaces involved in more than one contact interface, nodes
	// that have been marked as non free-draining are not reset to 
	// free-draining.
	for (np=0; np<2; ++np)
	{
		FESlidingSurfaceMP& s = (np == 0? m_ss : m_ms);
		
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
		if (s.m_bsolu) {
			// first, mark all nodes as free-draining (= neg. ID)
			// this is done by setting the dof's equation number
			// to a negative number
			for (i=0; i<s.Nodes(); ++i) 
			{
				for (j=0; j<MAX_CDOFS; ++j) {
					id = s.Node(i).m_ID[DOF_C+j];
					if (id >= 0) 
					{
						FENode& node = s.Node(i);
						// mark node as free-draining
						node.m_ID[DOF_C+j] = -id-2;
					}
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------

void FESlidingInterfaceMP::SetAmbient()
{	
	int i, j, np;
	
	// Set the pressure to zero for the free-draining nodes
	for (np=0; np<2; ++np)
	{
		FESlidingSurfaceMP& s = (np == 0? m_ss : m_ms);
		
		if (s.m_bporo) {
			// loop over all nodes
			for (i=0; i<s.Nodes(); ++i) 
			{
				if (s.Node(i).m_ID[DOF_P] < -1)
				{
					FENode& node = s.Node(i);
					// set the fluid pressure to ambient condition
					node.m_pt = m_ambp;
				}
			}
		}
		if (s.m_bsolu) {
			// loop over all nodes
			for (i=0; i<s.Nodes(); ++i) 
			{
				for (j=0; j<MAX_CDOFS; ++j) {
					if (s.Node(i).m_ID[DOF_C+j] < -1)
					{
						FENode& node = s.Node(i);
						// set the fluid pressure to ambient condition
						node.m_ct[j] = m_ambc[j];
					}
				}
			}
		}
	}
}
