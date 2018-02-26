//
//  FESlidingInterfaceBiphasic.cpp
//  FEBioMix
//
//  Created by Gerard Ateshian on 5/1/16.
//  Copyright © 2016 febio.org. All rights reserved.
//

#include "FESlidingInterfaceBiphasic.h"
#include "FEBiphasic.h"
#include "FECore/FEModel.h"
#include "FECore/FEAnalysis.h"
#include "FECore/FENormalProjection.h"
#include "FECore/log.h"

//-----------------------------------------------------------------------------
// Define sliding interface parameters
BEGIN_PARAMETER_LIST(FESlidingInterfaceBiphasic, FEContactInterface)
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
ADD_PARAMETER(m_breloc   , FE_PARAM_BOOL  , "node_reloc"         );
ADD_PARAMETER(m_mu       , FE_PARAM_DOUBLE, "fric_coeff"         );
ADD_PARAMETER(m_phi      , FE_PARAM_DOUBLE, "contact_frac"       );
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FESlidingSurfaceBiphasic::Data::Data()
{
    m_gap  = 0.0;
    m_Lmd  = 0.0;
    m_Lmt  = m_tr = vec3d(0,0,0);
    m_Lmp  = 0.0;
    m_epsn = 1.0;
    m_epsp = 1.0;
    m_Ln   = 0.0;
    m_pg   = 0.0;
    m_p1   = 0;
    m_nu   = m_s1 = m_dg = vec3d(0,0,0);
    m_rs   = m_rsp = vec2d(0,0);
    m_bstick = false;
    m_pme  = m_pmep = (FESurfaceElement*)0;
}

//-----------------------------------------------------------------------------
// FESlidingSurfaceBiphasic
//-----------------------------------------------------------------------------

FESlidingSurfaceBiphasic::FESlidingSurfaceBiphasic(FEModel* pfem) : FEBiphasicContactSurface(pfem)
{
    m_bporo = false;
    m_pfem = pfem;
}

//-----------------------------------------------------------------------------
bool FESlidingSurfaceBiphasic::Init()
{
    // initialize surface data first
    if (FEBiphasicContactSurface::Init() == false) return false;
    
    // allocate data structures
    int NE = Elements();
    m_Data.resize(NE);
    for (int i=0; i<NE; ++i)
    {
        FESurfaceElement& el = Element(i);
        int nint = el.GaussPoints();
        m_Data[i].resize(nint);
    }
    
    // allocate node normals and contact tractions
    m_nn.assign(Nodes(), vec3d(0,0,0));
    m_tn.assign(Nodes(), vec3d(0,0,0));
    
    // determine biphasic status
    m_poro.resize(Elements(),false);
    for (int i=0; i<Elements(); ++i)
    {
        // get the surface element
        FESurfaceElement& se = Element(i);
        
        // get the element this surface element belongs to
        FEElement* pe = m_pMesh->FindElementFromID(se.m_elem[0]);
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
void FESlidingSurfaceBiphasic::InitSlidingSurface()
{
    for (int i=0; i<Elements(); ++i)
    {
        FESurfaceElement& el = Element(i);
        int nint = el.GaussPoints();
        for (int j=0; j<nint; ++j)
        {
            // Store current surface projection values as previous
            // Set stick-slip switching counter to zero
            m_Data[i][j].m_rsp = m_Data[i][j].m_rs;
            m_Data[i][j].m_pmep = m_Data[i][j].m_pme;
        }
    }
}

//-----------------------------------------------------------------------------
//! Evaluate the nodal contact tractions by averaging values from surrounding
//! faces.  This function ensures that nodal contact tractions are always
//! compressive, so that they can be used to detect free-draining status.

void FESlidingSurfaceBiphasic::EvaluateNodalContactTractions()
{
    const int N = Nodes();
    
    // number of faces with non-zero contact pressure connected to this node
    vector<int> nfaces(N,0);
    
    // zero nodal contact tractions
    zero(m_tn);
    
    // loop over all elements
    for (int i=0; i<Elements(); ++i)
    {
        FESurfaceElement& el = Element(i);
        int ne = el.Nodes();
        
        // get the average contact pressure for that face
        vec3d tn(0,0,0);
        GetContactTraction(i, tn);
        
        if (tn.norm() > 0) {
            for (int j=0; j<ne; ++j)
            {
                m_tn[el.m_lnode[j]] += tn;
                ++nfaces[el.m_lnode[j]];
            }
        }
    }
    
    // get average over all contacting faces sharing that node
    for (int i=0; i<N; ++i)
        if (nfaces[i] > 0) m_tn[i] /= nfaces[i];
}

//-----------------------------------------------------------------------------
//! This function calculates the node normal. Due to the piecewise continuity
//! of the surface elements this normal is not uniquely defined so in order to
//! obtain a unique normal the normal is averaged for each node over all the
//! element normals at the node

void FESlidingSurfaceBiphasic::UpdateNodeNormals()
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
vec3d FESlidingSurfaceBiphasic::GetContactForce()
{
    return m_Ft;
}

//-----------------------------------------------------------------------------
double FESlidingSurfaceBiphasic::GetContactArea()
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
			if (data.m_Ln > 0)
			{
	            // get the base vectors
		        vec3d g[2];
			    CoBaseVectors(el, i, g);
            
				// normal (magnitude = area)
				vec3d n = g[0] ^ g[1];
            
				// gauss weight
				double w = el.GaussWeights()[i];
            
				// contact force
				a += n.norm()*w;
			}
        }
    }
    
    return a;
}

//-----------------------------------------------------------------------------
vec3d FESlidingSurfaceBiphasic::GetFluidForce()
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
        for (i=0; i<nseln; ++i) pn[i] = GetMesh()->Node(el.m_node[i]).get(m_dofP);
        
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
void FESlidingSurfaceBiphasic::Serialize(DumpStream& ar)
{
    if (ar.IsShallow())
    {
        FEContactSurface::Serialize(ar);
        if (ar.IsSaving())
        {
            ar << m_bporo;
            
            for (int i=0; i<(int) m_Data.size(); ++i)
            {
                vector<Data>& di = m_Data[i];
                int nint = (int) di.size();
                for (int j=0; j<nint; ++j)
                {
                    Data& d = di[j];
                    ar << d.m_gap;
                    ar << d.m_dg;
                    ar << d.m_nu;
                    ar << d.m_s1;
                    ar << d.m_rs;
                    ar << d.m_rsp;
                    ar << d.m_Lmd;
                    ar << d.m_Lmt;
                    ar << d.m_Lmp;
                    ar << d.m_epsn;
                    ar << d.m_epsp;
                    ar << d.m_pg;
                    ar << d.m_p1;
                    ar << d.m_Ln;
                    ar << d.m_bstick;
                    ar << d.m_tr;
                }
            }
            ar << m_poro;
            ar << m_nn;
        }
        else
        {
            ar >> m_bporo;
            
            for (int i=0; i<(int) m_Data.size(); ++i)
            {
                vector<Data>& di = m_Data[i];
                int nint = (int) di.size();
                for (int j=0; j<nint; ++j)
                {
                    Data& d = di[j];
                    ar >> d.m_gap;
                    ar >> d.m_dg;
                    ar >> d.m_nu;
                    ar >> d.m_s1;
                    ar >> d.m_rs;
                    ar >> d.m_rsp;
                    ar >> d.m_Lmd;
                    ar >> d.m_Lmt;
                    ar >> d.m_Lmp;
                    ar >> d.m_epsn;
                    ar >> d.m_epsp;
                    ar >> d.m_pg;
                    ar >> d.m_p1;
                    ar >> d.m_Ln;
                    ar >> d.m_bstick;
                    ar >> d.m_tr;
                }
            }
            ar >> m_poro;
            ar >> m_nn;
        }
    }
    else
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
                    ar << d.m_dg;
                    ar << d.m_nu;
                    ar << d.m_s1;
                    ar << d.m_rs;
                    ar << d.m_rsp;
                    ar << d.m_Lmd;
                    ar << d.m_Lmt;
                    ar << d.m_Lmp;
                    ar << d.m_epsn;
                    ar << d.m_epsp;
                    ar << d.m_pg;
                    ar << d.m_p1;
                    ar << d.m_Ln;
                    ar << d.m_bstick;
                    ar << d.m_tr;
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
                    ar >> d.m_dg;
                    ar >> d.m_nu;
                    ar >> d.m_s1;
                    ar >> d.m_rs;
                    ar >> d.m_rsp;
                    ar >> d.m_Lmd;
                    ar >> d.m_Lmt;
                    ar >> d.m_Lmp;
                    ar >> d.m_epsn;
                    ar >> d.m_epsp;
                    ar >> d.m_pg;
                    ar >> d.m_p1;
                    ar >> d.m_Ln;
                    ar >> d.m_bstick;
                    ar >> d.m_tr;
                }
            }
            ar >> m_poro;
            ar >> m_nn;
        }
    }
}

//-----------------------------------------------------------------------------
void FESlidingSurfaceBiphasic::GetContactGap(int nface, double& pg)
{
    FESurfaceElement& el = Element(nface);
    int ni = el.GaussPoints();
    pg = 0;
    for (int k=0; k<ni; ++k) pg += m_Data[nface][k].m_gap;
    pg /= ni;
}

//-----------------------------------------------------------------------------
void FESlidingSurfaceBiphasic::GetContactPressure(int nface, double& pg)
{
    FESurfaceElement& el = Element(nface);
    int ni = el.GaussPoints();
    pg = 0;
    for (int k=0; k<ni; ++k) pg += m_Data[nface][k].m_Ln;
    pg /= ni;
}

//-----------------------------------------------------------------------------
void FESlidingSurfaceBiphasic::GetContactTraction(int nface, vec3d& pt)
{
    FESurfaceElement& el = Element(nface);
    int ni = el.GaussPoints();
    pt = vec3d(0,0,0);
    for (int k=0; k<ni; ++k) pt += m_Data[nface][k].m_tr;
    pt /= ni;
}

//-----------------------------------------------------------------------------
void FESlidingSurfaceBiphasic::GetNodalContactGap(int nface, double* pg)
{
    FESurfaceElement& el = Element(nface);
    int ni = el.GaussPoints();
    double gi[FEElement::MAX_INTPOINTS];
    for (int k=0; k<ni; ++k) gi[k] = m_Data[nface][k].m_gap;
    el.project_to_nodes(gi, pg);
}

//-----------------------------------------------------------------------------
void FESlidingSurfaceBiphasic::GetNodalContactPressure(int nface, double* pg)
{
    FESurfaceElement& el = Element(nface);
    for (int k=0; k<el.Nodes(); ++k)
        pg[k] = m_tn[el.m_lnode[k]].norm();
}

//-----------------------------------------------------------------------------
void FESlidingSurfaceBiphasic::GetNodalPressureGap(int nface, double* pg)
{
    FESurfaceElement& el = Element(nface);
    int ni = el.GaussPoints();
    double gi[FEElement::MAX_INTPOINTS];
    for (int k=0; k<ni; ++k) gi[k] = m_Data[nface][k].m_pg;
    el.project_to_nodes(gi, pg);
}

//-----------------------------------------------------------------------------
void FESlidingSurfaceBiphasic::GetNodalContactTraction(int nface, vec3d* tn)
{
    FESurfaceElement& el = Element(nface);
    for (int k=0; k<el.Nodes(); ++k)
        tn[k] = m_tn[el.m_lnode[k]];
}

//-----------------------------------------------------------------------------
// FESlidingInterfaceBiphasic
//-----------------------------------------------------------------------------

FESlidingInterfaceBiphasic::FESlidingInterfaceBiphasic(FEModel* pfem) : FEContactInterface(pfem), m_ss(pfem), m_ms(pfem)
{
    static int count = 1;
    SetID(count++);
    
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
    m_breloc = false;
    m_mu = 0.0;
    m_phi = 0;
    
    m_naugmin = 0;
    m_naugmax = 10;
    
    m_dofP = pfem->GetDOFIndex("p");
    
    m_ss.SetSibling(&m_ms);
    m_ms.SetSibling(&m_ss);
}

//-----------------------------------------------------------------------------

FESlidingInterfaceBiphasic::~FESlidingInterfaceBiphasic()
{
}

//-----------------------------------------------------------------------------
bool FESlidingInterfaceBiphasic::Init()
{
    // initialize surface data
    if (m_ss.Init() == false) return false;
    if (m_ms.Init() == false) return false;
    
    return true;
}

//-----------------------------------------------------------------------------
void FESlidingInterfaceBiphasic::BuildMatrixProfile(FEGlobalMatrix& K)
{
    FEModel& fem = *GetFEModel();
    FEMesh& mesh = fem.GetMesh();
    
    // get the DOFS
    const int dof_X = fem.GetDOFIndex("x");
    const int dof_Y = fem.GetDOFIndex("y");
    const int dof_Z = fem.GetDOFIndex("z");
    const int dof_P = fem.GetDOFIndex("p");
    const int dof_RU = fem.GetDOFIndex("Ru");
    const int dof_RV = fem.GetDOFIndex("Rv");
    const int dof_RW = fem.GetDOFIndex("Rw");
    
    vector<int> lm(7*FEElement::MAX_NODES*2);
    
    int npass = (m_btwo_pass?2:1);
    for (int np=0; np<npass; ++np)
    {
        FESlidingSurfaceBiphasic& ss = (np == 0? m_ss : m_ms);
        
        int k, l;
        for (int j=0; j<ss.Elements(); ++j)
        {
            FESurfaceElement& se = ss.Element(j);
            int nint = se.GaussPoints();
            int* sn = &se.m_node[0];
            for (k=0; k<nint; ++k)
            {
                FESlidingSurfaceBiphasic::Data& pt = ss.m_Data[j][k];
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
                        lm[7*l  ] = id[dof_X];
                        lm[7*l+1] = id[dof_Y];
                        lm[7*l+2] = id[dof_Z];
                        lm[7*l+3] = id[dof_P];
                        lm[7*l+4] = id[dof_RU];
                        lm[7*l+5] = id[dof_RV];
                        lm[7*l+6] = id[dof_RW];
                    }
                    
                    for (l=0; l<nmeln; ++l)
                    {
                        vector<int>& id = mesh.Node(mn[l]).m_ID;
                        lm[7*(l+nseln)  ] = id[dof_X];
                        lm[7*(l+nseln)+1] = id[dof_Y];
                        lm[7*(l+nseln)+2] = id[dof_Z];
                        lm[7*(l+nseln)+3] = id[dof_P];
                        lm[7*(l+nseln)+4] = id[dof_RU];
                        lm[7*(l+nseln)+5] = id[dof_RV];
                        lm[7*(l+nseln)+6] = id[dof_RW];
                    }
                    
                    K.build_add(lm);
                }
            }
        }
    }
}

//-----------------------------------------------------------------------------
//! This function is called during the initialization
void FESlidingInterfaceBiphasic::Activate()
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
    Update(0, GetFEModel()->GetTime());
}

//-----------------------------------------------------------------------------
void FESlidingInterfaceBiphasic::CalcAutoPenalty(FESlidingSurfaceBiphasic& s)
{
    // loop over all surface elements
    for (int i=0; i<s.Elements(); ++i)
    {
        // get the surface element
        FESurfaceElement& el = s.Element(i);
        
        // calculate a penalty
        double eps = AutoPenalty(el, s);
        
        // assign to integation points of surface element
        int nint = el.GaussPoints();
        for (int j=0; j<nint; ++j)
        {
            FESlidingSurfaceBiphasic::Data& pt = s.m_Data[i][j];
            pt.m_epsn = eps;
        }
    }
}

//-----------------------------------------------------------------------------
void FESlidingInterfaceBiphasic::CalcAutoPressurePenalty(FESlidingSurfaceBiphasic& s)
{
    // loop over all surface elements
    for (int i=0; i<s.Elements(); ++i)
    {
        // get the surface element
        FESurfaceElement& el = s.Element(i);
        
        // calculate a penalty
        double eps = AutoPressurePenalty(el, s);
        
        // assign to integation points of surface element
        int nint = el.GaussPoints();
        for (int j=0; j<nint; ++j)
        {
            FESlidingSurfaceBiphasic::Data& pt = s.m_Data[i][j];
            pt.m_epsp = eps;
        }
    }
}

//-----------------------------------------------------------------------------

double FESlidingInterfaceBiphasic::AutoPressurePenalty(FESurfaceElement& el, FESlidingSurfaceBiphasic& s)
{
    // get the mesh
    FEMesh& m = GetFEModel()->GetMesh();
    
    // evaluate element surface normal at parametric center
    vec3d t[2];
    s.CoBaseVectors0(el, 0, 0, t);
    vec3d n = t[0] ^ t[1];
    n.unit();
    
    // get the element this surface element belongs to
    FEElement* pe = m.FindElementFromID(el.m_elem[0]);
    if (pe == 0) return 0.0;

    // get the material
    FEMaterial* pm = GetFEModel()->GetMaterial(pe->GetMatID());
        
    // see if this is a poro-elastic element
    FEBiphasic* biph = dynamic_cast<FEBiphasic*> (pm);
    if (biph == 0) return 0.0;

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
            
    double eps = n.x*(K[0][0]*n.x+K[0][1]*n.y+K[0][2]*n.z)
    +n.y*(K[1][0]*n.x+K[1][1]*n.y+K[1][2]*n.z)
    +n.z*(K[2][0]*n.x+K[2][1]*n.y+K[2][2]*n.z);
    
	// get the area of the surface element
	double A = s.FaceArea(el);

	// get the volume of the volume element
	double V = m.ElementVolume(*pe);

    return eps*A/V;
}

//-----------------------------------------------------------------------------
void FESlidingInterfaceBiphasic::ProjectSurface(FESlidingSurfaceBiphasic& ss, FESlidingSurfaceBiphasic& ms, bool bupseg, bool bmove)
{
    FEMesh& mesh = GetFEModel()->GetMesh();
    double R = m_srad*mesh.GetBoundingBox().radius();
    
    // initialize projection data
    FENormalProjection np(ms);
    np.SetTolerance(m_stol);
    np.SetSearchRadius(R);
    np.Init();
    
    // if we need to project the nodes onto the master surface,
    // let's do this first
    if (bmove)
    {
        int NN = ss.Nodes();
        int NE = ss.Elements();
        // first we need to calculate the node normals
        vector<vec3d> normal; normal.assign(NN, vec3d(0,0,0));
        for (int i=0; i<NE; ++i)
        {
            FESurfaceElement& el = ss.Element(i);
            int ne = el.Nodes();
            for (int j=0; j<ne; ++j)
            {
                vec3d r0 = ss.Node(el.m_lnode[ j         ]).m_rt;
                vec3d rp = ss.Node(el.m_lnode[(j+   1)%ne]).m_rt;
                vec3d rm = ss.Node(el.m_lnode[(j+ne-1)%ne]).m_rt;
                vec3d n = (rp - r0)^(rm - r0);
                normal[el.m_lnode[j]] += n;
            }
        }
        for (int i=0; i<NN; ++i) normal[i].unit();
        
        // loop over all nodes
        for (int i=0; i<NN; ++i)
        {
            FENode& node = ss.Node(i);
            
            // get the spatial nodal coordinates
            vec3d rt = node.m_rt;
            vec3d nu = normal[i];
            
            // project onto the master surface
            vec3d q;
            double rs[2] = {0,0};
            FESurfaceElement* pme = np.Project(rt, nu, rs);
            if (pme)
            {
                // the node could potentially be in contact
                // find the global location of the intersection point
                vec3d q = ms.Local2Global(*pme, rs[0], rs[1]);
                
                // calculate the gap function
                // NOTE: this has the opposite sign compared
                // to Gerard's notes.
                double gap = nu*(rt - q);
                
                if (gap>0) node.m_r0 = node.m_rt = q;
            }
        }
    }
    
    // loop over all integration points
#pragma omp parallel for
    for (int i=0; i<ss.Elements(); ++i)
    {
        FESurfaceElement& el = ss.Element(i);
        bool sporo = ss.m_poro[i];
        
        int ne = el.Nodes();
        int nint = el.GaussPoints();
        double ps[FEElement::MAX_INTPOINTS];
        // get the nodal pressures
        if (sporo)
        {
            for (int j=0; j<ne; ++j) ps[j] = mesh.Node(el.m_node[j]).get(m_dofP);
        }
        
        for (int j=0; j<nint; ++j)
        {
            // get the integration point data
            FESlidingSurfaceBiphasic::Data& pt = ss.m_Data[i][j];
            
            // calculate the global position of the integration point
            vec3d r = ss.Local2Global(el, j);
            
            // get the pressure at the integration point
            double p1 = 0;
            if (sporo) p1 = el.eval(ps, j);
            
            // calculate the normal at this integration point
            vec3d nu = ss.SurfaceNormal(el, j);
            
            // first see if the old intersected face is still good enough
            FESurfaceElement* pme = pt.m_pme;
            double rs[2] = {0,0};
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
                
                double Ln = pt.m_Lmd + eps*g;
                
                pt.m_gap = (g <= R? g : 0);
                
                if ((Ln >= 0) && (g <= R))
                {
                    
                    // calculate the pressure gap function
                    bool mporo = ms.m_poro[pme->m_lid];
                    if (sporo) {
                        pt.m_p1 = p1;
                        if (mporo) {
                            double pm[FEElement::MAX_NODES];
                            for (int k=0; k<pme->Nodes(); ++k) pm[k] = mesh.Node(pme->m_node[k]).get(m_dofP);
                            double p2 = pme->eval(pm, rs[0], rs[1]);
                            pt.m_pg = p1 - p2;
                        }
                    }
                }
                else
                {
                    pt.m_Lmd = 0;
                    pt.m_gap = 0;
                    pt.m_dg = pt.m_Lmt = vec3d(0,0,0);
                    if (bmove) pt.m_pme = 0;
                    if (sporo) {
                        pt.m_Lmp = 0;
                        pt.m_pg = 0;
                        pt.m_p1 = 0;
                    }
                }
            }
            else
            {
                // the node is not in contact
                pt.m_Lmd = 0;
                pt.m_gap = 0;
                pt.m_dg = pt.m_Lmt = vec3d(0,0,0);
                if (sporo) {
                    pt.m_Lmp = 0;
                    pt.m_pg = 0;
                    pt.m_p1 = 0;
                }
            }
        }
    }
}

//-----------------------------------------------------------------------------

void FESlidingInterfaceBiphasic::Update(int nsolve_iter, const FETimeInfo& tp)
{
    double R = m_srad*GetFEModel()->GetMesh().GetBoundingBox().radius();
    
    static int naug = 0;
    static int biter = 0;
    
    FEModel& fem = *GetFEModel();
    
    // get the iteration number
    // we need this number to see if we can do segment updates or not
    // also reset number of iterations after each augmentation
    FEAnalysis* pstep = fem.GetCurrentStep();
    FESolver* psolver = pstep->GetFESolver();
    if (psolver->m_niter == 0) {
        biter = 0;
        naug = psolver->m_naug;
    } else if (psolver->m_naug > naug) {
        biter = psolver->m_niter;
        naug = psolver->m_naug;
    }
    int niter = psolver->m_niter - biter;
    bool bupseg = ((m_nsegup == 0)? true : (niter <= m_nsegup));
    // get the logfile
    //	Logfile& log = GetLogfile();
    //	log.printf("seg_up iteration # %d\n", niter+1);
    
    // project the surfaces onto each other
    // this will update the gap functions as well
    static bool bfirst = true;
    ProjectSurface(m_ss, m_ms, bupseg, (m_breloc && bfirst));
    if (m_btwo_pass || m_ms.m_bporo) ProjectSurface(m_ms, m_ss, bupseg);
    bfirst = false;
    
    if (nsolve_iter == 0)
    {
        m_ss.InitSlidingSurface();
        if (m_btwo_pass) m_ms.InitSlidingSurface();
    }
    
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
    for (int np=0; np<npass; ++np)
    {
        FESlidingSurfaceBiphasic& ss = (np == 0? m_ss : m_ms);
        FESlidingSurfaceBiphasic& ms = (np == 0? m_ms : m_ss);
        
        // loop over all the nodes of the primary surface
        for (int n=0; n<ss.Nodes(); ++n) {
            FENode& node = ss.Node(n);
            int id = node.m_ID[m_dofP];
            if ((id < -1) && (ss.m_tn[n].norm() > 0))
            {
                // mark node as non-free-draining (= pos ID)
                node.m_ID[m_dofP] = -id-2;
            }
        }
        
        // loop over all nodes of the secondary surface
        // the secondary surface is trickier since we need
        // to look at the primary surface's projection
        if (ms.m_bporo) {
            FENormalProjection np(ss);
            np.SetTolerance(m_stol);
            np.SetSearchRadius(R);
            np.Init();
            
            for (int n=0; n<ms.Nodes(); ++n)
            {
                // get the node
                FENode& node = ms.Node(n);
                
                // project it onto the primary surface
                double rs[2] = {0,0};
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
                        // get the normal tractions at the nodes
                        double tn[FEElement::MAX_NODES];
                        for (int i=0; i<pse->Nodes(); ++i)
                            tn[i] = ss.m_tn[pse->m_lnode[i]].norm();
                        
                        // now evaluate the traction at the intersection point
                        double tp = pse->eval(tn, rs[0], rs[1]);
                        
                        // if tp > 0, mark node as non-free-draining. (= pos ID)
                        int id = node.m_ID[m_dofP];
                        if ((id < -1) && (tp > 0))
                        {
                            // mark as non free-draining
                            node.m_ID[m_dofP] = -id-2;
                        }
                    }
                }
            }
        }
    }
}

//-----------------------------------------------------------------------------
vec3d FESlidingInterfaceBiphasic::SlipTangent(FESlidingSurfaceBiphasic& ss, const int nel, const int nint, FESlidingSurfaceBiphasic& ms, double& dh, vec3d& r)
{
    vec3d s1(0,0,0);
    dh = 0;
    
    // get primary surface element
    FESurfaceElement& se = ss.Element(nel);
    
    // get integration point data
    FESlidingSurfaceBiphasic::Data& data = ss.m_Data[nel][nint];
    double g = data.m_gap;
    vec3d nu = data.m_nu;
    
    // find secondary surface element
    FESurfaceElement* pme = data.m_pme;
    
    // calculate previous positions
    vec3d x2p = ms.Local2GlobalP(*pme, data.m_rs[0], data.m_rs[1]);
    vec3d x1p = ss.Local2GlobalP(se, nint);
    
    // calculate dx2
    vec3d x2 = ms.Local2Global(*pme, data.m_rs[0], data.m_rs[1]);
    vec3d dx2 = x2 - x2p;
    
    // calculate dx1
    vec3d x1 = ss.Local2Global(se, nint);
    vec3d dx1 = x1 - x1p;
    
    // get current and previous covariant basis vectors
    vec3d gscov[2], gscovp[2];
    ss.CoBaseVectors(se, nint, gscov);
    ss.CoBaseVectorsP(se, nint, gscovp);
    
    // calculate delta gscov
    vec3d dgscov[2];
    dgscov[0] = gscov[0] - gscovp[0];
    dgscov[1] = gscov[1] - gscovp[1];
    
    // calculate m, J, Nhat
    vec3d m = ((dgscov[0] ^ gscov[1]) + (gscov[0] ^ dgscov[1]));
    double detJ = (gscov[0] ^ gscov[1]).norm();
    mat3d Nhat = (mat3dd(1) - (nu & nu));
    
    // calculate q
    vec3d c = Nhat*m*(1.0/detJ);
    
    // calculate slip direction s1
    double norm = (Nhat*(c*(-g) + dx1 - dx2)).norm();
    if (norm != 0)
    {
        s1 = (Nhat*(c*(-g) + dx1 - dx2))/norm;
        dh = norm;
        r = c*(-g) + dx1 - dx2;
    }
    
    return s1;
    
}

//-----------------------------------------------------------------------------
vec3d FESlidingInterfaceBiphasic::ContactTraction(FESlidingSurfaceBiphasic& ss, const int nel, const int n, FESlidingSurfaceBiphasic& ms, double& pn)
{
    vec3d s1(0,0,0);
    vec3d drdot(0,0,0);
    vec3d t(0,0,0);
    pn = 0;
    double tn = 0, ts = 0, tns = 0;
    
    // get the integration point data
    FESlidingSurfaceBiphasic::Data& data = ss.m_Data[nel][n];
    
    // penalty
    double eps = m_epsn*data.m_epsn;
    
    // normal gap
    double g = data.m_gap;
    
    // normal traction Lagrange multiplier
    double Lm = data.m_Lmd;
    
    // vector traction Lagrange multiplier
    vec3d Lt = data.m_Lmt;
    
    // get the primary surface element
    FESurfaceElement& se = ss.Element(nel);
    
    // get the normal at this integration point
    vec3d nu = data.m_nu;
    
    // get the fluid pressure at this integration point
    double p = data.m_p1;
    
    // get current and previous secondary elements
    FESurfaceElement* pme = data.m_pme;
    FESurfaceElement* pmep = data.m_pmep;
    
    data.m_bstick = false;
    
    if (pme)
    {
        // assume stick and calculate traction
        if (pmep)
        {
            // calculate current global position of the integration point
            vec3d xo = ss.Local2Global(se, n);
            
            // calculate current global position of the previous intersection point
            vec3d xt = ms.Local2Global(*pmep, data.m_rsp[0], data.m_rsp[1]);
            
            // vector and normal gaps
            vec3d dg = xt - xo;
            
            // calculate trial stick traction, normal component, shear component
            t = Lt + dg*eps;
            tn = t*nu;
            ts = (t - nu*tn).norm();
            tns = tn + (1-m_phi)*p;
            
            // check if stick
            if ( (tn < 0) && (ts < m_mu*fabs(tns)) )
            {
                // set boolean flag for stick
                data.m_bstick = true;
                
                // contact pressure
                pn = MBRACKET(-tn);
                
                // store the previous values as the current
                data.m_pme = data.m_pmep;
                data.m_rs = data.m_rsp;
                
                // recalculate gap
                data.m_dg = dg;
            }
            else
            {
                // recalculate contact pressure for slip
                pn = MBRACKET(Lm + eps*g);
                
                if (pn != 0)
                {
                    
                    double dh = 0;
                    
                    // slip direction
                    s1 = SlipTangent(ss, nel, n, ms, dh, drdot);
                    
                    // total traction
                    tns = -pn + (1-m_phi)*p;
//                    t = nu*(-pn) + s1*(m_mu*tns);
                    t = nu*(-pn) - s1*(m_mu*fabs(tns));
                    
                    // reset slip direction
                    data.m_s1 = s1;
                }
                else
                {
                    t = vec3d(0,0,0);
                }
                
            }
        }
        else
        {
            // assume slip upon first contact
            // calculate contact pressure for slip
            pn = MBRACKET(Lm + eps*g);
            
            if (pn != 0)
            {
                
                double dh = 0;
                
                // TODO: do we only want to consider frictionless for first contact?
                
                // slip direction
                s1 = SlipTangent(ss, nel, n, ms, dh, drdot);
                
                // calculate frictional traction
                tns = -pn + (1-m_phi)*p;
//                t = nu*(-pn) + s1*(m_mu*tns);
                t = nu*(-pn) - s1*(m_mu*fabs(tns));
                
                // reset slip direction
                data.m_s1 = s1;
            }
        }
    }
    
    return t;
    
}

//-----------------------------------------------------------------------------
void FESlidingInterfaceBiphasic::Residual(FEGlobalVector& R, const FETimeInfo& tp)
{
    vector<int> sLM, mLM, LM, en;
    vector<double> fe;
    const int MN = FEElement::MAX_NODES;
    double detJ[MN], w[MN], *Hs, Hm[MN];
    double N[4*MN*2]; // TODO: is the size correct?
    
    m_ss.m_Ft = vec3d(0,0,0);
    m_ms.m_Ft = vec3d(0,0,0);
    
    FEModel& fem = *GetFEModel();
    
    // if we're using the symmetric formulation
    // we need to multiply with the timestep
    double dt = fem.GetCurrentStep()->m_dt;
    
	m_ss.m_Ft = vec3d(0, 0, 0);
	m_ms.m_Ft = vec3d(0, 0, 0);

    // loop over the nr of passes
    int npass = (m_btwo_pass?2:1);
    for (int np=0; np<npass; ++np)
    {
        // get slave and master surface
        FESlidingSurfaceBiphasic& ss = (np == 0? m_ss : m_ms);
        FESlidingSurfaceBiphasic& ms = (np == 0? m_ms : m_ss);
        
        // loop over all slave elements
        for (int i=0; i<ss.Elements(); ++i)
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
            for (int j=0; j<nint; ++j)
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
            for (int j=0; j<nint; ++j)
            {
                // get the integration point data
                FESlidingSurfaceBiphasic::Data& pt = ss.m_Data[i][j];
                
                // calculate contact pressure and account for stick
                double pn;
                vec3d t = ContactTraction(ss, i, j, ms, pn);
                
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
                    for (int k=0; k<nseln; ++k)
                    {
                        LM[3*k  ] = sLM[3*k  ];
                        LM[3*k+1] = sLM[3*k+1];
                        LM[3*k+2] = sLM[3*k+2];
                    }
                    
                    for (int k=0; k<nmeln; ++k)
                    {
                        LM[3*(k+nseln)  ] = mLM[3*k  ];
                        LM[3*(k+nseln)+1] = mLM[3*k+1];
                        LM[3*(k+nseln)+2] = mLM[3*k+2];
                    }
                    
                    // build the en vector
                    en.resize(nseln+nmeln);
                    for (int k=0; k<nseln; ++k) en[k      ] = se.m_node[k];
                    for (int k=0; k<nmeln; ++k) en[k+nseln] = me.m_node[k];
                    
                    // get slave element shape functions
                    Hs = se.H(j);
                    
                    // get master element shape functions
                    double r = pt.m_rs[0];
                    double s = pt.m_rs[1];
                    me.shape_fnc(Hm, r, s);
                    
                    if (pn > 0) {
                        
                        // calculate the force vector
                        fe.resize(ndof);
                        zero(fe);
                        
                        for (int k=0; k<nseln; ++k)
                        {
                            N[3*k  ] = Hs[k]*t.x;
                            N[3*k+1] = Hs[k]*t.y;
                            N[3*k+2] = Hs[k]*t.z;
                        }
                        
                        for (int k=0; k<nmeln; ++k)
                        {
                            N[3*(k+nseln)  ] = -Hm[k]*t.x;
                            N[3*(k+nseln)+1] = -Hm[k]*t.y;
                            N[3*(k+nseln)+2] = -Hm[k]*t.z;
                        }
                        
                        for (int k=0; k<ndof; ++k) fe[k] += N[k]*detJ[j]*w[j];
                        
                        // calculate contact forces
                        for (int k=0; k<nseln; ++k)
                            ss.m_Ft += vec3d(fe[3*k], fe[3*k+1], fe[3*k+2]);
                        
                        for (int k = 0; k<nmeln; ++k)
                            ms.m_Ft += vec3d(fe[3*(k+nseln)], fe[3*(k+nseln)+1], fe[3*(k+nseln)+2]);
                        
                        // assemble the global residual
                        R.Assemble(en, LM, fe);
                        
                        // do the biphasic stuff
                        if (sporo && mporo)
                        {
                            // calculate nr of pressure dofs
                            int ndof = nseln + nmeln;
                            
                            // calculate the flow rate
                            double epsp = m_epsp*pt.m_epsp;
                            
                            double wn = pt.m_Lmp + epsp*pt.m_pg;
                            
                            // fill the LM
                            LM.resize(ndof);
                            for (int k=0; k<nseln; ++k) LM[k        ] = sLM[3*nseln+k];
                            for (int k=0; k<nmeln; ++k) LM[k + nseln] = mLM[3*nmeln+k];
                            
                            // fill the force array
                            fe.resize(ndof);
                            zero(fe);
                            for (int k=0; k<nseln; ++k) N[k      ] =  Hs[k];
                            for (int k=0; k<nmeln; ++k) N[k+nseln] = -Hm[k];
                            
                            for (int k=0; k<ndof; ++k) fe[k] += dt*wn*N[k]*detJ[j]*w[j];
                            
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
void FESlidingInterfaceBiphasic::StiffnessMatrix(FESolver* psolver, const FETimeInfo& tp)
{
    vector<int> sLM, mLM, LM, en;
    const int MN = FEElement::MAX_NODES;
    double detJ[MN], w[MN], *Hs, Hm[MN], pt[MN], dpr[MN], dps[MN];
    double N[4*MN*2];
    matrix ke;
    
    FEModel& fem = *GetFEModel();
    
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
        FESlidingSurfaceBiphasic& ss = (np == 0? m_ss : m_ms);
        FESlidingSurfaceBiphasic& ms = (np == 0? m_ms : m_ss);
        
        FEMesh& mesh = *ms.GetMesh();
        
        // loop over all slave elements
        for (int i=0; i<ss.Elements(); ++i)
        {
            // get ths slave element
            FESurfaceElement& se = ss.Element(i);
            
            bool sporo = ss.m_poro[i];
            
            // get nr of nodes and integration points
            int nseln = se.Nodes();
            int nint = se.GaussPoints();
            
            // nodal pressures
            double pn[MN];
            for (int j=0; j<nseln; ++j) pn[j] = ss.GetMesh()->Node(se.m_node[j]).get(m_dofP);
            
            // copy the LM vector
            ss.UnpackLM(se, sLM);
            
            // we calculate all the metrics we need before we
            // calculate the nodal forces
            for (int j=0; j<nint; ++j)
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
            for (int j=0; j<nint; ++j)
            {
                // get integration point data
                FESlidingSurfaceBiphasic::Data& pt = ss.m_Data[i][j];
                
                // calculate contact pressure and account for stick
                double pn;
                vec3d t = ContactTraction(ss, i, j, ms, pn);
                
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
                    for (int k=0; k<nmeln; ++k) pm[k] = ms.GetMesh()->Node(me.m_node[k]).get(m_dofP);
                    
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
                        
                        for (int k=0; k<nseln; ++k)
                        {
                            LM[4*k  ] = sLM[3*k  ];			// x-dof
                            LM[4*k+1] = sLM[3*k+1];			// y-dof
                            LM[4*k+2] = sLM[3*k+2];			// z-dof
                            LM[4*k+3] = sLM[3*nseln+k];		// p-dof
                        }
                        for (int k=0; k<nmeln; ++k)
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
                        
                        for (int k=0; k<nseln; ++k)
                        {
                            LM[3*k  ] = sLM[3*k  ];
                            LM[3*k+1] = sLM[3*k+1];
                            LM[3*k+2] = sLM[3*k+2];
                        }
                        
                        for (int k=0; k<nmeln; ++k)
                        {
                            LM[3*(k+nseln)  ] = mLM[3*k  ];
                            LM[3*(k+nseln)+1] = mLM[3*k+1];
                            LM[3*(k+nseln)+2] = mLM[3*k+2];
                        }
                    }
                    
                    // build the en vector
                    en.resize(nseln+nmeln);
                    for (int k=0; k<nseln; ++k) en[k      ] = se.m_node[k];
                    for (int k=0; k<nmeln; ++k) en[k+nseln] = me.m_node[k];
                    
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
                    
                    // penalty
                    double eps = m_epsn*pt.m_epsn;
                    
                    // only evaluate stiffness matrix if contact traction is non-zero
                    if (pn > 0) {

                        // evaluate basis vectors on both surfaces
                        vec3d gscov[2], gmcov[2];
                        ss.CoBaseVectors(se, j, gscov);
                        ms.CoBaseVectors(me, r, s, gmcov);
                        
                        // identity tensor
                        mat3d I = mat3dd(1);
                        
                        // evaluate Mc and Ac and combine them into As
                        double* Gsr = se.Gr(j);
                        double* Gss = se.Gs(j);
                        mat3d Ac[MN], As[MN];
                        mat3d gscovh[2];
                        gscovh[0].skew(gscov[0]); gscovh[1].skew(gscov[1]);
                        for (int k=0; k<nseln; ++k) {
                            Ac[k] = (gscovh[1]*Gsr[k] - gscovh[0]*Gss[k])/detJ[j];
                            As[k] = t & (Ac[k]*nu);
                        }
                        
                        double Gmr[MN], Gms[MN];
                        me.shape_deriv(Gmr, Gms, r, s);
                        
                        mat2d A;
                        A[0][0] = gscov[0]*gmcov[0]; A[0][1] = gscov[0]*gmcov[1];
                        A[1][0] = gscov[1]*gmcov[0]; A[1][1] = gscov[1]*gmcov[1];
                        mat2d a = A.inverse();
                        
                        // evaluate covariant basis vectors on primary surface at previous time step
                        vec3d gscovp[2];
                        ss.CoBaseVectorsP(se, j, gscovp);
                        
                        // calculate delta gscov
                        vec3d dgscov[2];
                        dgscov[0] = gscov[0] - gscovp[0];
                        dgscov[1] = gscov[1] - gscovp[1];
                        
                        // evaluate approximate contravariant basis vectors when gap != 0
                        vec3d gscnt[2], gmcnt[2];
                        gmcnt[0] = gscov[0]*a[0][0] + gscov[1]*a[0][1];
                        gmcnt[1] = gscov[0]*a[1][0] + gscov[1]*a[1][1];
                        gscnt[0] = gmcov[0]*a[0][0] + gmcov[1]*a[1][0];
                        gscnt[1] = gmcov[0]*a[0][1] + gmcov[1]*a[1][1];
                        
                        vec3d mb[MN];
                        for (int k=0; k<nmeln; ++k) {
                            mb[k] = gmcnt[0]*Gmr[k] + gmcnt[1]*Gms[k];
                        }
                        
                        // evaluate Gbc
                        matrix Gbc(nmeln,nseln);
                        for (int b=0; b<nmeln; ++b) {
                            for (int c=0; c<nseln; ++c) {
                                Gbc(b,c)
                                = (a[0][0]*Gmr[b]*Gsr[c]
                                   + a[0][1]*Gmr[b]*Gss[c]
                                   + a[1][0]*Gms[b]*Gsr[c]
                                   + a[1][1]*Gms[b]*Gss[c])*(-g);
                            }
                        }
                        
                        // if stick
                        if (pt.m_bstick)
                        {
                            double dtn = eps;
                            
                            // create the stiffness matrix
                            ke.resize(ndof, ndof); ke.zero();
                            
                            // --- S O L I D - S O L I D   C O N T A C T ---
                            
                            // a. I-term
                            //------------------------------------
                            
                            for (int k=0; k<nseln; ++k) N[k      ] =  Hs[k];
                            for (int k=0; k<nmeln; ++k) N[k+nseln] = -Hm[k];
                            
                            double tmp = knmult*dtn*detJ[j]*w[j];
                            for (int l=0; l<nseln+nmeln; ++l)
                            {
                                for (int k=0; k<nseln+nmeln; ++k)
                                {
                                    ke[k*ndpn  ][l*ndpn  ] -= -tmp*N[k]*N[l]*I[0][0];
                                    ke[k*ndpn  ][l*ndpn+1] -= -tmp*N[k]*N[l]*I[0][1];
                                    ke[k*ndpn  ][l*ndpn+2] -= -tmp*N[k]*N[l]*I[0][2];
                                    
                                    ke[k*ndpn+1][l*ndpn  ] -= -tmp*N[k]*N[l]*I[1][0];
                                    ke[k*ndpn+1][l*ndpn+1] -= -tmp*N[k]*N[l]*I[1][1];
                                    ke[k*ndpn+1][l*ndpn+2] -= -tmp*N[k]*N[l]*I[1][2];
                                    
                                    ke[k*ndpn+2][l*ndpn  ] -= -tmp*N[k]*N[l]*I[2][0];
                                    ke[k*ndpn+2][l*ndpn+1] -= -tmp*N[k]*N[l]*I[2][1];
                                    ke[k*ndpn+2][l*ndpn+2] -= -tmp*N[k]*N[l]*I[2][2];
                                }
                            }
                            
                            // b. A-term
                            //-------------------------------------
                            
                            tmp = knmult*detJ[j]*w[j];
                            // non-symmetric
                            for (int l=0; l<nseln; ++l)
                            {
                                for (int k=0; k<nseln+nmeln; ++k)
                                {
                                    ke[k*ndpn  ][l*ndpn  ] -= tmp*N[k]*As[l][0][0];
                                    ke[k*ndpn  ][l*ndpn+1] -= tmp*N[k]*As[l][0][1];
                                    ke[k*ndpn  ][l*ndpn+2] -= tmp*N[k]*As[l][0][2];
                                    
                                    ke[k*ndpn+1][l*ndpn  ] -= tmp*N[k]*As[l][1][0];
                                    ke[k*ndpn+1][l*ndpn+1] -= tmp*N[k]*As[l][1][1];
                                    ke[k*ndpn+1][l*ndpn+2] -= tmp*N[k]*As[l][1][2];
                                    
                                    ke[k*ndpn+2][l*ndpn  ] -= tmp*N[k]*As[l][2][0];
                                    ke[k*ndpn+2][l*ndpn+1] -= tmp*N[k]*As[l][2][1];
                                    ke[k*ndpn+2][l*ndpn+2] -= tmp*N[k]*As[l][2][2];
                                }
                            }
                        }
                        // if slip
                        else
                        {
                            double tn = -pn;
                            
                            // create the stiffness matrix
                            ke.resize(ndof, ndof); ke.zero();
                            
                            // obtain the slip direction s1 and inverse of spatial increment dh
                            double dh = 0, hd = 0;
                            vec3d dr(0,0,0);
                            vec3d s1 = SlipTangent(ss, i, j, ms, dh, dr);
                            
                            if (dh != 0)
                            {
                                hd = 1.0 / dh;
                            }
                            
                            // evaluate N and S tensors and approximations when gap != 0
                            mat3ds N1 = dyad(nu);
                            mat3d Nh1 = mat3dd(1) - (nu & nu);
                            mat3d Nb1 = mat3dd(1) - (gscov[0] & gscnt[0]) - (gscov[1] & gscnt[1]);
                            mat3d Nt1 = nu & (Nb1*nu);
                            mat3d S1 = s1 & nu;
                            mat3d Sh1 = (mat3dd(1) - (s1 & s1))*hd;
                            mat3d Sb1 = s1 & (Nb1*nu);
                            
                            // evaluate m, c, B, and R
                            // evaluate L1 from Mg and R
                            vec3d m = ((dgscov[0] ^ gscov[1]) + (gscov[0] ^ dgscov[1]));
                            vec3d c = Sh1*Nh1*m*(1/detJ[j]);
                            mat3d Mg = (mat3dd(1)*(nu * m) + (nu & m))*(1/detJ[j]);
                            mat3d B = (c & (Nb1*nu)) - Sh1*Nh1;
                            mat3d R = mat3dd(1)*(nu * dr) + (nu & dr);
                            mat3d L1 = Sh1*(Nh1*Mg - mat3dd(1) + R/(-g))*Nh1*(-g);
                            
                            // evaluate Mc and Ac and combine them into As
                            // evaluate s1 dyad (N1*mc - Ac*nu) + c dyad (N1*mc + Ac*nu)*g*hd as Pc
                            // evaluate Fc from Ac_bar (Ab)
                            double* Gsr = se.Gr(j);
                            double* Gss = se.Gs(j);
                            mat3d As[MN];
                            mat3d gscovh[2];
                            mat3d dgscovh[2];
                            mat3d Pc[MN];
                            mat3d Jc[MN];
                            gscovh[0].skew(gscov[0]); gscovh[1].skew(gscov[1]);
                            dgscovh[0].skew(dgscov[0]); dgscovh[1].skew(dgscov[1]);
                            for (int k=0; k<nseln; ++k) {
                                vec3d mc = gscnt[0]*Gsr[k] + gscnt[1]*Gss[k];
                                mat3d Mc = nu & mc;
                                mat3d Ac = (gscovh[1]*Gsr[k] - gscovh[0]*Gss[k])/detJ[j];
                                mat3d Ab = (dgscovh[1]*Gsr[k] - dgscovh[0]*Gss[k])/detJ[j];
                                Pc[k] = (s1 & (N1*mc - Ac*nu)) + ((c & (N1*mc + Ac*nu))*(-g));
                                As[k] = Ac + Mc*N1;
                                Jc[k] = (L1*Ac - (Sh1*Nh1*Ab*(-g)));
                            }
                            
                            // evaluate Mb
                            // evaluate s1 dyad mb and combine as Psb
                            mat3d Pb[MN];
                            for (int k=0; k<nmeln; ++k) {
                                mat3d Mb = (-nu) & mb[k] ;
                                Pb[k] = Mb - ((s1 & mb[k])*m_mu);
                            }
                            
                            // define T, Ttb
                            mat3d T = N1 + (S1*m_mu);
                            mat3d Ttb = Nt1 + (Sb1*m_mu);
                            
                            // --- S O L I D - S O L I D   C O N T A C T ---
                            
                            // a. NxN-term
                            //------------------------------------
                            
                            for (int k=0; k<nseln; ++k) N[k      ] =  Hs[k];
                            for (int k=0; k<nmeln; ++k) N[k+nseln] = -Hm[k];
                            
                            double tmp = knmult*detJ[j]*w[j];
                            for (int l=0; l<nseln+nmeln; ++l)
                            {
                                for (int k=0; k<nseln+nmeln; ++k)
                                {
                                    ke[k*ndpn  ][l*ndpn  ] -= -tmp*N[k]*N[l]*(eps*Ttb[0][0] + m_mu*tn*B[0][0]);
                                    ke[k*ndpn  ][l*ndpn+1] -= -tmp*N[k]*N[l]*(eps*Ttb[0][1] + m_mu*tn*B[0][1]);
                                    ke[k*ndpn  ][l*ndpn+2] -= -tmp*N[k]*N[l]*(eps*Ttb[0][2] + m_mu*tn*B[0][2]);
                                    
                                    ke[k*ndpn+1][l*ndpn  ] -= -tmp*N[k]*N[l]*(eps*Ttb[1][0] + m_mu*tn*B[1][0]);
                                    ke[k*ndpn+1][l*ndpn+1] -= -tmp*N[k]*N[l]*(eps*Ttb[1][1] + m_mu*tn*B[1][1]);
                                    ke[k*ndpn+1][l*ndpn+2] -= -tmp*N[k]*N[l]*(eps*Ttb[1][2] + m_mu*tn*B[1][2]);
                                    
                                    ke[k*ndpn+2][l*ndpn  ] -= -tmp*N[k]*N[l]*(eps*Ttb[2][0] + m_mu*tn*B[2][0]);
                                    ke[k*ndpn+2][l*ndpn+1] -= -tmp*N[k]*N[l]*(eps*Ttb[2][1] + m_mu*tn*B[2][1]);
                                    ke[k*ndpn+2][l*ndpn+2] -= -tmp*N[k]*N[l]*(eps*Ttb[2][2] + m_mu*tn*B[2][2]);
                                }
                            }
                            
                            // b. Na,Nb-term
                            //-------------------------------------
                            
                            tmp = knmult*tn*detJ[j]*w[j];
                            // non-symmetric
                            for (int l=0; l<nseln; ++l)
                            {
                                for (int k=0; k<nseln+nmeln; ++k)
                                {
                                    ke[k*ndpn  ][l*ndpn  ] -= -tmp*N[k]*(As[l][0][0] + m_mu*(Pc[l][0][0] - Jc[l][0][0]));
                                    ke[k*ndpn  ][l*ndpn+1] -= -tmp*N[k]*(As[l][0][1] + m_mu*(Pc[l][0][1] - Jc[l][0][1]));
                                    ke[k*ndpn  ][l*ndpn+2] -= -tmp*N[k]*(As[l][0][2] + m_mu*(Pc[l][0][2] - Jc[l][0][2]));
                                    
                                    ke[k*ndpn+1][l*ndpn  ] -= -tmp*N[k]*(As[l][1][0] + m_mu*(Pc[l][1][0] - Jc[l][1][0]));
                                    ke[k*ndpn+1][l*ndpn+1] -= -tmp*N[k]*(As[l][1][1] + m_mu*(Pc[l][1][1] - Jc[l][1][1]));
                                    ke[k*ndpn+1][l*ndpn+2] -= -tmp*N[k]*(As[l][1][2] + m_mu*(Pc[l][1][2] - Jc[l][1][2]));
                                    
                                    ke[k*ndpn+2][l*ndpn  ] -= -tmp*N[k]*(As[l][2][0] + m_mu*(Pc[l][2][0] - Jc[l][2][0]));
                                    ke[k*ndpn+2][l*ndpn+1] -= -tmp*N[k]*(As[l][2][1] + m_mu*(Pc[l][2][1] - Jc[l][2][1]));
                                    ke[k*ndpn+2][l*ndpn+2] -= -tmp*N[k]*(As[l][2][2] + m_mu*(Pc[l][2][2] - Jc[l][2][2]));
                                }
                            }
                            
                            // c. Nc,Nd-term
                            //---------------------------------------
                            
                            tmp = tn*knmult*detJ[j]*w[j];
                            // non-symmetric
                            for (int k=0; k<nmeln; ++k)
                            {
                                for (int l=0; l<nseln+nmeln; ++l)
                                {
                                    ke[(k+nseln)*ndpn  ][l*ndpn  ] -= tmp*N[l]*Pb[k][0][0];
                                    ke[(k+nseln)*ndpn  ][l*ndpn+1] -= tmp*N[l]*Pb[k][0][1];
                                    ke[(k+nseln)*ndpn  ][l*ndpn+2] -= tmp*N[l]*Pb[k][0][2];
                                    
                                    ke[(k+nseln)*ndpn+1][l*ndpn  ] -= tmp*N[l]*Pb[k][1][0];
                                    ke[(k+nseln)*ndpn+1][l*ndpn+1] -= tmp*N[l]*Pb[k][1][1];
                                    ke[(k+nseln)*ndpn+1][l*ndpn+2] -= tmp*N[l]*Pb[k][1][2];
                                    
                                    ke[(k+nseln)*ndpn+2][l*ndpn  ] -= tmp*N[l]*Pb[k][2][0];
                                    ke[(k+nseln)*ndpn+2][l*ndpn+1] -= tmp*N[l]*Pb[k][2][1];
                                    ke[(k+nseln)*ndpn+2][l*ndpn+2] -= tmp*N[l]*Pb[k][2][2];
                                }
                            }
                            
                            // c. Gbc-term
                            //---------------------------------------
                            
                            tmp = tn*knmult*detJ[j]*w[j];
                            for (int k=0; k<nmeln; ++k)
                            {
                                for (int l=0; l<nseln; ++l)
                                {
                                    mat3d gT = T*(Gbc[k][l]*tmp);
                                    ke[(k+nseln)*ndpn  ][l*ndpn  ] -= gT[0][0];
                                    ke[(k+nseln)*ndpn  ][l*ndpn+1] -= gT[0][1];
                                    ke[(k+nseln)*ndpn  ][l*ndpn+2] -= gT[0][2];
                                    
                                    ke[(k+nseln)*ndpn+1][l*ndpn  ] -= gT[1][0];
                                    ke[(k+nseln)*ndpn+1][l*ndpn+1] -= gT[1][1];
                                    ke[(k+nseln)*ndpn+1][l*ndpn+2] -= gT[1][2];
                                    
                                    ke[(k+nseln)*ndpn+2][l*ndpn  ] -= gT[2][0];
                                    ke[(k+nseln)*ndpn+2][l*ndpn+1] -= gT[2][1];
                                    ke[(k+nseln)*ndpn+2][l*ndpn+2] -= gT[2][2];
                                }
                            }
                        }
                        
                        // --- B I P H A S I C   S T I F F N E S S ---
                        if (sporo && mporo)
                        {
                            // get the nodal pressures on the master surface
                            double pm[MN];
                            for (int k=0; k<nmeln; ++k)
                                pm[k] = mesh.Node(me.m_node[k]).get(m_dofP);
                            
                            // evaluate q2
                            double dpr=0, dps = 0;
                            for (int k=0; k<nmeln; ++k) {
                                dpr += pm[k]*Gmr[k];
                                dps += pm[k]*Gms[k];
                            }
                            vec3d q2 = gmcnt[0]*dpr + gmcnt[1]*dps;
                            
                            // evaluate gc
                            vector<double> gc(nseln);
                            for (int k=0; k<nseln; ++k) {
                                    gc[k]
                                    = (a[0][0]*dpr*Gsr[k]
                                       + a[0][1]*dpr*Gss[k]
                                       + a[1][0]*dps*Gsr[k]
                                       + a[1][1]*dps*Gss[k])*(-g);
                                }
                            
                            double dt = fem.GetCurrentStep()->m_dt;
                            double tmp = dt*w[j]*detJ[j];
                            
                            double epsp = m_epsp*pt.m_epsp;
                            
                            
                            // --- S O L I D - P R E S S U R E   C O N T A C T ---
                            
                            // a. q-term
                            //-------------------------------------
                            for (int k=0; k<nseln+nmeln; ++k)
                                for (int l=0; l<nseln+nmeln; ++l)
                                {
                                    ke[4*k + 3][4*l  ] -= -tmp*epsp*N[k]*N[l]*q2.x;
                                    ke[4*k + 3][4*l+1] -= -tmp*epsp*N[k]*N[l]*q2.y;
                                    ke[4*k + 3][4*l+2] -= -tmp*epsp*N[k]*N[l]*q2.z;
                                }
                            
                            double wn = pt.m_Lmp + epsp*pt.m_pg;
                            
                            // b. A-term
                            //-------------------------------------
                            
                            for (int l=0; l<nseln; ++l) {
                                vec3d Acn = Ac[l]*nu;
                                for (int k=0; k<nseln+nmeln; ++k)
                                {
                                    ke[4*k + 3][4*l  ] -= tmp*wn*N[k]*Acn.x;
                                    ke[4*k + 3][4*l+1] -= tmp*wn*N[k]*Acn.y;
                                    ke[4*k + 3][4*l+2] -= tmp*wn*N[k]*Acn.z;
                                }
                            }
                            
                            // c. m-term
                            //---------------------------------------
                            
                            for (int k=0; k<nmeln; ++k) {
                                for (int l=0; l<nseln+nmeln; ++l)
                                {
                                    ke[4*(k+nseln) + 3][4*l  ] -= -tmp*wn*N[l]*mb[k].x;
                                    ke[4*(k+nseln) + 3][4*l+1] -= -tmp*wn*N[l]*mb[k].y;
                                    ke[4*(k+nseln) + 3][4*l+2] -= -tmp*wn*N[l]*mb[k].z;
                                }
                            }
                            
                            // d. gc-term
                            //-------------------------------------
                            for (int k=0; k<nseln+nmeln; ++k)
                                for (int l=0; l<nseln; ++l)
                                {
                                    ke[4*k + 3][4*l  ] -= tmp*epsp*N[k]*gc[l]*nu.x;
                                    ke[4*k + 3][4*l+1] -= tmp*epsp*N[k]*gc[l]*nu.y;
                                    ke[4*k + 3][4*l+2] -= tmp*epsp*N[k]*gc[l]*nu.z;
                                }
                            
                            // e. Gbc-term
                            //---------------------------------------
                            
                            for (int k=0; k<nmeln; ++k) {
                                for (int l=0; l<nseln; ++l)
                                {
                                    ke[4*(k+nseln) + 3][4*l  ] -= tmp*wn*Gbc[k][l]*nu.x;
                                    ke[4*(k+nseln) + 3][4*l+1] -= tmp*wn*Gbc[k][l]*nu.y;
                                    ke[4*(k+nseln) + 3][4*l+2] -= tmp*wn*Gbc[k][l]*nu.z;
                                }
                            }
                            
                            // --- P R E S S U R E - P R E S S U R E   C O N T A C T ---
                            
                            // calculate the N-vector
                            for (int k=0; k<nseln; ++k)
                            {
                                N[ndpn*k  ] = 0;
                                N[ndpn*k+1] = 0;
                                N[ndpn*k+2] = 0;
                                N[ndpn*k+3] = Hs[k];
                            }
                            
                            for (int k=0; k<nmeln; ++k)
                            {
                                N[ndpn*(k+nseln)  ] = 0;
                                N[ndpn*(k+nseln)+1] = 0;
                                N[ndpn*(k+nseln)+2] = 0;
                                N[ndpn*(k+nseln)+3] = -Hm[k];
                            }
                            
                            for (int k=0; k<ndof; ++k)
                                for (int l=0; l<ndof; ++l) ke[k][l] -= tmp*epsp*N[k]*N[l];
                            
                        }
                        
                        // assemble the global stiffness
                        psolver->AssembleStiffness(en, LM, ke);
                    }
                }
            }
        }
    }
}

//-----------------------------------------------------------------------------
void FESlidingInterfaceBiphasic::UpdateContactPressures()
{
    int npass = (m_btwo_pass?2:1);
    const int MN = FEElement::MAX_NODES;
    const int MI = FEElement::MAX_INTPOINTS;
    for (int np=0; np<npass; ++np)
    {
        FESlidingSurfaceBiphasic& ss = (np == 0? m_ss : m_ms);
        FESlidingSurfaceBiphasic& ms = (np == 0? m_ms : m_ss);
        
        // loop over all elements of the primary surface
#pragma omp parallel for
        for (int n=0; n<ss.Elements(); ++n)
        {
            FESurfaceElement& el = ss.Element(n);
            int nint = el.GaussPoints();
            
            // get the normal tractions at the integration points
            for (int i=0; i<nint; ++i)
            {
                // get integration point data
                FESlidingSurfaceBiphasic::Data& sd = ss.m_Data[n][i];
                
                // evaluate traction on slave surface
                double eps = m_epsn*sd.m_epsn;
                if (sd.m_bstick) {
                    // if stick, evaluate total traction
                    sd.m_tr = sd.m_Lmt + sd.m_dg*eps;
                    // then derive normal component
                    sd.m_Ln = -sd.m_tr*sd.m_nu;
                }
                else {
                    // if slip, evaluate normal traction
                    double Ln = sd.m_Lmd + eps*sd.m_gap;
                    sd.m_Ln = MBRACKET(Ln);
                    // then derive total traction
                    sd.m_tr = -(sd.m_nu + sd.m_s1*m_mu)*sd.m_Ln;
                }
                
                FESurfaceElement* pme = sd.m_pme;
                
                if (m_btwo_pass && pme)
                {
                    // get master element data
                    vector<FESlidingSurfaceBiphasic::Data>& mdv = ms.m_Data[pme->m_lid];
                    int mint = pme->GaussPoints();
                    double pi[MI];
                    vec3d ti[MI];
                    for (int j=0; j<mint; ++j)
                    {
                        FESlidingSurfaceBiphasic::Data& md = mdv[j];
                        // evaluate traction on master surface
                        double eps = m_epsn*md.m_epsn;
                        if (md.m_bstick) {
                            // if stick, evaluate total traction
                            ti[j] = md.m_Lmt + md.m_dg*eps;
                            // then derive normal component
                            pi[j] = -ti[j]*md.m_nu;
                        }
                        else {
                            // if slip, evaluate normal traction
                            double Ln = md.m_Lmd + eps*md.m_gap;
                            pi[j] = MBRACKET(Ln);
                            // then derive total traction
                            ti[j] = -(md.m_nu + md.m_s1*m_mu)*pi[j];
                        }
                    }
                    // project the data to the nodes
                    double pn[MN];
                    vec3d tn[MN];
                    pme->project_to_nodes(pi, pn);
                    pme->project_to_nodes(ti, tn);
                    // now evaluate the traction at the intersection point
                    double Ln = pme->eval(pn, sd.m_rs[0], sd.m_rs[1]);
                    vec3d trac = pme->eval(tn, sd.m_rs[0], sd.m_rs[1]);
                    sd.m_Ln += MBRACKET(Ln);
                    // tractions on master-slave are opposite, so subtract
                    sd.m_tr -= trac;
                }
            }
        }
        ss.EvaluateNodalContactTractions();
    }
}

//-----------------------------------------------------------------------------
bool FESlidingInterfaceBiphasic::Augment(int naug, const FETimeInfo& tp)
{
    // make sure we need to augment
    if (!m_blaugon) return true;
    
    int i;
    double Ln, Lp;
    bool bconv = true;
    
    bool bporo = (m_ss.m_bporo && m_ms.m_bporo);
    int NS = (int)m_ss.m_Data.size();
    int NM = (int)m_ms.m_Data.size();
    
    // --- c a l c u l a t e   i n i t i a l   n o r m s ---
    // a. normal component
    double normL0 = 0, normP = 0, normDP = 0;
    for (int i=0; i<NS; ++i)
    {
        vector<FESlidingSurfaceBiphasic::Data>& sd = m_ss.m_Data[i];
        for (int j=0; j<(int)sd.size(); ++j)
        {
            FESlidingSurfaceBiphasic::Data& ds = sd[j];
            if (ds.m_bstick)
                normL0 += ds.m_Lmt*ds.m_Lmt;
            else
                normL0 += ds.m_Lmd*ds.m_Lmd;
        }
    }
    for (int i=0; i<NM; ++i)
    {
        vector<FESlidingSurfaceBiphasic::Data>& md = m_ms.m_Data[i];
        for (int j=0; j<(int)md.size(); ++j)
        {
            FESlidingSurfaceBiphasic::Data& dm = md[j];
            if (dm.m_bstick)
                normL0 += dm.m_Lmt*dm.m_Lmt;
            else
                normL0 += dm.m_Lmd*dm.m_Lmd;
        }
    }
    
    // b. gap component
    // (is calculated during update)
    double maxgap = 0;
    double maxpg = 0;
    
    // update Lagrange multipliers
    double normL1 = 0, epsp;
    for (i=0; i<NS; ++i)
    {
        vector<FESlidingSurfaceBiphasic::Data>& sd = m_ss.m_Data[i];
        for (int j=0; j<(int)sd.size(); ++j)
        {
            FESlidingSurfaceBiphasic::Data& ds = sd[j];
            
            // update Lagrange multipliers on slave surface
            double eps = m_epsn*ds.m_epsn;
            if (ds.m_bstick) {
                // if stick, augment total traction
                ds.m_Lmt += ds.m_dg*eps;
                // then derive normal component
                ds.m_Lmd = -ds.m_Lmt*ds.m_nu;
                normL1 += ds.m_Lmt*ds.m_Lmt;
            }
            else {
                // if slip, augment normal traction
                Ln = ds.m_Lmd + eps*ds.m_gap;
                ds.m_Lmd = MBRACKET(Ln);
                // then derive total traction
                ds.m_Lmt = -(ds.m_nu + ds.m_s1*m_mu)*ds.m_Lmd;
                normL1 += ds.m_Lmd*ds.m_Lmd;
            }
            
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
        vector<FESlidingSurfaceBiphasic::Data>& md = m_ms.m_Data[i];
        for (int j=0; j<(int)md.size(); ++j)
        {
            FESlidingSurfaceBiphasic::Data& dm = md[j];
            
            // update Lagrange multipliers on master surface
            double eps = m_epsn*dm.m_epsn;
            if (dm.m_bstick) {
                // if stick, augment total traction
                dm.m_Lmt += dm.m_dg*eps;
                // then derive normal component
                dm.m_Lmd = -dm.m_Lmt*dm.m_nu;
                normL1 += dm.m_Lmt*dm.m_Lmt;
            }
            else {
                // if slip, augment normal traction
                Ln = dm.m_Lmd + eps*dm.m_gap;
                dm.m_Lmd = MBRACKET(Ln);
                // then derive total traction
                dm.m_Lmt = -(dm.m_nu + dm.m_s1*m_mu)*dm.m_Lmd;
                normL1 += dm.m_Lmd*dm.m_Lmd;
            }
            
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
    
    // normP should be a measure of the fluid pressure at the
    // contact interface.  However, since it could be zero,
    // use an average measure of the contact traction instead.
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
    
    felog.printf(" sliding interface # %d\n", GetID());
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
void FESlidingInterfaceBiphasic::Serialize(DumpStream &ar)
{
    // serialize contact data
    FEContactInterface::Serialize(ar);
    
    // serialize contact surface data
    m_ms.Serialize(ar);
    m_ss.Serialize(ar);
}

//-----------------------------------------------------------------------------

void FESlidingInterfaceBiphasic::MarkFreeDraining()
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
        FESlidingSurfaceBiphasic& s = (np == 0? m_ss : m_ms);
        
        if (s.m_bporo) {
            // first, mark all nodes as free-draining (= neg. ID)
            // this is done by setting the dof's equation number
            // to a negative number
            for (i=0; i<s.Nodes(); ++i) 
            {
                id = s.Node(i).m_ID[m_dofP];
                if (id >= 0) 
                {
                    FENode& node = s.Node(i);
                    // mark node as free-draining
                    node.m_ID[m_dofP] = -id-2;
                }
            }
        }
    }
}

//-----------------------------------------------------------------------------
void FESlidingInterfaceBiphasic::SetFreeDraining()
{	
    int i, np;
    
    // Set the pressure to zero for the free-draining nodes
    for (np=0; np<2; ++np)
    {
        FESlidingSurfaceBiphasic& s = (np == 0? m_ss : m_ms);
        
        if (s.m_bporo) {
            // loop over all nodes
            for (i=0; i<s.Nodes(); ++i) 
            {
                if (s.Node(i).m_ID[m_dofP] < -1)
                {
                    FENode& node = s.Node(i);
                    // set the fluid pressure to zero
                    node.set(m_dofP, 0);
                }
            }
        }
    }
}
