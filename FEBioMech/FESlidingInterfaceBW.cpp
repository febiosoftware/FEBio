#include "stdafx.h"
#include "FESlidingInterfaceBW.h"
#include "FECore/FENormalProjection.h"
#include "FECore/FEModel.h"
#include "FECore/FEAnalysis.h"
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
ADD_PARAMETER(m_breloc   , FE_PARAM_BOOL  , "node_reloc"         );
ADD_PARAMETER(m_mu       , FE_PARAM_DOUBLE, "fric_coeff"         );
ADD_PARAMETER(m_bsmaug   , FE_PARAM_BOOL  , "smooth_aug"         );
ADD_PARAMETER(m_bflipm   , FE_PARAM_BOOL  , "flip_master"        );
ADD_PARAMETER(m_bflips   , FE_PARAM_BOOL  , "flip_slave"         );
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FESlidingSurfaceBW::Data::Data()
{
    m_gap = 0.0;
    m_Lmd = 0.0;
    m_Ln = 0.0;
    m_epsn = 1.0;
    m_Lmt = vec3d(0,0,0);
    m_nu = m_s1 = m_dg = vec3d(0,0,0);
    m_tr = vec3d(0,0,0);
    m_rs = m_rsp = vec2d(0,0);
    m_bstick = false;
    
    m_pme = m_pmep = (FESurfaceElement*)0;
}

//-----------------------------------------------------------------------------
// FESlidingSurfaceBW
//-----------------------------------------------------------------------------

FESlidingSurfaceBW::FESlidingSurfaceBW(FEModel* pfem) : FEContactSurface(pfem)
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
void FESlidingSurfaceBW::InitSlidingSurface()
{
    for (int i=0; i<Elements(); ++i)
    {
        FESurfaceElement& el = Element(i);
        int nint = el.GaussPoints();
        for (int j=0; j<nint; ++j)
        {
            // Store current surface projection values as previous
            m_Data[i][j].m_rsp = m_Data[i][j].m_rs;
            m_Data[i][j].m_pmep = m_Data[i][j].m_pme;
        }
    }
}

//-----------------------------------------------------------------------------
void FESlidingSurfaceBW::Serialize(DumpStream& ar)
{
    // We can serialize the base-class data
    FEContactSurface::Serialize(ar);
    
    if (ar.IsShallow())
    {
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
                    ar << d.m_Lmt;
                    ar << d.m_Lmd;
                    ar << d.m_epsn;
                    ar << d.m_Ln;
                    ar << d.m_bstick;
                    ar << d.m_tr;
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
                    ar >> d.m_dg;
                    ar >> d.m_nu;
                    ar >> d.m_s1;
                    ar >> d.m_rs;
                    ar >> d.m_rsp;
                    ar >> d.m_Lmt;
                    ar >> d.m_Lmd;
                    ar >> d.m_epsn;
                    ar >> d.m_Ln;
                    ar >> d.m_bstick;
                    ar >> d.m_tr;
                }
            }
        }
    }
    else
    {
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
                    ar << d.m_Lmt;
                    ar << d.m_Lmd;
                    ar << d.m_epsn;
                    ar << d.m_Ln;
                    ar << d.m_bstick;
                    ar << d.m_tr;
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
                    ar >> d.m_dg;
                    ar >> d.m_nu;
                    ar >> d.m_s1;
                    ar >> d.m_rs;
                    ar >> d.m_rsp;
                    ar >> d.m_Lmt;
                    ar >> d.m_Lmd;
                    ar >> d.m_epsn;
                    ar >> d.m_Ln;
                    ar >> d.m_bstick;
                    ar >> d.m_tr;
                }
            }
        }
    }
}

//-----------------------------------------------------------------------------
vec3d FESlidingSurfaceBW::GetContactForce()
{
    return m_Ft;
}

//-----------------------------------------------------------------------------
double FESlidingSurfaceBW::GetContactArea()
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
void FESlidingSurfaceBW::GetContactGap(int nface, double& pg)
{
    FESurfaceElement& el = Element(nface);
    int ni = el.GaussPoints();
    pg = 0;
    for (int k=0; k<ni; ++k) pg += m_Data[nface][k].m_gap;
    pg /= ni;
}

//-----------------------------------------------------------------------------
void FESlidingSurfaceBW::GetVectorGap(int nface, vec3d& pg)
{
    FESurfaceElement& el = Element(nface);
    int ni = el.GaussPoints();
    pg = vec3d(0,0,0);
    for (int k=0; k<ni; ++k) pg += m_Data[nface][k].m_dg;
    pg /= ni;
}

//-----------------------------------------------------------------------------
void FESlidingSurfaceBW::GetContactPressure(int nface, double& pg)
{
    FESurfaceElement& el = Element(nface);
    int ni = el.GaussPoints();
    pg = 0;
    for (int k=0; k<ni; ++k) pg += m_Data[nface][k].m_Ln;
    pg /= ni;
}

//-----------------------------------------------------------------------------
void FESlidingSurfaceBW::GetContactTraction(int nface, vec3d& pt)
{
    FESurfaceElement& el = Element(nface);
    int ni = el.GaussPoints();
    pt = vec3d(0,0,0);
    for (int k=0; k<ni; ++k) pt += m_Data[nface][k].m_tr;
    pt /= ni;
}

//-----------------------------------------------------------------------------
void FESlidingSurfaceBW::GetNodalContactGap(int nface, double* pg)
{
    FESurfaceElement& el = Element(nface);
    int ni = el.GaussPoints();
    double gi[FEElement::MAX_INTPOINTS];
    for (int k=0; k<ni; ++k) gi[k] = m_Data[nface][k].m_gap;
    el.project_to_nodes(gi, pg);
}

//-----------------------------------------------------------------------------
void FESlidingSurfaceBW::GetNodalVectorGap(int nface, vec3d* pg)
{
    FESurfaceElement& el = Element(nface);
    int ni = el.GaussPoints();
    vec3d gi[FEElement::MAX_INTPOINTS];
    for (int k=0; k<ni; ++k) gi[k] = m_Data[nface][k].m_dg;
    el.project_to_nodes(gi, pg);
}

//-----------------------------------------------------------------------------
void FESlidingSurfaceBW::GetNodalContactPressure(int nface, double* pg)
{
    FESurfaceElement& el = Element(nface);
    int ni = el.GaussPoints();
    double ti[FEElement::MAX_INTPOINTS];
    for (int k=0; k<ni; ++k) ti[k] = m_Data[nface][k].m_Ln;
    el.project_to_nodes(ti, pg);
}

//-----------------------------------------------------------------------------
void FESlidingSurfaceBW::GetNodalContactTraction(int nface, vec3d* pt)
{
    FESurfaceElement& el = Element(nface);
    int ni = el.GaussPoints();
    vec3d ti[FEElement::MAX_INTPOINTS];
    for (int k=0; k<ni; ++k) ti[k] = m_Data[nface][k].m_tr;
    el.project_to_nodes(ti, pt);
}

//-----------------------------------------------------------------------------
void FESlidingSurfaceBW::GetStickStatus(int nface, double& pg)
{
    FESurfaceElement& el = Element(nface);
    int ni = el.GaussPoints();
    pg = 0;
    for (int k=0; k<ni; ++k)
        if (m_Data[nface][k].m_bstick) pg += 1.0;
    pg /= ni;
}

//-----------------------------------------------------------------------------
// FESlidingInterfaceBW
//-----------------------------------------------------------------------------

FESlidingInterfaceBW::FESlidingInterfaceBW(FEModel* pfem) : FEContactInterface(pfem), m_ss(pfem), m_ms(pfem)
{
    static int count = 1;
    SetID(count++);
    
    // initial values
    m_knmult = 0;
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
    m_breloc = false;
    m_bsmaug = false;
    m_mu = 0.0;
    
    m_naugmin = 0;
    m_naugmax = 10;
    
    m_bfreeze = false;
    m_bflipm = m_bflips = false;

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
    
    
    // check friction and tension parameters
    // since they cannot be used simultaneously
    if ((m_mu != 0) && m_btension)
        return fecore_error("The tension option cannot be used with friction in sliding-elastic.");
    
    // when the friction coefficient is non-zero include higher-order terms in stiffness matrix
    if (m_mu != 0) m_knmult = 1;
    
    // initialize surface data
    if (m_ss.Init() == false) return false;
    if (m_ms.Init() == false) return false;
    
    if (m_bflips) m_ss.Invert();
    if (m_bflipm) m_ms.Invert();

    return true;
}

//-----------------------------------------------------------------------------
void FESlidingInterfaceBW::Activate()
{
    // don't forget to call the base class!
    FEContactInterface::Activate();
    
    // calculate the penalty
    if (m_bautopen)
    {
        CalcAutoPenalty(m_ss);
        CalcAutoPenalty(m_ms);
    }
    
    // update sliding interface data
    Update(0, GetFEModel()->GetTime());
}

//-----------------------------------------------------------------------------
//! build the matrix profile for use in the stiffness matrix
void FESlidingInterfaceBW::BuildMatrixProfile(FEGlobalMatrix& K)
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
    
    vector<int> lm(6*FEElement::MAX_NODES*2);
    
    int npass = (m_btwo_pass?2:1);
    for (int np=0; np<npass; ++np)
    {
        FESlidingSurfaceBW& ss = (np == 0? m_ss : m_ms);
        
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
                    FESurfaceElement& me = *pe;
                    int* mn = &me.m_node[0];
                    
                    assign(lm, -1);
                    
                    int nseln = se.Nodes();
                    int nmeln = me.Nodes();
                    
                    for (l=0; l<nseln; ++l)
                    {
                        vector<int>& id = mesh.Node(sn[l]).m_ID;
                        lm[6*l  ] = id[dof_X];
                        lm[6*l+1] = id[dof_Y];
                        lm[6*l+2] = id[dof_Z];
                        lm[6*l+3] = id[dof_RU];
                        lm[6*l+4] = id[dof_RV];
                        lm[6*l+5] = id[dof_RW];
                    }
                    
                    for (l=0; l<nmeln; ++l)
                    {
                        vector<int>& id = mesh.Node(mn[l]).m_ID;
                        lm[6*(l+nseln)  ] = id[dof_X];
                        lm[6*(l+nseln)+1] = id[dof_Y];
                        lm[6*(l+nseln)+2] = id[dof_Z];
                        lm[6*(l+nseln)+3] = id[dof_RU];
                        lm[6*(l+nseln)+4] = id[dof_RV];
                        lm[6*(l+nseln)+5] = id[dof_RW];
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
        
        // calculate a penalty
        double eps = AutoPenalty(el, s);
        
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
void FESlidingInterfaceBW::ProjectSurface(FESlidingSurfaceBW& ss, FESlidingSurfaceBW& ms, bool bupseg, bool bmove)
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
            if (pme == 0 && bupseg) pme = np.Project(r, nu, rs);
            
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
                    data.m_Lmd = 0;
                    data.m_pme = 0;
                    data.m_gap = 0;
                    data.m_dg = data.m_Lmt = vec3d(0,0,0);
                }
            }
            else
            {
                // the node is not in contact
                data.m_Lmd = 0;
                data.m_gap = 0;
                data.m_dg = data.m_Lmt = vec3d(0,0,0);
            }
        }
    }
}

//-----------------------------------------------------------------------------

void FESlidingInterfaceBW::Update(int nsolve_iter, const FETimeInfo& tp)
{
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
    bfirst = false;
    if (m_btwo_pass) ProjectSurface(m_ms, m_ss, bupseg);
    
    if (nsolve_iter == 0)
    {
        m_ss.InitSlidingSurface();
        if (m_btwo_pass) m_ms.InitSlidingSurface();
        m_bfreeze = false;
    }
    
    // Update the net contact pressures
    UpdateContactPressures();
    
    if (niter == 0) m_bfreeze = false;
    
    return;
}

//-----------------------------------------------------------------------------
vec3d FESlidingInterfaceBW::SlipTangent(FESlidingSurfaceBW& ss, const int nel, const int nint, FESlidingSurfaceBW& ms, double& dh, vec3d& r)
{
    vec3d s1(0,0,0);
    dh = 0;
    
    // get primary surface element
    FESurfaceElement& se = ss.Element(nel);
    
    // get integration point data
    FESlidingSurfaceBW::Data& data = ss.m_Data[nel][nint];
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
    double norm = (Nhat*(c*(-g)*m_knmult + dx1 - dx2)).norm();
    if (norm != 0)
    {
        s1 = (Nhat*(c*(-g)*m_knmult + dx1 - dx2))/norm;
        dh = norm;
        r = c*(-g)*m_knmult + dx1 - dx2;
    }
    
    return s1;
    
}

//-----------------------------------------------------------------------------
vec3d FESlidingInterfaceBW::ContactTraction(FESlidingSurfaceBW& ss, const int nel, const int n, FESlidingSurfaceBW& ms, double& pn)
{
    vec3d s1(0,0,0);
    vec3d drdot(0,0,0);
    vec3d t(0,0,0);
    pn = 0;
    double tn = 0, ts = 0;
    
    // get the integration point data
    FESlidingSurfaceBW::Data& data = ss.m_Data[nel][n];
    
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
    
    // calculate the normal at this integration point
    vec3d nu = data.m_nu;
    
    // get current and previous secondary elements
    FESurfaceElement* pme = data.m_pme;
    FESurfaceElement* pmep = data.m_pmep;
    
    // if we just returned from an augmentation, do not update stick or slip status
    if (m_bfreeze && pme) {
        if (data.m_bstick) {
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
            
            // contact pressure
            pn = m_btension ? (-tn) : MBRACKET(-tn);
            
            // store the previous values as the current
            data.m_pme = data.m_pmep;
            data.m_rs = data.m_rsp;
            
            // recalculate gap
            data.m_dg = dg;
        }
        else {
            // recalculate contact pressure for slip
            pn = m_btension ? (Lm + eps*g) : MBRACKET(Lm + eps*g);
            
            if (pn != 0)
            {
                
                double dh = 0;
                
                // slip direction
                s1 = FESlidingInterfaceBW::SlipTangent(ss, nel, n, ms, dh, drdot);
                
                // total traction
                t = (nu + s1*(m_mu))*(-pn);
                
                // reset slip direction
                data.m_s1 = s1;
            }
            else
            {
                t = vec3d(0,0,0);
            }
        }
    }
    // update contact tractions
    else {
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
                
                // check if stick
                if ( (tn < 0) && (ts < m_mu*fabs(tn)) )
                {
                    // set boolean flag for stick
                    data.m_bstick = true;
                    
                    // contact pressure
                    pn = m_btension ? (-tn) : MBRACKET(-tn);
                    
                    // store the previous values as the current
                    data.m_pme = data.m_pmep;
                    data.m_rs = data.m_rsp;
                    
                    // recalculate gap
                    data.m_dg = dg;
                }
                else
                {
                    // recalculate contact pressure for slip
                    pn = m_btension ? (Lm + eps*g) : MBRACKET(Lm + eps*g);
                    
                    if (pn != 0)
                    {
                        
                        double dh = 0;
                        
                        // slip direction
                        s1 = FESlidingInterfaceBW::SlipTangent(ss, nel, n, ms, dh, drdot);
                        
                        // total traction
                        t = (nu + s1*(m_mu))*(-pn);
                        
                        // reset slip direction
                        data.m_s1 = s1;
                        data.m_bstick = false;
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
                pn = m_btension ? (Lm + eps*g) : MBRACKET(Lm + eps*g);
                
                if (pn != 0)
                {
                    
                    double dh = 0;
                    
                    // slip direction
                    s1 = FESlidingInterfaceBW::SlipTangent(ss, nel, n, ms, dh, drdot);
                    
                    // calculate frictional traction
                    t = (nu + s1*(m_mu))*(-pn);
                    
                    // reset slip direction
                    data.m_s1 = s1;
                    data.m_bstick = false;
                }
            }
        }
    }
    
    return t;
    
}

//-----------------------------------------------------------------------------
void FESlidingInterfaceBW::Residual(FEGlobalVector& R, const FETimeInfo& tp)
{
    const int MN = FEElement::MAX_NODES;
    
    vector<int> sLM, mLM, LM, en;
    vector<double> fe;
    double detJ[MN], w[MN], Hm[MN];
    double N[MN*6];
    
    m_ss.m_Ft = vec3d(0,0,0);
    m_ms.m_Ft = vec3d(0,0,0);
    
    // loop over the nr of passes
    int npass = (m_btwo_pass?2:1);
    for (int np=0; np<npass; ++np)
    {
        // get slave and master surface
        FESlidingSurfaceBW& ss = (np == 0? m_ss : m_ms);
        FESlidingSurfaceBW& ms = (np == 0? m_ms : m_ss);
        
        // loop over all slave elements
        //#pragma omp parallel for private(sLM, mLM, LM, en, fe, detJ, w, Hm, N)
        for (int i=0; i<ss.Elements(); ++i)
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
                // get integration point data
                FESlidingSurfaceBW::Data& data = ss.m_Data[i][j];
                
                // calculate contact pressure and account for stick
                double pn;
                vec3d t = ContactTraction(ss, i, j, ms, pn);
                
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
                    double *Hs = se.H(j);
                    
                    // get master element shape functions
                    double r = data.m_rs[0];
                    double s = data.m_rs[1];
                    me.shape_fnc(Hm, r, s);
                    
                    if (pn != 0) {
                        
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
                        {
                            ss.m_Ft += vec3d(fe[k*3], fe[k*3+1], fe[k*3+2]);
                        }
                        
                        for (int k = 0; k<nmeln; ++k)
                        {
                            ms.m_Ft += vec3d(fe[(k + nseln) * 3], fe[(k + nseln) * 3 + 1], fe[(k + nseln) * 3 + 2]);
                        }
                        
                        // assemble the global residual
                        R.Assemble(en, LM, fe);
                    }
                }
            }
        }
    }
}

//-----------------------------------------------------------------------------
void FESlidingInterfaceBW::StiffnessMatrix(FESolver* psolver, const FETimeInfo& tp)
{
    // see how many reformations we've had to do so far
    int nref = psolver->m_nref;
    
    const int MN = FEElement::MAX_NODES;
    
    double detJ[MN], w[MN], Hm[MN];
    double N[MN*6];
    vector<int> sLM, mLM, LM, en;
    matrix ke;
    
    // do single- or two-pass
    int npass = (m_btwo_pass?2:1);
    for (int np=0; np < npass; ++np)
    {
        // get the slave and master surface
        FESlidingSurfaceBW& ss = (np == 0? m_ss : m_ms);
        FESlidingSurfaceBW& ms = (np == 0? m_ms : m_ss);
        
        // loop over all slave elements
        //#pragma omp parallel for private(detJ, w, Hm, N, sLM, mLM, LM, en, ke)
        for (int i=0; i<ss.Elements(); ++i)
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
            for (int j=0; j<nint; ++j)
            {
                // get integration point data
                FESlidingSurfaceBW::Data& data = ss.m_Data[i][j];
                
                // calculate contact pressure and account for stick
                double pn;
                vec3d t = ContactTraction(ss, i, j, ms, pn);
                
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
                    
                    // slave shape functions
                    double* Hs = se.H(j);
                    
                    // master shape functions
                    double r = data.m_rs[0];
                    double s = data.m_rs[1];
                    me.shape_fnc(Hm, r, s);
                    
                    // get slave normal vector
                    vec3d nu = data.m_nu;
                    
                    // gap function
                    double g = data.m_gap;
                    
                    // penalty
                    double eps = m_epsn*data.m_epsn;
                    
                    // only evaluate stiffness matrix if contact traction is non-zero
                    if (pn != 0)
                    {
                        // if stick
                        if (data.m_bstick)
                        {
                            double dtn = eps;
                            
                            // create the stiffness matrix
                            ke.resize(ndof, ndof); ke.zero();
                            
                            // evaluate basis vectors on primary surface
                            vec3d gscov[2];
                            ss.CoBaseVectors(se, j, gscov);
                            
                            // identity tensor
                            mat3d I = mat3dd(1);
                            
                            // evaluate Mc and Ac and combine them into As
                            double* Gsr = se.Gr(j);
                            double* Gss = se.Gs(j);
                            mat3d As[MN];
                            mat3d gscovh[2];
                            gscovh[0].skew(gscov[0]); gscovh[1].skew(gscov[1]);
                            for (int k=0; k<nseln; ++k) {
                                mat3d Ac = (gscovh[1]*Gsr[k] - gscovh[0]*Gss[k])/detJ[j];
                                As[k] = t & (Ac*nu);
                            }
                            
                            // --- S O L I D - S O L I D   C O N T A C T ---
                            
                            // a. I-term
                            //------------------------------------
                            
                            for (int k=0; k<nseln; ++k) N[k      ] =  Hs[k];
                            for (int k=0; k<nmeln; ++k) N[k+nseln] = -Hm[k];
                            
                            double tmp = dtn*detJ[j]*w[j];
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
                            
                            tmp = detJ[j]*w[j];
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
                            
                            // assemble the global stiffness
                            //					#pragma omp critical
                            {
                                psolver->AssembleStiffness(en, LM, ke);
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
                            vec3d s1 = FESlidingInterfaceBW::SlipTangent(ss, i, j, ms, dh, dr);
                            
                            if (dh != 0)
                            {
                                hd = 1.0 / dh;
                            }
                            
                            // evaluate basis vectors on both surfaces
                            vec3d gscov[2], gmcov[2];
                            ss.CoBaseVectors(se, j, gscov);
                            ms.CoBaseVectors(me, r, s, gmcov);
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
                            
                            // evaluate contravariant basis vectors
                            vec3d gscnt[2], gmcnt[2];
                            if (m_knmult == 0)
                            {
                                // evaluate true contravariant basis vectors when gap = 0
                                ss.ContraBaseVectors(se, j, gscnt);
                                ms.ContraBaseVectors(me, r, s, gmcnt);
                            }
                            else
                            {
                                // evaluate approximate contravariant basis vectors when gap != 0
                                gmcnt[0] = gscov[0]*a[0][0] + gscov[1]*a[0][1];
                                gmcnt[1] = gscov[0]*a[1][0] + gscov[1]*a[1][1];
                                gscnt[0] = gmcov[0]*a[0][0] + gmcov[1]*a[1][0];
                                gscnt[1] = gmcov[0]*a[0][1] + gmcov[1]*a[1][1];
                            }
                            
                            // evaluate N and S tensors and approximations when gap != 0
                            mat3ds N1 = dyad(nu);
                            mat3d Nh1 = mat3dd(1) - (nu & nu);
                            mat3d Nb1 = mat3dd(1) - (gscov[0] & gscnt[0])*m_knmult - (gscov[1] & gscnt[1])*m_knmult;
                            mat3d Nt1 = nu & (Nb1*nu);
                            mat3d S1 = s1 & nu;
                            mat3d Sh1 = (mat3dd(1) - (s1 & s1))*hd;
                            mat3d Sb1 = s1 & (Nb1*nu);
                            
                            // evaluate m, c, B, and R
                            // evaluate L1 from Mg and R
                            vec3d m = ((dgscov[0] ^ gscov[1]) + (gscov[0] ^ dgscov[1]));
                            vec3d c = Sh1*Nh1*m*(1/detJ[j]);
                            mat3d Mg = (mat3dd(1)*(nu * m) + (nu & m))*(1/detJ[j]);
                            mat3d B = (c & (Nb1*nu))*m_knmult - Sh1*Nh1;
                            mat3d R = mat3dd(1)*(nu * dr) + (nu & dr);
                            mat3d L1 = Sh1*((Nh1*Mg - mat3dd(1))*(-g)*m_knmult + R)*Nh1;
                            
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
                                Pc[k] = (s1 & (N1*mc*m_knmult - Ac*nu)) + ((c & (N1*mc + Ac*nu))*(-g)*m_knmult);
                                As[k] = Ac + Mc*N1*m_knmult;
                                Jc[k] = (L1*Ac - (Sh1*Nh1*Ab*(-g)*m_knmult));
                            }
                            
                            // evaluate Mb
                            // evaluate s1 dyad mb and combine as Psb
                            double Gmr[MN], Gms[MN];
                            me.shape_deriv(Gmr, Gms, r, s);
                            mat3d Pb[MN];
                            for (int k=0; k<nmeln; ++k) {
                                vec3d n(0,0,0);
                                if (m_knmult == 0)
                                {
                                    n = gmcnt[0] ^ gmcnt[1];
                                    n.unit();
                                } else
                                {
                                    n = -nu;
                                }
                                vec3d mb = gmcnt[0]*Gmr[k] + gmcnt[1]*Gms[k];
                                mat3d Mb = n & mb;
                                Pb[k] = Mb - ((s1 & mb)*m_mu);
                            }
                            
                            // evaluate Gbc
                            matrix Gbc(nmeln,nseln);
                            for (int b=0; b<nmeln; ++b) {
                                for (int c=0; c<nseln; ++c) {
                                    Gbc(b,c)
                                    = (a[0][0]*Gmr[b]*Gsr[c]
                                       + a[0][1]*Gmr[b]*Gss[c]
                                       + a[1][0]*Gms[b]*Gsr[c]
                                       + a[1][1]*Gms[b]*Gss[c])*(-g)*m_knmult;
                                }
                            }
                            
                            // define T, Ttb
                            mat3d T = N1 + (S1*m_mu);
                            mat3d Ttb = Nt1 + (Sb1*m_mu);
                            
                            // --- S O L I D - S O L I D   C O N T A C T ---
                            
                            // a. NxN-term
                            //------------------------------------
                            
                            for (int k=0; k<nseln; ++k) N[k      ] =  Hs[k];
                            for (int k=0; k<nmeln; ++k) N[k+nseln] = -Hm[k];
                            
                            double tmp = detJ[j]*w[j];
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
                            
                            tmp = tn*detJ[j]*w[j];
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
                            
                            tmp = tn*detJ[j]*w[j];
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
                            
                            tmp = tn*detJ[j]*w[j];
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
                            
                            // assemble the global stiffness
                            //					#pragma omp critical
                            {
                                psolver->AssembleStiffness(en, LM, ke);
                            }
                        }
                        
                    }
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
    const int MI = FEElement::MAX_INTPOINTS;
    for (int np=0; np<npass; ++np)
    {
        FESlidingSurfaceBW& ss = (np == 0? m_ss : m_ms);
        FESlidingSurfaceBW& ms = (np == 0? m_ms : m_ss);
        
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
                FESlidingSurfaceBW::Data& sd = ss.m_Data[n][i];
                
                double pn = 0;
                
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
                    sd.m_Ln = m_btension ? Ln : MBRACKET(Ln);
                    // then derive total traction
                    sd.m_tr = -(sd.m_nu + sd.m_s1*m_mu)*sd.m_Ln;
                }
                
                FESurfaceElement* pme = sd.m_pme;
                
                if (m_btwo_pass && pme)
                {
                    // get master element data
                    vector<FESlidingSurfaceBW::Data>& mdv = ms.m_Data[pme->m_lid];
                    int mint = pme->GaussPoints();
                    double pi[MI];
                    vec3d ti[MI];
                    for (int j=0; j<mint; ++j)
                    {
                        FESlidingSurfaceBW::Data& md = mdv[j];
                        
                        pn = 0;
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
                            pi[j] = m_btension ? Ln : MBRACKET(Ln);
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
                    sd.m_Ln += (m_btension ? Ln : MBRACKET(Ln));
                    // tractions on master-slave are opposite, so subtract
                    sd.m_tr -= trac;
                }
            }
        }
    }
}

//-----------------------------------------------------------------------------
bool FESlidingInterfaceBW::Augment(int naug, const FETimeInfo& tp)
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
            if (ds.m_bstick)
                normL0 += ds.m_Lmt*ds.m_Lmt;
            else
                normL0 += ds.m_Lmd*ds.m_Lmd;
        }
    }
    for (int i=0; i<NM; ++i)
    {
        vector<FESlidingSurfaceBW::Data>& md = m_ms.m_Data[i];
        for (int j=0; j<(int)md.size(); ++j)
        {
            FESlidingSurfaceBW::Data& dm = md[j];
            if (dm.m_bstick)
                normL0 += dm.m_Lmt*dm.m_Lmt;
            else
                normL0 += dm.m_Lmd*dm.m_Lmd;
            
        }
    }
    
    // b. gap component
    // (is calculated during update)
    double maxgap = 0;
    
    // update Lagrange multipliers
    double normL1 = 0;
    for (int i=0; i<m_ss.Elements(); ++i) {
        FESurfaceElement& el = m_ss.Element(i);
        vec3d tn[FEElement::MAX_INTPOINTS];
        if (m_bsmaug) m_ss.GetGPSurfaceTraction(i, tn);
        for (int j=0; j<el.GaussPoints(); ++j) {
            FESlidingSurfaceBW::Data& data = m_ss.m_Data[i][j];
            // update Lagrange multipliers on slave surface
            if (data.m_bstick) {
                // if stick, augment total traction
                if (m_bsmaug) {
                    data.m_Lmt = tn[j];
                    if (m_btwo_pass) data.m_Lmt /= 2;
                }
                else {
                    double eps = m_epsn*data.m_epsn;
                    data.m_Lmt += data.m_dg*eps;
                }
                // then derive normal component
                data.m_Lmd = -data.m_Lmt*data.m_nu;
                Ln = data.m_Lmd;
                normL1 += data.m_Lmt*data.m_Lmt;
                
                if (m_btension)
                    maxgap = max(maxgap,data.m_dg.norm());
                else if (Ln > 0) maxgap = max(maxgap,data.m_dg.norm());
            }
            else {
                // if slip, augment normal traction
                if (m_bsmaug) {
                    Ln = -(tn[j]*data.m_nu);
                    data.m_Lmd = m_btension ? Ln : MBRACKET(Ln);
                    if (m_btwo_pass) data.m_Lmd /= 2;
                }
                else {
                    double eps = m_epsn*data.m_epsn;
                    Ln = data.m_Lmd + eps*data.m_gap;
                    data.m_Lmd = m_btension ? Ln : MBRACKET(Ln);
                }
                // then derive total traction
                data.m_Lmt = -(data.m_nu + data.m_s1*m_mu)*data.m_Lmd;
                normL1 += data.m_Lmd*data.m_Lmd;
                
                if (m_btension)
                    maxgap = max(maxgap,fabs(data.m_gap));
                else if (Ln > 0) maxgap = max(maxgap,fabs(data.m_gap));
            }
        }
    }
    
    for (int i=0; i<m_ms.Elements(); ++i) {
        FESurfaceElement& el = m_ms.Element(i);
        vec3d tn[FEElement::MAX_INTPOINTS];
        if (m_bsmaug) m_ms.GetGPSurfaceTraction(i, tn);
        for (int j=0; j<el.GaussPoints(); ++j) {
            FESlidingSurfaceBW::Data& data = m_ms.m_Data[i][j];
            // update Lagrange multipliers on master surface
            if (data.m_bstick) {
                // if stick, augment total traction
                if (m_bsmaug) {
                    data.m_Lmt = tn[j];
                    if (m_btwo_pass) data.m_Lmt /= 2;
                }
                else {
                    double eps = m_epsn*data.m_epsn;
                    data.m_Lmt += data.m_dg*eps;
                }
                // then derive normal component
                data.m_Lmd = -data.m_Lmt*data.m_nu;
                Ln = data.m_Lmd;
                normL1 += data.m_Lmt*data.m_Lmt;
                
                if (m_btension)
                    maxgap = max(maxgap,fabs(data.m_dg.norm()));
                else if (Ln > 0) maxgap = max(maxgap,fabs(data.m_dg.norm()));
            }
            else {
                // if slip, augment normal traction
                if (m_bsmaug) {
                    Ln = -(tn[j]*data.m_nu);
                    data.m_Lmd = m_btension ? Ln : MBRACKET(Ln);
                    if (m_btwo_pass) data.m_Lmd /= 2;
                }
                else {
                    double eps = m_epsn*data.m_epsn;
                    Ln = data.m_Lmd + eps*data.m_gap;
                    data.m_Lmd = m_btension ? Ln : MBRACKET(Ln);
                }
                // then derive total traction
                data.m_Lmt = -(data.m_nu + data.m_s1*m_mu)*data.m_Lmd;
                normL1 += data.m_Lmd*data.m_Lmd;
                
                if (m_btension)
                    maxgap = max(maxgap,fabs(data.m_gap));
                else if (Ln > 0) maxgap = max(maxgap,fabs(data.m_gap));
            }
        }
    }
    
    // calculate relative norms
    double lnorm = (normL1 != 0 ? fabs((normL1 - normL0) / normL1) : fabs(normL1 - normL0));
    
    // check convergence
    if ((m_gtol > 0) && (maxgap > m_gtol)) bconv = false;
    if ((m_atol > 0) && (lnorm > m_atol)) bconv = false;
    
    if (naug < m_naugmin ) bconv = false;
    if (naug >= m_naugmax) bconv = true;
    
    felog.printf(" sliding interface # %d\n", GetID());
    felog.printf("                        CURRENT        REQUIRED\n");
    felog.printf("    D multiplier : %15le", lnorm); if (m_atol > 0) felog.printf("%15le\n", m_atol); else felog.printf("       ***\n");
    
    felog.printf("    maximum gap  : %15le", maxgap);
    if (m_gtol > 0) felog.printf("%15le\n", m_gtol); else felog.printf("       ***\n");
    
    ProjectSurface(m_ss, m_ms, true);
    if (m_btwo_pass) ProjectSurface(m_ms, m_ss, true);
    
    m_bfreeze = true;
    
    return bconv;
}

//-----------------------------------------------------------------------------
void FESlidingInterfaceBW::Serialize(DumpStream &ar)
{
    // serialize contact data
    FEContactInterface::Serialize(ar);
    
    // serialize contact surface data
    m_ms.Serialize(ar);
    m_ss.Serialize(ar);
}
