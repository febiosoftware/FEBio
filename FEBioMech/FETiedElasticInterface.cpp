/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/



#include "stdafx.h"
#include "FETiedElasticInterface.h"
#include "FECore/FEModel.h"
#include "FECore/FEAnalysis.h"
#include "FECore/FENormalProjection.h"
#include <FECore/FELinearSystem.h>
#include "FECore/log.h"

//-----------------------------------------------------------------------------
// Define sliding interface parameters
BEGIN_FECORE_CLASS(FETiedElasticInterface, FEContactInterface)
	ADD_PARAMETER(m_laugon   , "laugon"             )->setLongName("Enforcement method")->setEnums("PENALTY\0AUGLAG\0");
	ADD_PARAMETER(m_atol     , "tolerance"          );
	ADD_PARAMETER(m_gtol     , "gaptol"             );
	ADD_PARAMETER(m_epsn     , "penalty"            );
	ADD_PARAMETER(m_bautopen , "auto_penalty"       );
    ADD_PARAMETER(m_bupdtpen , "update_penalty"     );
	ADD_PARAMETER(m_btwo_pass, "two_pass"           );
	ADD_PARAMETER(m_knmult   , "knmult"             );
	ADD_PARAMETER(m_stol     , "search_tol"         );
	ADD_PARAMETER(m_bsymm    , "symmetric_stiffness");
	ADD_PARAMETER(m_srad     , "search_radius"      );
	ADD_PARAMETER(m_naugmin  , "minaug"             );
	ADD_PARAMETER(m_naugmax  , "maxaug"             );
	ADD_PARAMETER(m_bflips   , "flip_primary"       );
	ADD_PARAMETER(m_bflipm   , "flip_secondary"     );
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FETiedElasticSurface::Data::Data()
{
    m_Gap = vec3d(0,0,0);
    m_dg = vec3d(0,0,0);
    m_nu = vec3d(0,0,0);
    m_rs = vec2d(0,0);
    m_Lmd = vec3d(0,0,0);
    m_tr = vec3d(0,0,0);
    m_epsn = 1.0;
}

//-----------------------------------------------------------------------------
void FETiedElasticSurface::Data::Serialize(DumpStream& ar)
{
	FEContactMaterialPoint::Serialize(ar);
	ar & m_Gap;
	ar & m_dg;
	ar & m_nu;
	ar & m_rs;
	ar & m_Lmd;
	ar & m_epsn;
	ar & m_tr;
}

void FETiedElasticSurface::Data::Init()
{
    FEContactMaterialPoint::Init();
    m_Gap = vec3d(0, 0, 0);
    m_dg = vec3d(0, 0, 0);
    m_nu = vec3d(0, 0, 0);
    m_rs = vec2d(0, 0);
    m_Lmd = vec3d(0, 0, 0);
    m_tr = vec3d(0, 0, 0);
    m_epsn = 1.0;
}

//-----------------------------------------------------------------------------
// Helpers for the frame-invariant formulation (see Theory Manual, Section 7.6)
//-----------------------------------------------------------------------------
namespace {

// nodal positions of a surface element at current (prev=false) or previous (prev=true) time
void TiedNodalPositions(FESurface& s, FESurfaceElement& el, bool prev, vec3d* y)
{
    FEMesh& m = *s.GetMesh();
    int ne = el.Nodes();
    bool bshb = s.IsShellBottom();
    for (int i = 0; i < ne; ++i)
    {
        FENode& nd = m.Node(el.m_node[i]);
        if (bshb) y[i] = (prev ? nd.sp() : nd.st());
        else      y[i] = (prev ? nd.m_rp : nd.m_rt);
    }
}

// Kinematics of a tied integration point.
// The contact forces derive from the potential (per unit reference area)
//    Pi = lam.dg + eps/2 dg.dg,   dg = g - Rn.Gap,   lam = Rn.Lmd
// where Rn is the corotational frame (FEContactFrame) of the primary surface at
// the start of the time step. The initial gap and the multiplier are stored by
// their components in this frame, so they rotate with the interface from one
// time step to the next (frame invariance), while acting as dead vectors within
// a time step (which keeps the tangent simple and the formulation stable).
struct TiedKinematics
{
    FEContactFrame Rn;  // corotational frame at start of time step
    vec3d   gap;    // full gap x2 - x1
    vec3d   dg;     // gap - Rn.Gap
    vec3d   t;      // traction (per unit reference area) = Rn.Lmd + eps*dg
    double  J;      // current area Jacobian
    double  J0;     // reference area Jacobian (integration weight)
};

void TiedEvalKinematics(FESurface& ss, FESurfaceElement& se, int j, FESurface& ms,
                        FETiedElasticSurface::Data& pt, double eps, bool prev, TiedKinematics& k)
{
    vec3d ys[FEElement::MAX_NODES], ym[FEElement::MAX_NODES];
    TiedNodalPositions(ss, se, prev, ys);

    // primary point and base vectors
    double* H = se.H(j);
    double* Gr = se.Gr(j);
    double* Gs = se.Gs(j);
    vec3d r(0, 0, 0), g1(0, 0, 0), g2(0, 0, 0);
    for (int i = 0; i < se.Nodes(); ++i) { r += ys[i]*H[i]; g1 += ys[i]*Gr[i]; g2 += ys[i]*Gs[i]; }
    k.J = (g1 ^ g2).norm();

    // tied point on secondary surface
    FESurfaceElement& me = *pt.m_pme;
    TiedNodalPositions(ms, me, prev, ym);
    double Hm[FEElement::MAX_NODES];
    me.shape_fnc(Hm, pt.m_rs[0], pt.m_rs[1]);
    vec3d q(0, 0, 0);
    for (int i = 0; i < me.Nodes(); ++i) q += ym[i]*Hm[i];

    // frame at start of time step
    k.Rn.Evaluate(ss, se, j, true);

    k.gap = q - r;
    k.dg  = k.gap - k.Rn.R(pt.m_Gap);
    k.t   = k.Rn.R(pt.m_Lmd) + k.dg*eps;

    // reference area Jacobian
    vec3d G[2];
    ss.CoBaseVectors0(se, j, G);
    k.J0 = (G[0] ^ G[1]).norm();
}

} // namespace

//-----------------------------------------------------------------------------
// FETiedElasticSurface
//-----------------------------------------------------------------------------

FETiedElasticSurface::FETiedElasticSurface(FEModel* pfem) : FEContactSurface(pfem)
{

}

//-----------------------------------------------------------------------------
bool FETiedElasticSurface::Init()
{
    // initialize surface data first
    if (FEContactSurface::Init() == false) return false;
    
    // allocate node normals
    m_nn.assign(Nodes(), vec3d(0,0,0));
    
    // initialize nodal force vector
    m_Fn.assign(Nodes(), vec3d(0, 0, 0));

    return true;
}

//-----------------------------------------------------------------------------
//! This function calculates the node normal. Due to the piecewise continuity
//! of the surface elements this normal is not uniquely defined so in order to
//! obtain a unique normal the normal is averaged for each node over all the
//! element normals at the node

void FETiedElasticSurface::UpdateNodeNormals()
{
    int N = Nodes(), i, j, ne, jp1, jm1;
    vec3d y[FEElement::MAX_NODES], n;
    
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
//! create material point data
FEMaterialPoint* FETiedElasticSurface::CreateMaterialPoint()
{
	return new FETiedElasticSurface::Data;
}

//-----------------------------------------------------------------------------
void FETiedElasticSurface::Serialize(DumpStream& ar)
{
	FEContactSurface::Serialize(ar);
	ar & m_nn & m_Fn;
}

//-----------------------------------------------------------------------------
void FETiedElasticSurface::GetVectorGap(int nface, vec3d& pg)
{
    FESurfaceElement& el = Element(nface);
    int ni = el.GaussPoints();
    pg = vec3d(0,0,0);
	for (int k = 0; k < ni; ++k)
	{
		Data& data = static_cast<Data&>(*el.GetMaterialPoint(k));
		pg += data.m_dg;
	}
    pg /= ni;
}

//-----------------------------------------------------------------------------
void FETiedElasticSurface::GetContactTraction(int nface, vec3d& pt)
{
    FESurfaceElement& el = Element(nface);
    int ni = el.GaussPoints();
    pt = vec3d(0,0,0);
	for (int k = 0; k < ni; ++k)
	{
		Data& data = static_cast<Data&>(*el.GetMaterialPoint(k));
		pt += data.m_tr;
	}
    pt /= ni;
}

//-----------------------------------------------------------------------------
vec3d FETiedElasticSurface::GetContactForce()
{
    // initialize contact force
    vec3d f(0, 0, 0);

    // loop over all elements of the primary surface
    for (int i = 0; i < m_Fn.size(); ++i) f += m_Fn[i];

    return f;
}

//-----------------------------------------------------------------------------
//! evaluate net contact area
double FETiedElasticSurface::GetContactArea()
{
    // initialize contact area
    double a = 0;

    // loop over all elements of the primary surface
    for (int n = 0; n < Elements(); ++n)
    {
        FESurfaceElement& el = Element(n);
        int nint = el.GaussPoints();

        // evaluate the contact force for that element
        for (int i = 0; i < nint; ++i)
        {
            // get data for this integration point
            Data& data = static_cast<Data&>(*el.GetMaterialPoint(i));
            double T = data.m_tr.norm2();
            if (data.m_pme && (T != 0.0))
            {
                // get the base vectors
                vec3d g[2];
                CoBaseVectors(el, i, g);

                // normal (magnitude = area)
                vec3d n = g[0] ^ g[1];

                // gauss weight
                double w = el.GaussWeights()[i];

                // contact force
                a += n.norm() * w;
            }
        }
    }

    return a;
}

//-----------------------------------------------------------------------------
// FETiedElasticInterface
//-----------------------------------------------------------------------------

FETiedElasticInterface::FETiedElasticInterface(FEModel* pfem) : FEContactInterface(pfem), m_ss(pfem), m_ms(pfem)
{
    static int count = 1;
    SetID(count++);
    
    // initial values
    m_knmult = 1;
    m_atol = 0.1;
    m_epsn = 1;
    m_btwo_pass = false;
    m_stol = 0.01;
    m_bsymm = true;
    m_srad = 1.0;
    m_gtol = -1;    // we use augmentation tolerance by default
    m_bautopen = false;
    m_bupdtpen = false;
    
    m_naugmin = 0;
    m_naugmax = 10;

	m_bflips = false;
	m_bflipm = false;

    // set parents
    m_ss.SetContactInterface(this);
    m_ms.SetContactInterface(this);

    m_ss.SetSibling(&m_ms);
    m_ms.SetSibling(&m_ss);
}

//-----------------------------------------------------------------------------

FETiedElasticInterface::~FETiedElasticInterface()
{
}

//-----------------------------------------------------------------------------
bool FETiedElasticInterface::Init()
{
    // initialize surface data
    if (m_ss.Init() == false) return false;
    if (m_ms.Init() == false) return false;

	// Flip secondary and primary surfaces, if requested.
	// Note that we turn off those flags because otherwise we keep flipping, each time we get here (e.g. in optimization)
	// TODO: Of course, we shouldn't get here more than once. I think we also get through the FEModel::Reset, so I'll have
	//       look into that. 
	if (m_bflips) { m_ss.Invert(); m_bflips = false; }
	if (m_bflipm) { m_ms.Invert(); m_bflipm = false; }
    
    return true;
}

//-----------------------------------------------------------------------------
//! build the matrix profile for use in the stiffness matrix
void FETiedElasticInterface::BuildMatrixProfile(FEGlobalMatrix& K)
{
    FEMesh& mesh = GetMesh();
    
    // get the DOFS
    const int dof_X = GetDOFIndex("x");
    const int dof_Y = GetDOFIndex("y");
    const int dof_Z = GetDOFIndex("z");
    const int dof_RU = GetDOFIndex("Ru");
    const int dof_RV = GetDOFIndex("Rv");
    const int dof_RW = GetDOFIndex("Rw");
    
    const int ndpn = 6;
    vector<int> lm(ndpn*FEElement::MAX_NODES*2);
    
    int npass = (m_btwo_pass?2:1);
    for (int np=0; np<npass; ++np)
    {
        FETiedElasticSurface& ss = (np == 0? m_ss : m_ms);
        
        int ni = 0, k, l;
        for (int j=0; j<ss.Elements(); ++j)
        {
            FESurfaceElement& se = ss.Element(j);
            int nint = se.GaussPoints();
            int* sn = &se.m_node[0];
            for (k=0; k<nint; ++k, ++ni)
            {
				FETiedElasticSurface::Data& pt = static_cast<FETiedElasticSurface::Data&>(*se.GetMaterialPoint(k));
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
                        lm[ndpn*l  ] = id[dof_X];
                        lm[ndpn*l+1] = id[dof_Y];
                        lm[ndpn*l+2] = id[dof_Z];
                        lm[ndpn*l+3] = id[dof_RU];
                        lm[ndpn*l+4] = id[dof_RV];
                        lm[ndpn*l+5] = id[dof_RW];
                    }
                    
                    for (l=0; l<nmeln; ++l)
                    {
                        vector<int>& id = mesh.Node(mn[l]).m_ID;
                        lm[ndpn*(l+nseln)  ] = id[dof_X];
                        lm[ndpn*(l+nseln)+1] = id[dof_Y];
                        lm[ndpn*(l+nseln)+2] = id[dof_Z];
                        lm[ndpn*(l+nseln)+3] = id[dof_RU];
                        lm[ndpn*(l+nseln)+4] = id[dof_RV];
                        lm[ndpn*(l+nseln)+5] = id[dof_RW];
                    }
                    
                    K.build_add(lm);
                }
            }
        }
    }
}

//-----------------------------------------------------------------------------
void FETiedElasticInterface::UpdateAutoPenalty()
{
    // calculate the penalty
    if (m_bautopen)
    {
        CalcAutoPenalty(m_ss);
        if (m_btwo_pass) CalcAutoPenalty(m_ms);
    }
}

//-----------------------------------------------------------------------------
void FETiedElasticInterface::Activate()
{
    // don't forget to call the base class
    FEContactInterface::Activate();
    
    UpdateAutoPenalty();
    
    // project the surfaces onto each other
    // this will evaluate the gap functions in the reference configuration
    InitialProjection(m_ss, m_ms);
    if (m_btwo_pass) InitialProjection(m_ms, m_ss);
}

//-----------------------------------------------------------------------------
void FETiedElasticInterface::CalcAutoPenalty(FETiedElasticSurface& s)
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
			FETiedElasticSurface::Data& pt = static_cast<FETiedElasticSurface::Data&>(*el.GetMaterialPoint(j));
			pt.m_epsn = eps;
        }
    }
}

//-----------------------------------------------------------------------------
// Perform initial projection between tied surfaces in reference configuration
void FETiedElasticInterface::InitialProjection(FETiedElasticSurface& ss, FETiedElasticSurface& ms)
{
    FESurfaceElement* pme;
    vec3d r, nu;
    double rs[2];
    
    // initialize projection data
    FENormalProjection np(ms);
    np.SetTolerance(m_stol);
    np.SetSearchRadius(m_srad);
    np.Init();

	// let's count the number of contact pairs we find.
	int contacts = 0;
    
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
            
            // find the intersection point with the secondary surface
            pme = np.Project2(r, nu, rs);
            
			FETiedElasticSurface::Data& pt = static_cast<FETiedElasticSurface::Data&>(*el.GetMaterialPoint(j));
			pt.m_pme = pme;
            pt.m_rs[0] = rs[0];
            pt.m_rs[1] = rs[1];
            if (pme)
            {
                // the node could potentially be in contact
                // find the global location of the intersection point
                vec3d q = ms.Local2Global(*pme, rs[0], rs[1]);
                
                // calculate the gap function and store its components
                // in the corotational frame of the primary surface
                FEContactFrame f;
                f.Evaluate(ss, el, j, true);
                pt.m_Gap = f.Rt(q - r);

				contacts++;
            }
            else
            {
                // the node is not in contact
                pt.m_Gap = vec3d(0,0,0);
            }
        }
    }

	// if we found no contact pairs, let's report this since this is probably not the user's intention
	if (contacts == 0)
	{
		std::string name = GetName();
		feLogWarning("No contact pairs found for tied interface \"%s\".\nThis contact interface may not have any effect.", name.c_str());
	}
}

//-----------------------------------------------------------------------------
// Evaluate gap functions for position and fluid pressure
void FETiedElasticInterface::ProjectSurface(FETiedElasticSurface& ss, FETiedElasticSurface& ms)
{
    FESurfaceElement* pme;
    vec3d r;
    
    // loop over all integration points
    for (int i=0; i<ss.Elements(); ++i)
    {
        FESurfaceElement& el = ss.Element(i);
        
        int nint = el.GaussPoints();
        
        for (int j=0; j<nint; ++j)
        {
			FETiedElasticSurface::Data& pt = static_cast<FETiedElasticSurface::Data&>(*el.GetMaterialPoint(j));

            // calculate the global position of the integration point
            r = ss.Local2Global(el, j);
            
            // calculate the normal at this integration point
            pt.m_nu = ss.SurfaceNormal(el, j);
            
            // if this node is tied, evaluate gap functions
            pme = pt.m_pme;
            if (pme)
            {
                // find the global location of the intersection point
                vec3d q = ms.Local2Global(*pme, pt.m_rs[0], pt.m_rs[1]);
                
                // calculate the gap function, relative to the initial gap
                // rotated with the corotational frame (frame invariant)
                vec3d g = q - r;
                FEContactFrame f;
                f.Evaluate(ss, el, j, true);
                pt.m_dg = g - f.R(pt.m_Gap);
				pt.m_gap = pt.m_dg.norm();
            }
            else
            {
                // the node is not tied
                pt.m_dg = vec3d(0,0,0);
				pt.m_gap = 0.0;
            }
        }
    }
}

//-----------------------------------------------------------------------------

void FETiedElasticInterface::Update()
{
    // project the surfaces onto each other
    // this will update the gap functions as well
    ProjectSurface(m_ss, m_ms);
    if (m_btwo_pass) ProjectSurface(m_ms, m_ss);

}

//-----------------------------------------------------------------------------
void FETiedElasticInterface::LoadVector(FEGlobalVector& R, const FETimeInfo& tp)
{
    vector<int> sLM, mLM, LM, en;
    vector<double> fe;
    double Hm[FEElement::MAX_NODES];
    

    // zero nodal forces
    m_ss.m_Fn.assign(m_ss.Nodes(), vec3d(0, 0, 0));
    m_ms.m_Fn.assign(m_ms.Nodes(), vec3d(0, 0, 0));

    // Update auto-penalty if requested
    if (m_bupdtpen && (GetFEModel()->GetCurrentStep()->GetFESolver()->m_niter == 0))
        UpdateAutoPenalty();
    
    // loop over the nr of passes
    int npass = (m_btwo_pass?2:1);
    for (int np=0; np<npass; ++np)
    {
        // get primary and secondary surface
        FETiedElasticSurface& ss = (np == 0? m_ss : m_ms);
        FETiedElasticSurface& ms = (np == 0? m_ms : m_ss);
        
        // loop over all primary elements
        for (int i=0; i<ss.Elements(); ++i)
        {
            FESurfaceElement& se = ss.Element(i);
            int nseln = se.Nodes();
            int nint = se.GaussPoints();
            ss.UnpackLM(se, sLM);
            
            for (int j=0; j<nint; ++j)
            {
				FETiedElasticSurface::Data& pt = static_cast<FETiedElasticSurface::Data&>(*se.GetMaterialPoint(j));
                FESurfaceElement* pme = pt.m_pme;
                if (pme == nullptr) continue;
                
                FESurfaceElement& me = *pme;
                int nmeln = me.Nodes();
                ms.UnpackLM(me, mLM);
                
                int ndof = 3*(nseln + nmeln);
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
                en.resize(nseln+nmeln);
                for (int k=0; k<nseln; ++k) en[k      ] = se.m_node[k];
                for (int k=0; k<nmeln; ++k) en[k+nseln] = me.m_node[k];
                
                // shape functions and derivatives
                double* Hs = se.H(j);
                me.shape_fnc(Hm, pt.m_rs[0], pt.m_rs[1]);
                double w = se.GaussWeights()[j];
                
                // penalty
                double eps = m_epsn*pt.m_epsn;
                
                // kinematics, traction and angular momentum correction at t_{n+1}
                TiedKinematics kn;
                TiedEvalKinematics(ss, se, j, ms, pt, eps, false, kn);
                // report traction per current area
                pt.m_tr = kn.t*(kn.J0/kn.J);
                
                // NOTE: as in all other FEBio contact interfaces, the tied constraint
                // forces are evaluated on the end-of-step configuration t_{n+1}
                // (constraint forces at t_{n+1} in the generalized-alpha method, cf.
                // Arnold & Bruls, 2007). This is objective, and was found to be more
                // robust with augmented Lagrangian than alpha-weighted force interpolation.
                vec3d Jt = kn.t*kn.J0;
                
                // nodal forces
                //   f1_a =  J0*N_a t
                //   f2_b = -J0*N_b t
                // These are minus the gradient of the potential
                //   Pi = J0*(lam.dg + eps/2 dg.dg)
                fe.assign(ndof, 0.0);
                for (int k=0; k<nseln; ++k)
                {
                    vec3d f = Jt*(Hs[k]*w);
                    fe[3*k  ] = f.x;
                    fe[3*k+1] = f.y;
                    fe[3*k+2] = f.z;
                }
                for (int k=0; k<nmeln; ++k)
                {
                    vec3d f = -Jt*(Hm[k]*w);
                    fe[3*(k+nseln)  ] = f.x;
                    fe[3*(k+nseln)+1] = f.y;
                    fe[3*(k+nseln)+2] = f.z;
                }
                
                for (int k = 0; k < nseln; ++k) ss.m_Fn[se.m_lnode[k]] += vec3d(fe[3 * k], fe[3 * k + 1], fe[3 * k + 2]);
                for (int k = 0; k < nmeln; ++k) ms.m_Fn[me.m_lnode[k]] += vec3d(fe[3 * nseln + 3 * k], fe[3 * nseln + 3 * k + 1], fe[3 * nseln + 3 * k + 2]);
                
                // assemble the global residual
                R.Assemble(en, LM, fe);
            }
        }
    }
}

//-----------------------------------------------------------------------------
void FETiedElasticInterface::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp)
{
    vector<int> sLM, mLM, LM, en;
    double Hm[FEElement::MAX_NODES];
    FEElementMatrix ke;
    
    
    // do single- or two-pass
    int npass = (m_btwo_pass?2:1);
    for (int np=0; np < npass; ++np)
    {
        FETiedElasticSurface& ss = (np == 0? m_ss : m_ms);
        FETiedElasticSurface& ms = (np == 0? m_ms : m_ss);
        
        for (int i=0; i<ss.Elements(); ++i)
        {
            FESurfaceElement& se = ss.Element(i);
            int nseln = se.Nodes();
            int nint = se.GaussPoints();
            ss.UnpackLM(se, sLM);
            
            for (int j=0; j<nint; ++j)
            {
				FETiedElasticSurface::Data& pt = static_cast<FETiedElasticSurface::Data&>(*se.GetMaterialPoint(j));
                FESurfaceElement* pme = pt.m_pme;
                if (pme == nullptr) continue;
                
                FESurfaceElement& me = *pme;
                int nmeln = me.Nodes();
                ms.UnpackLM(me, mLM);
                
                const int ndpn = 3;
                int ndof = ndpn*(nseln + nmeln);
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
                en.resize(nseln+nmeln);
                for (int k=0; k<nseln; ++k) en[k      ] = se.m_node[k];
                for (int k=0; k<nmeln; ++k) en[k+nseln] = me.m_node[k];
                
                double* Hs = se.H(j);
                me.shape_fnc(Hm, pt.m_rs[0], pt.m_rs[1]);
                double w = se.GaussWeights()[j];
                double eps = m_epsn*pt.m_epsn;
                
                TiedKinematics k;
                TiedEvalKinematics(ss, se, j, ms, pt, eps, false, k);
                const double W = k.J0;      // integration weight (reference area)
                
                // Tangent T = d(f)/d(u) = -Hessian of Pi (symmetric), with
                // f1_a = W*N_a*t, f2_b = -W*N_b*t and d(t)/d(u1_c) = -eps*N_c*I,
                // d(t)/d(u2_d) = eps*N_d*I. FEBio stiffness ke = -w*T.
                ke.resize(ndof, ndof); ke.zero();
                double c = w*W*eps;
                for (int a=0; a<nseln+nmeln; ++a)
                {
                    double Na = (a < nseln ? Hs[a] : -Hm[a-nseln]);
                    for (int b=0; b<nseln+nmeln; ++b)
                    {
                        double Nb = (b < nseln ? Hs[b] : -Hm[b-nseln]);
                        for (int l=0; l<3; ++l) ke[ndpn*a + l][ndpn*b + l] += c*Na*Nb;
                    }
                }
                
                // assemble the global stiffness
				ke.SetNodes(en);
				ke.SetIndices(LM);
				LS.Assemble(ke);
            }
        }
    }
}

//-----------------------------------------------------------------------------
bool FETiedElasticInterface::Augment(int naug, const FETimeInfo& tp)
{
    // make sure we need to augment
    if (m_laugon != FECore::AUGLAG_METHOD) return true;
    
    int i;
    vec3d Ln;
    bool bconv = true;
    
    int NS = m_ss.Elements();
	int NM = m_ms.Elements();
    
    // --- c a l c u l a t e   i n i t i a l   n o r m s ---
    // a. normal component
    double normL0 = 0, normP = 0, normDP = 0;
    for (int i=0; i<NS; ++i)
    {
		FESurfaceElement& el = m_ss.Element(i);
        for (int j=0; j<el.GaussPoints(); ++j)
        {
			FETiedElasticSurface::Data& ds = static_cast<FETiedElasticSurface::Data&>(*el.GetMaterialPoint(j));
			normL0 += ds.m_Lmd*ds.m_Lmd;
        }
    }
    for (int i=0; i<NM; ++i)
    {
		FESurfaceElement& el = m_ms.Element(i);
		for (int j=0; j<el.GaussPoints(); ++j)
        {
			FETiedElasticSurface::Data& dm = static_cast<FETiedElasticSurface::Data&>(*el.GetMaterialPoint(j));
			normL0 += dm.m_Lmd*dm.m_Lmd;
        }
    }
    
    // b. gap component
    // (is calculated during update)
    double maxgap = 0;
    
    // update Lagrange multipliers
    double normL1 = 0, eps;
    for (i=0; i<NS; ++i)
    {
		FESurfaceElement& el = m_ss.Element(i);
		for (int j = 0; j<el.GaussPoints(); ++j)
		{
			FETiedElasticSurface::Data& ds = static_cast<FETiedElasticSurface::Data&>(*el.GetMaterialPoint(j));

            // update Lagrange multipliers on primary surface
            // (components in corotational frame at start of step: Lmd += eps*Rn^T.dg)
            if (ds.m_pme)
            {
                eps = m_epsn*ds.m_epsn;
                FEContactFrame f;
                f.Evaluate(m_ss, el, j, true);
                ds.m_Lmd = ds.m_Lmd + f.Rt(ds.m_dg)*eps;
            }
            
            normL1 += ds.m_Lmd*ds.m_Lmd;
            
            maxgap = max(maxgap,sqrt(ds.m_dg*ds.m_dg));
        }
    }
    
	for (int i = 0; i<NM; ++i)
	{
		FESurfaceElement& el = m_ms.Element(i);
		for (int j = 0; j<el.GaussPoints(); ++j)
		{
			FETiedElasticSurface::Data& dm = static_cast<FETiedElasticSurface::Data&>(*el.GetMaterialPoint(j));

            // update Lagrange multipliers on secondary surface
            // (components in corotational frame at start of step: Lmd += eps*Rn^T.dg)
            if (dm.m_pme)
            {
                eps = m_epsn*dm.m_epsn;
                FEContactFrame f;
                f.Evaluate(m_ms, el, j, true);
                dm.m_Lmd = dm.m_Lmd + f.Rt(dm.m_dg)*eps;
            }
            
            normL1 += dm.m_Lmd*dm.m_Lmd;
            
            maxgap = max(maxgap,sqrt(dm.m_dg*dm.m_dg));
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
    
    if ((m_atol > 0) && (lnorm > m_atol)) bconv = false;
    if ((m_atol > 0) && (pnorm > m_atol)) bconv = false;
    
    if (naug < m_naugmin ) bconv = false;
    if (naug >= m_naugmax) bconv = true;
    
    feLog(" sliding interface # %d\n", GetID());
    feLog("                        CURRENT        REQUIRED\n");
    feLog("    D multiplier : %15le", lnorm); if (m_atol > 0) feLog("%15le\n", m_atol); else feLog("       ***\n");
    
    feLog("    maximum gap  : %15le", maxgap);
    if (m_gtol > 0) feLog("%15le\n", m_gtol); else feLog("       ***\n");
    
    return bconv;
}

//-----------------------------------------------------------------------------
void FETiedElasticInterface::Serialize(DumpStream &ar)
{
    // store contact data
    FEContactInterface::Serialize(ar);
    
    // store contact surface data
    m_ms.Serialize(ar);
    m_ss.Serialize(ar);
    
    // serialize pointers
    if (ar.IsShallow() == false)
    {
		SerializeElementPointers(m_ss, m_ms, ar);
    }
}
