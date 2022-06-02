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
                        lm[ndpn*l+4] = id[dof_RU];
                        lm[ndpn*l+5] = id[dof_RV];
                        lm[ndpn*l+6] = id[dof_RW];
                    }
                    
                    for (l=0; l<nmeln; ++l)
                    {
                        vector<int>& id = mesh.Node(mn[l]).m_ID;
                        lm[ndpn*(l+nseln)  ] = id[dof_X];
                        lm[ndpn*(l+nseln)+1] = id[dof_Y];
                        lm[ndpn*(l+nseln)+2] = id[dof_Z];
                        lm[ndpn*(l+nseln)+4] = id[dof_RU];
                        lm[ndpn*(l+nseln)+5] = id[dof_RV];
                        lm[ndpn*(l+nseln)+6] = id[dof_RW];
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
                
                // calculate the gap function
                pt.m_Gap = q - r;
            }
            else
            {
                // the node is not in contact
                pt.m_Gap = vec3d(0,0,0);
            }
        }
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
                
                // calculate the gap function
                vec3d g = q - r;
                pt.m_dg = g - pt.m_Gap;
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
    int i, j, k;
    vector<int> sLM, mLM, LM, en;
    vector<double> fe;
    const int MN = FEElement::MAX_NODES;
    double detJ[MN], w[MN], *Hs, Hm[MN];
    double N[8*MN];

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
				FETiedElasticSurface::Data& pt = static_cast<FETiedElasticSurface::Data&>(*se.GetMaterialPoint(j));

                // get the secondary element
                FESurfaceElement* pme = pt.m_pme;
                if (pme)
                {
                    // get the secondary element
                    FESurfaceElement& me = *pme;
                    
                    // get the nr of secondary element nodes
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
                    
                    // get primary element shape functions
                    Hs = se.H(j);
                    
                    // get secondary element shape functions
                    double r = pt.m_rs[0];
                    double s = pt.m_rs[1];
                    me.shape_fnc(Hm, r, s);
                    
                    // gap function
                    vec3d dg = pt.m_dg;
                    
                    // lagrange multiplier
                    vec3d Lm = pt.m_Lmd;
                    
                    // penalty
                    double eps = m_epsn*pt.m_epsn;
                    
                    // contact traction
                    vec3d t = Lm + dg*eps;
                    pt.m_tr = t;
                    
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

                    for (int k = 0; k < nseln; ++k) ss.m_Fn[se.m_lnode[k]] += vec3d(fe[3 * k], fe[3 * k + 1], fe[3 * k + 2]);
                    for (int k = 0; k < nmeln; ++k) ms.m_Fn[me.m_lnode[k]] += vec3d(fe[3 * nseln + 3 * k], fe[3 * nseln + 3 * k + 1], fe[3 * nseln + 3 * k + 2]);

                    // assemble the global residual
                    R.Assemble(en, LM, fe);
                }
            }
        }
    }
}

//-----------------------------------------------------------------------------
void FETiedElasticInterface::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp)
{
    int i, j, k, l;
    vector<int> sLM, mLM, LM, en;
    const int MN = FEElement::MAX_NODES;
    double detJ[MN], w[MN], *Hs, Hm[MN];
    FEElementMatrix ke;
    
    // see how many reformations we've had to do so far
    int nref = LS.GetSolver()->m_nref;
    
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
            feLog("Higher order stiffness terms included.\n");
        }
        else knmult = 0;
    }
    
    // do single- or two-pass
    int npass = (m_btwo_pass?2:1);
    for (int np=0; np < npass; ++np)
    {
        // get the primary and secondary surface
        FETiedElasticSurface& ss = (np == 0? m_ss : m_ms);
        FETiedElasticSurface& ms = (np == 0? m_ms : m_ss);
        
        // loop over all primary elements
        for (i=0; i<ss.Elements(); ++i)
        {
            // get ths primary element
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
				FETiedElasticSurface::Data& pt = static_cast<FETiedElasticSurface::Data&>(*se.GetMaterialPoint(j));

                // get the secondary element
                FESurfaceElement* pme = pt.m_pme;
                if (pme)
                {
                    FESurfaceElement& me = *pme;
                    
                    // get the nr of secondary nodes
                    int nmeln = me.Nodes();
                    
                    // copy the LM vector
                    ms.UnpackLM(me, mLM);
                    
                    int ndpn;    // number of dofs per node
                    int ndof;    // number of dofs in stiffness matrix
                    
                    // calculate degrees of freedom for elastic-on-elastic contact
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

                    // build the en vector
                    en.resize(nseln+nmeln);
                    for (k=0; k<nseln; ++k) en[k      ] = se.m_node[k];
                    for (k=0; k<nmeln; ++k) en[k+nseln] = me.m_node[k];
                    
                    // primary shape functions
                    Hs = se.H(j);
                    
                    // secondary shape functions
                    double r = pt.m_rs[0];
                    double s = pt.m_rs[1];
                    me.shape_fnc(Hm, r, s);
                    
                    // get primary normal vector
                    vec3d nu = pt.m_nu;
                    
                    // gap function
                    vec3d dg = pt.m_dg;
                    
                    // lagrange multiplier
                    vec3d Lm = pt.m_Lmd;
                    
                    // penalty
                    double eps = m_epsn*pt.m_epsn;
                    
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
                    
                    vec3d as[FEElement::MAX_NODES];
                    mat3d As[FEElement::MAX_NODES];
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
                    // assemble the global stiffness
					ke.SetNodes(en);
					ke.SetIndices(LM);
					LS.Assemble(ke);
                }
            }
        }
    }
}

//-----------------------------------------------------------------------------
bool FETiedElasticInterface::Augment(int naug, const FETimeInfo& tp)
{
    // make sure we need to augment
    if (m_laugon != 1) return true;
    
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
            eps = m_epsn*ds.m_epsn;
            ds.m_Lmd = ds.m_Lmd + ds.m_dg*eps;
            
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
            eps = m_epsn*dm.m_epsn;
            dm.m_Lmd = dm.m_Lmd + dm.m_dg*eps;
            
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
