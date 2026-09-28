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
#include "FETiedFluidFSISolidTraction.h"
#include "FEFluidMaterial.h"
#include "FEFluidDomain3D.h"
#include <FECore/FENormalProjection.h>
#include <FECore/FEClosestPointProjection.h>
#include <FECore/log.h>
#include <FECore/DumpStream.h>
#include <FECore/FEGlobalMatrix.h>
#include <FECore/FELinearSystem.h>
#include <FECore/FEModel.h>
#include "FEBioFluid.h"

//-----------------------------------------------------------------------------
// Define sliding interface parameters
BEGIN_FECORE_CLASS(FETiedFluidFSISolidTraction, FEContactInterface)
	ADD_PARAMETER(m_laugon   , "laugon")->setLongName("Enforcement method")->setEnums("PENALTY\0AUGLAG\0");
	ADD_PARAMETER(m_atol     , "tolerance"          );
	ADD_PARAMETER(m_gtol     , "gaptol"             );
	ADD_PARAMETER(m_epsn     , "penalty"   );
	ADD_PARAMETER(m_bautopen , "auto_penalty"       );
	ADD_PARAMETER(m_stol     , "search_tol"         );
	ADD_PARAMETER(m_srad     , "search_radius"      );
	ADD_PARAMETER(m_naugmin  , "minaug"             );
	ADD_PARAMETER(m_naugmax  , "maxaug"             );
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
// FETiedFluidFSISolidTraction
//-----------------------------------------------------------------------------

FETiedFluidFSISolidTraction::FETiedFluidFSISolidTraction(FEModel* pfem) : FEContactInterface(pfem), m_s1(pfem), m_s2(pfem), m_dofU(pfem)
{
    static int count = 1;
    SetID(count++);
    
    // initial values
    m_atol = 0.1;
    m_epsn = 1;
    m_stol = 0.01;
    m_srad = 1.0;
    m_gtol = -1;    // we use augmentation tolerance by default
    m_bautopen = false;
    
    m_naugmin = 0;
    m_naugmax = 10;
    
    m_pfluid = nullptr;
    m_psolid = nullptr;
    
    m_solid = 0;

    // set parents
    m_s1.SetContactInterface(this);
    m_s2.SetContactInterface(this);

    m_s1.SetSibling(&m_s2);
    m_s2.SetSibling(&m_s1);
}

//-----------------------------------------------------------------------------

FETiedFluidFSISolidTraction::~FETiedFluidFSISolidTraction()
{
}

//-----------------------------------------------------------------------------
bool FETiedFluidFSISolidTraction::Init()
{
    // initialize surface data
    if (m_s1.Init() == false) return false;
    if (m_s2.Init() == false) return false;
    
    // get the DOFS
    if (m_dofU.AddVariable(FEBioFluid::GetVariableName(FEBioFluid::DISPLACEMENT)) == false) return false;
    
    // Get the fluid-FSI material on either side of the interface, and get the solid material on the other side
    FEFluidFSI* pf1 = GetFluidFSIMaterial(m_s1);
    FEFluidFSI* pf2 = GetFluidFSIMaterial(m_s2);
    FESolidMaterial* ps1 = GetSolidMaterial(m_s1);
    FESolidMaterial* ps2 = GetSolidMaterial(m_s2);
    if ((pf1 == nullptr) && (pf2 == nullptr)) {
        feLogError("Tied fluid-FSI-solid traction interface %d: could not identify a unique fluid-FSI material on either surface.\n", GetID());
        return false;
    }
    else if ((ps1 == nullptr) && (ps2 == nullptr)) {
        feLogError("Tied fluid-FSI-solid traction interface %d: could not identify a unique solid material on either surface.\n", GetID());
        return false;
    }
    m_pfluid = (pf1 == nullptr) ? pf2 : pf1;
    m_psolid = (ps1 == nullptr) ? ps2 : ps1;
    m_solid = (ps1 == nullptr) ? 2 : 1;

    return true;
}

//-----------------------------------------------------------------------------
//! Return the fluid material shared by all elements attached to this surface.
//! Returns nullptr if the elements are not all backed by the same fluid material.
FEFluidFSI* FETiedFluidFSISolidTraction::GetFluidFSIMaterial(FETiedElasticSurface& s)
{
    FEFluidFSI* pfluid = nullptr;
    for (int i=0; i<s.Elements(); ++i)
    {
        FESurfaceElement& el = s.Element(i);
        if (el.m_elem[0].pe == nullptr) return nullptr;
        FEMaterial* pmat = GetFEModel()->GetMaterial(el.m_elem[0].pe->GetMatID());
        FEFluidFSI* pf = dynamic_cast<FEFluidFSI*>(pmat);
        if (pf == nullptr) return nullptr;
        if (pfluid == nullptr) pfluid = pf;
        else if (pfluid != pf) return nullptr;
    }
    return pfluid;
}

//-----------------------------------------------------------------------------
//! Return the fluid material shared by all elements attached to this surface.
//! Returns nullptr if the elements are not all backed by the same fluid material.
FESolidMaterial* FETiedFluidFSISolidTraction::GetSolidMaterial(FETiedElasticSurface& s)
{
    FESolidMaterial* psolid = nullptr;
    for (int i=0; i<s.Elements(); ++i)
    {
        FESurfaceElement& el = s.Element(i);
        if (el.m_elem[0].pe == nullptr) return nullptr;
        FEMaterial* pmat = GetFEModel()->GetMaterial(el.m_elem[0].pe->GetMatID());
        FESolidMaterial* ps = dynamic_cast<FESolidMaterial*>(pmat);
        if (ps == nullptr) return nullptr;
        if (psolid == nullptr) psolid = ps;
        else if (psolid != ps) return nullptr;
    }
    return psolid;
}

//-----------------------------------------------------------------------------
//! build the matrix profile for use in the stiffness matrix
void FETiedFluidFSISolidTraction::BuildMatrixProfile(FEGlobalMatrix& K)
{
    FEMesh& mesh = GetMesh();
    
    const int ndpn = 3;
    
    vector<int> lm(ndpn*FEElement::MAX_NODES*2);
    
    int npass = 2;
    for (int np=0; np<npass; ++np)
    {
        FETiedElasticSurface& s1 = (np == 0? m_s1 : m_s2);
        
        int ni = 0, k, l;
        for (int j=0; j<s1.Elements(); ++j)
        {
            FESurfaceElement& se = s1.Element(j);
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
                    
                    int neln1 = se.Nodes();
                    int neln2 = me.Nodes();
                    
                    for (l=0; l<neln1; ++l)
                    {
                        vector<int>& id = mesh.Node(sn[l]).m_ID;
                        lm[ndpn*l  ] = id[m_dofU[0]];
                        lm[ndpn*l+1] = id[m_dofU[1]];
                        lm[ndpn*l+2] = id[m_dofU[2]];
                    }
                    
                    for (l=0; l<neln2; ++l)
                    {
                        vector<int>& id = mesh.Node(mn[l]).m_ID;
                        lm[ndpn*(l+neln1)  ] = id[m_dofU[0]];
                        lm[ndpn*(l+neln1)+1] = id[m_dofU[1]];
                        lm[ndpn*(l+neln1)+2] = id[m_dofU[2]];
                    }
                    
                    K.build_add(lm);
                }
            }
        }
    }
}

//-----------------------------------------------------------------------------
void FETiedFluidFSISolidTraction::Activate()
{
    // don't forget to call the base class
    FEContactInterface::Activate();
    
    // calculate the penalty
    // regardless whether this is a two-pass analysis or not, the penalty is calculated based on the solid material properties.
    if (m_bautopen) {
        if (m_solid == 1) {
            CalcAutoPenalty(m_s1);
        }
        else {
            CalcAutoPenalty(m_s2);
        }
    }

    // project the surfaces onto each other
    // this will evaluate the gap functions in the reference configuration
    // always perform a two-pass analysis
    InitialProjection(m_s1, m_s2);
    InitialProjection(m_s2, m_s1);
    
    if (m_bautopen) {
        // set the penalty parameter on the fluid-FSI surface to match that of the opposing solid surface
        if (m_solid == 1) SetAutoPenalty(m_s1, m_s2);
        else SetAutoPenalty(m_s2, m_s1);
    }
}

//-----------------------------------------------------------------------------
void FETiedFluidFSISolidTraction::CalcAutoPenalty(FETiedElasticSurface& s)
{
    // loop over all surface elements
    for (int i=0; i<s.Elements(); ++i)
    {
        // get the surface element
        FESurfaceElement& el = s.Element(i);
        
        // calculate a penalty
        double eps = AutoPenalty(el, s);
        
        // assign to integration points of surface element
        int nint = el.GaussPoints();
        for (int j=0; j<nint; ++j)
        {
            FETiedElasticSurface::Data& pt = static_cast<FETiedElasticSurface::Data&>(*el.GetMaterialPoint(j));
			pt.m_epsn = eps;
        }
    }
}

//-----------------------------------------------------------------------------
void FETiedFluidFSISolidTraction::SetAutoPenalty(FETiedElasticSurface& s1, FETiedElasticSurface& s2)
{
    // loop over all surface elements
    for (int i=0; i<s1.Elements(); ++i)
    {
        // get the surface element
        FESurfaceElement& el = s1.Element(i);
        
        // assign to integration points of surface element
        int nint = el.GaussPoints();
        for (int j=0; j<nint; ++j)
        {
            FETiedElasticSurface::Data& data = static_cast<FETiedElasticSurface::Data&>(*el.GetMaterialPoint(j));
            if (data.m_pme) {
                FETiedElasticSurface::Data& pt = static_cast<FETiedElasticSurface::Data&>(*data.m_pme->GetMaterialPoint(j));
                pt.m_epsn = data.m_epsn;
            }
        }
    }
}

//-----------------------------------------------------------------------------
// Perform initial projection between tied surfaces in reference configuration
void FETiedFluidFSISolidTraction::InitialProjection(FETiedElasticSurface& s1, FETiedElasticSurface& s2)
{
    FEMesh& mesh = GetMesh();
    FESurfaceElement* pme;
    vec3d r, nu;
    double rs[2];
    
    // initialize projection data
    FENormalProjection np(s2);
    np.SetTolerance(m_stol);
    np.SetSearchRadius(m_srad);
    np.Init();
    
    FEClosestPointProjection cp(s1);
    cp.SetTolerance(m_stol);
    cp.SetSearchRadius(m_srad);
    cp.Init();
    vec3d sq;
    vec2d srs;

    // projection diagnostics
    int nproj = 0, nfail = 0;
    double maxgap = 0;
    
    // loop over all integration points
    int n = 0;
    for (int i=0; i<s1.Elements(); ++i)
    {
        FESurfaceElement& el = s1.Element(i);
        
        int nint = el.GaussPoints();
        
        for (int j=0; j<nint; ++j, ++n)
        {
            // calculate the global position of the integration point
            r = s1.Local2Global(el, j);
            
            // calculate the normal at this integration point
            nu = s1.SurfaceNormal(el, j);
            
            // find the intersection point with the secondary surface
            pme = np.Project2(r, nu, rs);
            
            FETiedElasticSurface::Data& pt = static_cast<FETiedElasticSurface::Data&>(*el.GetMaterialPoint(j));
			pt.m_pme = pme;
            pt.m_nu = nu;
            pt.m_rs[0] = rs[0];
            pt.m_rs[1] = rs[1];
            if (pme)
            {
                // the node could potentially be in contact
                // find the global location of the intersection point
                vec3d q = s2.Local2Global(*pme, rs[0], rs[1]);
                
                // calculate the gap function
                pt.m_Gap = q - r;
                ++nproj;
                maxgap = max(maxgap, pt.m_Gap.norm());
            }
            else
            {
                // the integration point could not be projected onto the opposing surface
                pt.m_Gap = vec3d(0,0,0);
                ++nfail;
            }
        }
    }
    
    // report the outcome of the projection. Integration points that could not be
    // projected are left untied: they revert to the natural boundary conditions
    feLog(" tied elastic interface # %d:\n", GetID());
    feLog("    tied integration points  : %d\n", nproj);
    feLog("    maximum initial gap      : %15le\n", maxgap);
    if (nfail > 0) {
        feLogWarning("Tied elastic interface %d: %d integration point(s) could not be projected onto the\n"
                     "opposing surface. These points remain untied and behave as a frictionless\n"
                     "impermeable wall. Consider increasing search_radius or search_tol.", GetID(), nfail);
    }
}

//-----------------------------------------------------------------------------
// Evaluate gap functions for fluid velocity and fluid pressure
void FETiedFluidFSISolidTraction::ProjectSurface(FETiedElasticSurface& s1, FETiedElasticSurface& s2)
{
    FEMesh& mesh = GetMesh();
    FESurfaceElement* pme;
    vec3d r;
    double alpha = GetFEModel()->GetTime().alphaf;
    
    vec3d  ut[FEElement::MAX_NODES], up[FEElement::MAX_NODES], u1;
    
    // loop over all integration points
    for (int i=0; i<s1.Elements(); ++i)
    {
        FESurfaceElement& el = s1.Element(i);
        
        int ne = el.Nodes();
        int nint = el.GaussPoints();
        
        // get the nodal velocities and dilatations
        for (int j=0; j<ne; ++j) {
            FENode& node = mesh.Node(el.m_node[j]);
            ut[j] = node.get_vec3d(m_dofU[0], m_dofU[1], m_dofU[2]);
            up[j] = node.get_vec3d_prev(m_dofU[0], m_dofU[1], m_dofU[2]);
        }

        for (int j=0; j<nint; ++j)
        {
            FETiedElasticSurface::Data& pt = static_cast<FETiedElasticSurface::Data&>(*el.GetMaterialPoint(j));

            // calculate the global position of the integration point
            r = s1.Local2Global(el, j);
            
            // get the velocity and dilatation at the integration point
            u1 = el.eval(ut, j)*alpha + el.eval(up, j)*(1-alpha);
            
            // if this node is tied, evaluate gap functions
            pme = pt.m_pme;
            if (pme)
            {
                // calculate the vectorial gap function
                vec3d umt[FEElement::MAX_NODES], ump[FEElement::MAX_NODES];
                for (int k=0; k<pme->Nodes(); ++k) {
                    FENode& node = mesh.Node(pme->m_node[k]);
                    umt[k] = node.get_vec3d(m_dofU[0], m_dofU[1], m_dofU[2]);
                    ump[k] = node.get_vec3d_prev(m_dofU[0], m_dofU[1], m_dofU[2]);
                }
                vec3d u2 = pme->eval(umt, pt.m_rs[0], pt.m_rs[1])*alpha + pme->eval(ump, pt.m_rs[0], pt.m_rs[1])*(1-alpha);
                pt.m_Gap = u2 - u1;

                // penalty factors
                double epsn = m_epsn*pt.m_epsn;
                
                // viscous traction and normal velocity (augmented Lagrangian form)
                pt.m_tr = pt.m_Lmd + pt.m_Gap*epsn;
            }
            else
            {
                // the node is not tied
                pt.m_Gap = vec3d(0,0,0);
                pt.m_tr = vec3d(0,0,0);
            }
        }
    }
}

//-----------------------------------------------------------------------------

void FETiedFluidFSISolidTraction::Update()
{
    // project the surfaces onto each other
    // this will update the gap functions as well
    ProjectSurface(m_s1, m_s2);
    ProjectSurface(m_s2, m_s1);
}

//-----------------------------------------------------------------------------
void FETiedFluidFSISolidTraction::LoadVector(FEGlobalVector& R, const FETimeInfo& tp)
{
    vector<int> LM1, LM2, LM, en;
    vector<double> fe;
    const int MI = FEElement::MAX_INTPOINTS;
    const int MN = FEElement::MAX_NODES;
    double detJ[MI], w[MI], *H1, H2[MN];
    vec3d f1[MN], f2[MN];
    double w1[MN], w2[MN];

    // loop over the nr of passes
    int npass = 2;
    for (int np=0; np<npass; ++np)
    {
        // get primary and secondary surfaces
        FETiedElasticSurface& s1 = (np == 0? m_s1 : m_s2);
        FETiedElasticSurface& s2 = (np == 0? m_s2 : m_s1);
        
        // loop over all elements of primary surface
        for (int i=0; i<s1.Elements(); ++i)
        {
            // get the surface element
            FESurfaceElement& se1 = s1.Element(i);
            
            // get the nr of nodes and integration points
            int neln1 = se1.Nodes();
            int nint1 = se1.GaussPoints();
            
            // copy the LM vector; we'll need it later
            s1.UnpackLM(se1, LM1);
            
            // we calculate all the metrics we need before we
            // calculate the nodal forces
            for (int j=0; j<nint1; ++j)
            {
                // get the base vectors
                vec3d g[2];
                s1.CoBaseVectors(se1, j, g);
                
                // jacobians: J = |g0xg1|
                detJ[j] = (g[0] ^ g[1]).norm();
                
                // integration weights
                w[j] = se1.GaussWeights()[j];
            }
            
            // loop over all integration points
            // note that we are integrating over the current surface
            for (int j=0; j<nint1; ++j)
            {
                FETiedElasticSurface::Data& pt = static_cast<FETiedElasticSurface::Data&>(*se1.GetMaterialPoint(j));

                // get the secondary surface element
                FESurfaceElement* pme = pt.m_pme;
                if (pme)
                {
                    // get the secondary surface element
                    FESurfaceElement& se2 = *pme;
                    
                    // get the nr of secondary element nodes
                    int neln2 = se2.Nodes();
                    
                    // copy LM vector
                    s2.UnpackLM(se2, LM2);
                    
                    // calculate degrees of freedom
                    const int ndpn = 3;
                    int ndof = ndpn*(neln1 + neln2);
                    
                    // build the LM vector
                    LM.resize(ndof);
                    for (int a=0; a<neln1; ++a)
                        for (int k=0; k< ndpn; ++k) LM[ndpn*a+k] = LM1[ndpn*a+k];
                    
                    for (int b=0; b<neln2; ++b)
                        for (int k=0; k< ndpn; ++k) LM[ndpn*(b+neln1)+k] = LM2[ndpn*b+k];
                    
                    // build the en vector
                    en.resize(neln1+neln2);
                    for (int a=0; a<neln1; ++a) en[a      ] = se1.m_node[a];
                    for (int b=0; b<neln2; ++b) en[b+neln1] = se2.m_node[b];
                    
                    // get primary element shape functions
                    H1 = se1.H(j);
                    
                    // get secondary element shape functions
                    double r = pt.m_rs[0];
                    double s = pt.m_rs[1];
                    se2.shape_fnc(H2, r, s);
                    
                    // contact traction
                    // (evaluated in ProjectSurface, called from Update)
                    vec3d tr = pt.m_tr;
                    
                    // calculate the force vector
                    fe.resize(ndof);
                    zero(fe);
                    
                    for (int a=0; a<neln1; ++a) f1[a] = tr*H1[a];
                    for (int b=0; b<neln2; ++b) f2[b] = -tr*H2[b];
                    
                    for (int a=0; a<neln1; ++a)
                    {
                        fe[ndpn*a  ] += f1[a].x*detJ[j]*w[j];
                        fe[ndpn*a+1] += f1[a].y*detJ[j]*w[j];
                        fe[ndpn*a+2] += f1[a].z*detJ[j]*w[j];
                    }
                    for (int b = 0; b<neln2; ++b) {
                        fe[ndpn*(b+neln1)  ] += f2[b].x*detJ[j]*w[j];
                        fe[ndpn*(b+neln1)+1] += f2[b].y*detJ[j]*w[j];
                        fe[ndpn*(b+neln1)+2] += f2[b].z*detJ[j]*w[j];
                    }

                    // assemble the global residual
                    R.Assemble(en, LM, fe);
                }
            }
        }
    }
}

//-----------------------------------------------------------------------------
void FETiedFluidFSISolidTraction::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp)
{
    vector<int> LM1, LM2, LM, en;
    const int MI = FEElement::MAX_INTPOINTS;
    const int MN = FEElement::MAX_NODES;
    double detJ[MI], w[MI], *H1, H2[MN], dH1r, dH1s;
    FEElementMatrix ke;
    mat3da Ac[MI];
    vec3d nu[MI];
    
    double alpha = tp.alphaf;
    
    // do single- or two-pass
    int npass = 2;
    for (int np=0; np < npass; ++np)
    {
		// get primary and secondary surfaces
		FETiedElasticSurface& s1 = (np == 0? m_s1 : m_s2);
        FETiedElasticSurface& s2 = (np == 0? m_s2 : m_s1);
        
        // loop over all elements of primary surface
        for (int i=0; i<s1.Elements(); ++i)
        {
            // get the next element
            FESurfaceElement& se1 = s1.Element(i);
            
            // get nr of nodes and integration points
            int neln1 = se1.Nodes();
            int nint1 = se1.GaussPoints();
            
            // copy the LM vector
            s1.UnpackLM(se1, LM1);
            
            // we calculate all the metrics we need before we
            // calculate the nodal forces
            for (int c=0; c<nint1; ++c)
            {
                // get the base vectors
                vec3d g[2];
                s1.CoBaseVectors(se1, c, g);
                
                nu[c] = g[0] ^ g[1];

                // jacobians: J = |g0xg1|
                detJ[c] = nu[c].unit();
                
                // integration weights
                w[c] = se1.GaussWeights()[c];
                
                // primary shape function derivatives
                dH1r = se1.gr(c);
                dH1s = se1.gs(c);
                vec3d omega = g[1]*dH1r - g[0]*dH1s;
                Ac[c] = mat3da(omega);
            }
            
            // loop over all integration points
            for (int j=0; j<nint1; ++j)
            {
                FETiedElasticSurface::Data& pt = static_cast<FETiedElasticSurface::Data&>(*se1.GetMaterialPoint(j));

                // get the secondary element
                FESurfaceElement* pme = pt.m_pme;
                if (pme)
                {
                    FESurfaceElement& se2 = *pme;
                    
                    // get the nr of secondary nodes
                    int neln2 = se2.Nodes();
                    
                    // copy the LM vector
                    s2.UnpackLM(se2, LM2);
                    
                    int ndpn;    // number of dofs per node
                    int ndof;    // number of dofs in stiffness matrix
                    
                    // calculate degrees of freedom for elastic-on-elastic contact
                    ndpn = 3;
                    ndof = ndpn*(neln1 + neln2);
                    
                    // build the LM vector
                    LM.resize(ndof);
                    
                    for (int a=0; a<neln1; ++a)
                        for (int k=0; k<ndpn; ++k) LM[ndpn*a+k] = LM1[ndpn*a+k];
                    
                    for (int b=0; b<neln2; ++b)
                        for (int k=0; k<ndpn; ++k) LM[ndpn*(b+neln1)+k] = LM2[ndpn*b+k];

                    // build the en vector
                    en.resize(neln1+neln2);
                    for (int a=0; a<neln1; ++a) en[a      ] = se1.m_node[a];
                    for (int b=0; b<neln2; ++b) en[b+neln1] = se2.m_node[b];
                    
                    // primary shape functions
                    H1 = se1.H(j);
                    
                    // secondary shape functions
                    double r = pt.m_rs[0];
                    double s = pt.m_rs[1];
                    se2.shape_fnc(H2, r, s);
                    
                    // penalty
                    double epsn = m_epsn*pt.m_epsn;
                    
                    // create the stiffness matrix
                    ke.resize(ndof, ndof); ke.zero();
                    
                    //------------------------------------
                    
                    // NOTE: The dilatation blocks k11..k22 have the opposite sign of the
                    // velocity blocks K11..K22, because the dilatation gap is now
                    // pi = J(2) - J(1) whereas the velocity gap is g = v(2) - v(1) and
                    // both enter the residual with the same sign. Both the velocity and
                    // the dilatation contributions to ke are then positive semi-definite,
                    // as a penalty contribution must be.
                    
                    for (int a=0; a<neln1; ++a) {
                        for (int c=0; c<neln1; ++c)
                        {
                            mat3d K11(((pt.m_Gap*H1[a] & (Ac[c]*nu[c]))-mat3dd(H1[a]*H1[c]*detJ[j]))*(epsn*alpha*w[j]));
                            ke[ndpn*a  ][ndpn*c  ] -= K11(0,0); ke[ndpn*a  ][ndpn*c+1] -= K11(0,1); ke[ndpn*a  ][ndpn*c+2] -= K11(0,2);
                            ke[ndpn*a+1][ndpn*c  ] -= K11(1,0); ke[ndpn*a+1][ndpn*c+1] -= K11(1,1); ke[ndpn*a+1][ndpn*c+2] -= K11(1,2);
                            ke[ndpn*a+2][ndpn*c  ] -= K11(2,0); ke[ndpn*a+2][ndpn*c+1] -= K11(2,1); ke[ndpn*a+2][ndpn*c+2] -= K11(2,2);
                        }
                        for (int d=0; d<neln2; ++d)
                        {
                            mat3dd K12(epsn*H1[a]*H2[d]*detJ[j]*w[j]*alpha);
                            ke[ndpn*a    ][ndpn*(neln1+d)    ] -= K12.xx();
                            ke[ndpn*a + 1][ndpn*(neln1+d) + 1] -= K12.yy();
                            ke[ndpn*a + 2][ndpn*(neln1+d) + 2] -= K12.zz();
                        }
                    }

                    for (int b=0; b<neln2; ++b) {
                        for (int c=0; c<neln1; ++c)
                        {
                            mat3d K21(mat3dd(H2[b]*H1[c]*detJ[j])-(pt.m_Gap*H2[b] & (Ac[c]*nu[c]))*(epsn*alpha*w[j]));
                            ke[ndpn*(neln1+b)  ][ndpn*c  ] -= K21(0,0); ke[ndpn*(neln1+b)  ][ndpn*c+1] -= K21(0,1); ke[ndpn*(neln1+b)  ][ndpn*c+2] -= K21(0,2);
                            ke[ndpn*(neln1+b)+1][ndpn*c  ] -= K21(1,0); ke[ndpn*(neln1+b)+1][ndpn*c+1] -= K21(1,1); ke[ndpn*(neln1+b)+1][ndpn*c+2] -= K21(1,2);
                            ke[ndpn*(neln1+b)+2][ndpn*c  ] -= K21(2,0); ke[ndpn*(neln1+b)+2][ndpn*c+1] -= K21(2,1); ke[ndpn*(neln1+b)+2][ndpn*c+2] -= K21(2,2);
                        }
                        for (int d=0; d<neln2; ++d)
                        {
                            mat3dd K22(-epsn*H2[b]*H2[d]*detJ[j]*w[j]*alpha);
                            ke[ndpn*(neln1+b)    ][ndpn*(neln1+d)    ] -= K22.xx();
                            ke[ndpn*(neln1+b) + 1][ndpn*(neln1+d) + 1] -= K22.yy();
                            ke[ndpn*(neln1+b) + 2][ndpn*(neln1+d) + 2] -= K22.zz();
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
bool FETiedFluidFSISolidTraction::Augment(int naug, const FETimeInfo& tp)
{
    // make sure we need to augment
	if (m_laugon != FECore::AUGLAG_METHOD) return true;

    bool bconv = true;
    
    int N1 = m_s1.Elements();
    int N2 = m_s2.Elements();
    
    // --- c a l c u l a t e   i n i t i a l   n o r m s ---
    // a. normal component
    double normL0 = 0;
    for (int i=0; i<N1; ++i)
    {
		FESurfaceElement& s1 = m_s1.Element(i);
        for (int j=0; j<s1.GaussPoints(); ++j)
        {
			FETiedElasticSurface::Data& d1 = static_cast<FETiedElasticSurface::Data&>(*s1.GetMaterialPoint(j));
			normL0 += d1.m_Lmd*d1.m_Lmd;
        }
    }
    for (int i=0; i<N2; ++i)
    {
		FESurfaceElement& s2 = m_s2.Element(i);
        for (int j=0; j<s2.GaussPoints(); ++j)
        {
            FETiedElasticSurface::Data& d2 = static_cast<FETiedElasticSurface::Data&>(*s2.GetMaterialPoint(j));
			normL0 += d2.m_Lmd*d2.m_Lmd;
        }
    }
    normL0 = sqrt(normL0);
    
    // b. gap component
    // (is calculated during update)
    double maxgap = 0;
    
    // update Lagrange multipliers
    double normL1 = 0, normJ1 = 0, epsn;
    for (int i=0; i<N1; ++i)
    {
		FESurfaceElement& s1 = m_s1.Element(i);
		for (int j = 0; j<s1.GaussPoints(); ++j)
		{
            FETiedElasticSurface::Data& d1 = static_cast<FETiedElasticSurface::Data&>(*s1.GetMaterialPoint(j));

            if (d1.m_pme) {
                // update Lagrange multipliers on primary surface
                epsn = m_epsn*d1.m_epsn;
                d1.m_Lmd = d1.m_Lmd + d1.m_Gap*epsn;
                maxgap = max(maxgap,sqrt(d1.m_Gap*d1.m_Gap));
                normL1 += d1.m_Lmd*d1.m_Lmd;
                
                // keep the reported traction and normal velocity consistent with the
                // updated multipliers (they are re-evaluated on the next Update)
                d1.m_tr = d1.m_Lmd + d1.m_Gap*epsn;
            }
        }
    }
    
    for (int i=0; i<N2; ++i)
    {
		FESurfaceElement& s2 = m_s2.Element(i);
		for (int j = 0; j<s2.GaussPoints(); ++j)
		{
            FETiedElasticSurface::Data& d2 = static_cast<FETiedElasticSurface::Data&>(*s2.GetMaterialPoint(j));

            if (d2.m_pme) {
                // update Lagrange multipliers on secondary surface
                epsn = m_epsn*d2.m_epsn;
                d2.m_Lmd = d2.m_Lmd + d2.m_Gap*epsn;
                maxgap = max(maxgap,sqrt(d2.m_Gap*d2.m_Gap));
                normL1 += d2.m_Lmd*d2.m_Lmd;
                
                d2.m_tr = d2.m_Lmd + d2.m_Gap*epsn;
            }
        }
    }
    normL1 = sqrt(normL1);
    
    // calculate relative norms
    double lnorm = (normL1 != 0 ? fabs((normL1 - normL0) / normL1) : fabs(normL1 - normL0));
    
    // check convergence
    if ((m_gtol > 0) && (maxgap > m_gtol)) bconv = false;
    
    if ((m_atol > 0) && (lnorm > m_atol)) bconv = false;
    
    if (naug < m_naugmin ) bconv = false;
    if (naug >= m_naugmax) bconv = true;
    
    feLog(" tied fluid-FSI-solid traction interface # %d\n", GetID());
    feLog("                                CURRENT        REQUIRED\n");
    feLog("    solid multiplier : %15le", lnorm); if (m_atol > 0) feLog("%15le\n", m_atol); else feLog("       ***\n");
    feLog("    maximum gap      : %15le", maxgap);
    if (m_gtol > 0) feLog("%15le\n", m_gtol); else feLog("       ***\n");

    return bconv;
}

//-----------------------------------------------------------------------------
void FETiedFluidFSISolidTraction::Serialize(DumpStream &ar)
{
    // store contact data
    FEContactInterface::Serialize(ar);
    
    // store contact surface data
    m_s1.Serialize(ar);
    m_s2.Serialize(ar);

	// serialize element pointers
	SerializeElementPointers(m_s1, m_s2, ar);
	SerializeElementPointers(m_s2, m_s1, ar);
    
    if (ar.IsShallow()) return;
    ar & m_pfluid & m_psolid;
    ar & m_dofU;
}
