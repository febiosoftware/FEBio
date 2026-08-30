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
#include "FETiedFluidInterface.h"
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
BEGIN_FECORE_CLASS(FETiedFluidInterface, FEContactInterface)
	ADD_PARAMETER(m_laugon   , "laugon")->setLongName("Enforcement method")->setEnums("PENALTY\0AUGLAG\0");
	ADD_PARAMETER(m_atol     , "tolerance"          );
	ADD_PARAMETER(m_gtol     , "vtol"               );
	ADD_PARAMETER(m_etol     , "etol"               );
	ADD_PARAMETER(m_epst     , "traction_penalty"   );
	ADD_PARAMETER(m_bautopen , "auto_penalty"       );
	ADD_PARAMETER(m_btwo_pass, "two_pass"           );
	ADD_PARAMETER(m_stol     , "search_tol"         );
	ADD_PARAMETER(m_epsn     , "dilatation_penalty"   );
	ADD_PARAMETER(m_srad     , "search_radius"      );
	ADD_PARAMETER(m_naugmin  , "minaug"             );
	ADD_PARAMETER(m_naugmax  , "maxaug"             );
    ADD_PARAMETER(m_bfreedofs, "free_dofs"          );
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FETiedFluidSurface::Data::Data()
{
    m_Gap = vec3d(0,0,0);
    m_vg  = vec3d(0,0,0);
    m_nu  = vec3d(0,0,0);
    m_rs  = vec2d(0,0);
    m_Lmd = vec3d(0,0,0);
    m_tv  = vec3d(0,0,0);
    m_Lmp = 0.0;
    m_epst= 1.0;
    m_epsn= 1.0;
    m_Jg  = 0.0;
    m_vn = 0;
    m_pme = (FESurfaceElement*)0;
}

//-----------------------------------------------------------------------------
void FETiedFluidSurface::Data::Serialize(DumpStream& ar)
{
	FEContactMaterialPoint::Serialize(ar);
	ar & m_Gap;
	ar & m_vg;
	ar & m_nu;
	ar & m_rs;
	ar & m_Lmd;
    ar & m_tv;
    ar & m_Lmp;
	ar & m_epst;
	ar & m_epsn;
	ar & m_Jg;
    ar & m_vn;
}

//-----------------------------------------------------------------------------
// FETiedFluidSurface
//-----------------------------------------------------------------------------

FETiedFluidSurface::FETiedFluidSurface(FEModel* pfem) : FEContactSurface(pfem), m_dofWE(pfem)
{
}

//-----------------------------------------------------------------------------
bool FETiedFluidSurface::Init()
{
    // initialize surface data first
    if (FEContactSurface::Init() == false) return false;
    
	// set the dof list
	if (m_dofWE.AddVariable(FEBioFluid::GetVariableName(FEBioFluid::RELATIVE_FLUID_VELOCITY)) == false) return false;
    if (m_dofWE.AddVariable(FEBioFluid::GetVariableName(FEBioFluid::FLUID_DILATATION)) == false) return false;

    return true;
}

//-----------------------------------------------------------------------------
void FETiedFluidSurface::UnpackLM(FEElement& el, vector<int>& lm)
{
    int N = el.Nodes();
    lm.resize(N*4);
    for (int i=0; i<N; ++i)
    {
        int n = el.m_node[i];
        FENode& node = m_pMesh->Node(n);
        vector<int>& id = node.m_ID;
        
        lm[4*i  ] = id[m_dofWE[0]];
        lm[4*i+1] = id[m_dofWE[1]];
        lm[4*i+2] = id[m_dofWE[2]];
        lm[4*i+3] = id[m_dofWE[3]];
    }
}

//-----------------------------------------------------------------------------
//! create material point data
FEMaterialPoint* FETiedFluidSurface::CreateMaterialPoint()
{
	return new FETiedFluidSurface::Data;
}

//-----------------------------------------------------------------------------
void FETiedFluidSurface::GetVelocityGap(int nface, vec3d& vg)
{
    FESurfaceElement& el = Element(nface);
    int ni = el.GaussPoints();
    vg = vec3d(0,0,0);
	for (int k = 0; k < ni; ++k)
	{
		Data& d = static_cast<Data&>(*el.GetMaterialPoint(k));
		vg += d.m_vg;
	}
    vg /= ni;
}

//-----------------------------------------------------------------------------
void FETiedFluidSurface::GetDilatationGap(int nface, double& Jg)
{
    FESurfaceElement& el = Element(nface);
    int ni = el.GaussPoints();
    Jg = 0;
	for (int k = 0; k < ni; ++k)
	{
		Data& d = static_cast<Data&>(*el.GetMaterialPoint(k));
		Jg += d.m_Jg;
	}
    Jg /= ni;
}

//-----------------------------------------------------------------------------
void FETiedFluidSurface::GetViscousTraction(int nface, vec3d& tv)
{
    FESurfaceElement& el = Element(nface);
    int ni = el.GaussPoints();
    tv = vec3d(0,0,0);
	for (int k = 0; k < ni; ++k)
	{
		Data& d = static_cast<Data&>(*el.GetMaterialPoint(k));
		tv += d.m_tv;
	}
    tv /= ni;
}

//-----------------------------------------------------------------------------
void FETiedFluidSurface::GetNormalVelocity(int nface, double& vn)
{
    FESurfaceElement& el = Element(nface);
    int ni = el.GaussPoints();
    vn = 0;
	for (int k = 0; k < ni; ++k)
	{
		Data& d = static_cast<Data&>(*el.GetMaterialPoint(k));
		vn += d.m_vn;
	}
    vn /= ni;
}

//-----------------------------------------------------------------------------
double FETiedFluidSurface::GetArea(FESurfaceElement& el)
{
    int ni = el.GaussPoints();
    vec3d rt[FEElement::MAX_NODES];
    GetNodalCoordinates(el, 1.0, rt);
    double* gw = el.GaussWeights();
    double area = 0;
    for (int i=0; i<ni;++i) {
        FEMaterialPoint& mp = *el.GetMaterialPoint(i);
        vec3d dxr = el.eval_deriv1(rt, mp.m_index);
        vec3d dxs = el.eval_deriv2(rt, mp.m_index);
        // normal and area element
        vec3d n = dxr ^ dxs;
        double da = n.unit();
        area += da*gw[i];
    }
    return area;
}

//-----------------------------------------------------------------------------
double FETiedFluidSurface::GetVolume(FESolidElement& el)
{
    int ni = el.GaussPoints();
    double* gw = el.GaussWeights();
    double volume = 0;
    for (int i=0; i<ni; ++i) {
        double detJ = 1./(el.m_J0i[i].det());
        volume += detJ*gw[i];
    }
    return volume;
}

//-----------------------------------------------------------------------------
// FETiedFluidInterface
//-----------------------------------------------------------------------------

FETiedFluidInterface::FETiedFluidInterface(FEModel* pfem) : FEContactInterface(pfem), m_s1(pfem), m_s2(pfem), m_dofWE(pfem)
{
    static int count = 1;
    SetID(count++);
    
    // initial values
    m_atol = 0.1;
    m_epst = 1;
    m_epsn = 1;
    m_btwo_pass = false;
    m_stol = 0.01;
    m_srad = 1.0;
    m_gtol = -1;    // we use augmentation tolerance by default
    m_etol = -1;    // we use augmentation tolerance by default
    m_bautopen = false;
    m_bfreedofs = false;
    
    m_naugmin = 0;
    m_naugmax = 10;
    
    // set parents
    m_s1.SetContactInterface(this);
    m_s2.SetContactInterface(this);

    m_s1.SetSibling(&m_s2);
    m_s2.SetSibling(&m_s1);
}

//-----------------------------------------------------------------------------

FETiedFluidInterface::~FETiedFluidInterface()
{
}

//-----------------------------------------------------------------------------
bool FETiedFluidInterface::Init()
{
    // initialize surface data
    if (m_s1.Init() == false) return false;
    if (m_s2.Init() == false) return false;
    
    int N = m_s1.Nodes() + m_s2.Nodes();
    
    // get the DOFS
    m_dofWE.AddVariable(FEBioFluid::GetVariableName(FEBioFluid::RELATIVE_FLUID_VELOCITY));
    m_dofWE.AddVariable(FEBioFluid::GetVariableName(FEBioFluid::FLUID_DILATATION));
    
    return true;
}

//-----------------------------------------------------------------------------
//! build the matrix profile for use in the stiffness matrix
void FETiedFluidInterface::BuildMatrixProfile(FEGlobalMatrix& K)
{
    FEMesh& mesh = GetMesh();
    
    vector<int> lm(4*FEElement::MAX_NODES*2);
    
    int npass = (m_btwo_pass?2:1);
    for (int np=0; np<npass; ++np)
    {
        FETiedFluidSurface& s1 = (np == 0? m_s1 : m_s2);
        
        int ni = 0, k, l;
        for (int j=0; j<s1.Elements(); ++j)
        {
            FESurfaceElement& se = s1.Element(j);
            int nint = se.GaussPoints();
            int* sn = &se.m_node[0];
            for (k=0; k<nint; ++k, ++ni)
            {
                FETiedFluidSurface::Data& pt = static_cast<FETiedFluidSurface::Data&>(*se.GetMaterialPoint(k));
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
                        lm[4*l  ] = id[m_dofWE[0]];
                        lm[4*l+1] = id[m_dofWE[1]];
                        lm[4*l+2] = id[m_dofWE[2]];
                        lm[4*l+3] = id[m_dofWE[3]];
                    }
                    
                    for (l=0; l<neln2; ++l)
                    {
                        vector<int>& id = mesh.Node(mn[l]).m_ID;
                        lm[4*(l+neln1)  ] = id[m_dofWE[0]];
                        lm[4*(l+neln1)+1] = id[m_dofWE[1]];
                        lm[4*(l+neln1)+2] = id[m_dofWE[2]];
                        lm[4*(l+neln1)+3] = id[m_dofWE[3]];
                    }
                    
                    K.build_add(lm);
                }
            }
        }
    }
}

//-----------------------------------------------------------------------------
void FETiedFluidInterface::Activate()
{
    // don't forget to call the base class
    FEContactInterface::Activate();
    
    // calculate the penalty
    if (m_bautopen)
    {
        CalcAutoViscousTractionPenalty(m_s1);
        CalcAutoNormalVelocityPenalty(m_s1);
    }
    
    // project the surfaces onto each other
    // this will evaluate the gap functions in the reference configuration
    InitialProjection(m_s1, m_s2);
    if (m_btwo_pass) InitialProjection(m_s2, m_s1);
}

//-----------------------------------------------------------------------------
void FETiedFluidInterface::CalcAutoViscousTractionPenalty(FETiedFluidSurface& s)
{
    // loop over all surface elements
    for (int i=0; i<s.Elements(); ++i)
    {
        // get the surface element
        FESurfaceElement& el = s.Element(i);
        
        // calculate a penalty
        double eps = AutoViscousTractionPenalty(el, s);
        
        // assign to integration points of surface element
        int nint = el.GaussPoints();
        for (int j=0; j<nint; ++j)
        {
			FETiedFluidSurface::Data& pt = static_cast<FETiedFluidSurface::Data&>(*el.GetMaterialPoint(j));
			pt.m_epst = eps;
        }
    }
}

//-----------------------------------------------------------------------------
double FETiedFluidInterface::AutoViscousTractionPenalty(FESurfaceElement& el, FETiedFluidSurface& s)
{
    // get the solid element attached to the surface element
    FESolidElement& sel = static_cast<FESolidElement&>(*el.m_elem[0].pe);
    // get the fluid material for thast solid element
    FEMaterial* pmat = GetFEModel()->GetMaterial(sel.GetMatID());
    FEFluidMaterial* pfluid = dynamic_cast<FEFluidMaterial*>(pmat);
    if (pfluid == nullptr) return 0;
    // get the viscous part of this fluid
    FEViscousFluid* pvfluid = pfluid->GetViscous();
    
    // evaluate the viscosity for each material point ang get its average
    double eta = 0;
    int nint = sel.GaussPoints();
    for (int i=0; i<nint; ++i) {
        FEMaterialPoint* mp = sel.GetMaterialPoint(i);
        eta += pvfluid->ShearViscosity(*mp);
    }
    eta /= nint;
    
    // get the element thickness
    double area = s.GetArea(el);
    double vol = s.GetVolume(sel);
    double h = vol/area;
    
    return eta/h;
}

//-----------------------------------------------------------------------------
void FETiedFluidInterface::CalcAutoNormalVelocityPenalty(FETiedFluidSurface& s)
{
    // loop over all surface elements
    for (int i=0; i<s.Elements(); ++i)
    {
        // get the surface element
        FESurfaceElement& el = s.Element(i);
        
        // calculate a penalty
        double eps = AutoNormalVelocityPenalty(el, s);
        
        // assign to integration points of surface element
        int nint = el.GaussPoints();
        for (int j=0; j<nint; ++j)
        {
            FETiedFluidSurface::Data& pt = static_cast<FETiedFluidSurface::Data&>(*el.GetMaterialPoint(j));
            pt.m_epsn = eps;
        }
    }
}

//-----------------------------------------------------------------------------
double FETiedFluidInterface::AutoNormalVelocityPenalty(FESurfaceElement& el, FETiedFluidSurface& s)
{
    // get the solid element attached to the surface element
    FESolidElement& sel = static_cast<FESolidElement&>(*el.m_elem[0].pe);
    // get the fluid material for thast solid element
    FEMaterial* pmat = GetFEModel()->GetMaterial(sel.GetMatID());
    FEFluidMaterial* pfluid = dynamic_cast<FEFluidMaterial*>(pmat);
    if (pfluid == nullptr) return 0;
    // get the viscous part of this fluid
    FEViscousFluid* pvfluid = pfluid->GetViscous();
    
    // evaluate the viscosity and bulk modulus for each material point ang get its average
    double eta = 0;
    double k = 0;
    int nint = sel.GaussPoints();
    for (int i=0; i<nint; ++i) {
        FEMaterialPoint* mp = sel.GetMaterialPoint(i);
        eta += pvfluid->ShearViscosity(*mp);
    }
    eta /= nint;
    double tau = eta/pfluid->m_k;
    
    // get the element thickness
    double area = s.GetArea(el);
    double vol = s.GetVolume(sel);
    double h = vol/area;
    
    return h/tau;
}

//-----------------------------------------------------------------------------
// Perform initial projection between tied surfaces in reference configuration
void FETiedFluidInterface::InitialProjection(FETiedFluidSurface& s1, FETiedFluidSurface& s2)
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
            
			FETiedFluidSurface::Data& pt = static_cast<FETiedFluidSurface::Data&>(*el.GetMaterialPoint(j));
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
                
                // free the nodal dofs of pme if requested
                if (m_bfreedofs) {
                    for (int k=0; k<pme->Nodes(); ++k) {
                        FENode& node = mesh.Node(pme->m_node[k]);
                        FESurfaceElement* sme = cp.Project(node.m_rt, sq, srs);
                        if (sme) {
                            for (int l=0; l<m_dofWE.Size(); ++l)
                                if (node.get_bc(m_dofWE[l]) != DOF_OPEN) node.set_bc(m_dofWE[l], DOF_OPEN);
                        }
                    }
                }
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
// Evaluate gap functions for fluid velocity and fluid pressure
void FETiedFluidInterface::ProjectSurface(FETiedFluidSurface& s1, FETiedFluidSurface& s2)
{
    FEMesh& mesh = GetMesh();
    FESurfaceElement* pme;
    vec3d r;
    double alpha = GetFEModel()->GetTime().alphaf;
    
    vec3d  vt[FEElement::MAX_NODES], vp[FEElement::MAX_NODES], v1;
    double et[FEElement::MAX_NODES], ep[FEElement::MAX_NODES], e1;
    
    // loop over all integration points
    for (int i=0; i<s1.Elements(); ++i)
    {
        FESurfaceElement& el = s1.Element(i);
        
        int ne = el.Nodes();
        int nint = el.GaussPoints();
        
        // get the nodal velocities and dilatations
        for (int j=0; j<ne; ++j) {
            FENode& node = mesh.Node(el.m_node[j]);
            vt[j] = node.get_vec3d(m_dofWE[0], m_dofWE[1], m_dofWE[2]);
            vp[j] = node.get_vec3d_prev(m_dofWE[0], m_dofWE[1], m_dofWE[2]);
            et[j] = node.get(m_dofWE[3]);
            ep[j] = node.get_prev(m_dofWE[3]);
        }

        for (int j=0; j<nint; ++j)
        {
			FETiedFluidSurface::Data& pt = static_cast<FETiedFluidSurface::Data&>(*el.GetMaterialPoint(j));

            // calculate the global position of the integration point
            r = s1.Local2Global(el, j);
            
            // get the velocity and dilatation at the integration point
            v1 = el.eval(vt, j)*alpha + el.eval(vp, j)*(1-alpha);
            e1 = el.eval(et, j)*alpha + el.eval(ep, j)*(1-alpha);
            
            // if this node is tied, evaluate gap functions
            pme = pt.m_pme;
            if (pme)
            {
                // calculate the velocity gap function
                vec3d vmt[FEElement::MAX_NODES], vmp[FEElement::MAX_NODES];
                double emt[FEElement::MAX_NODES], emp[FEElement::MAX_NODES];
                for (int k=0; k<pme->Nodes(); ++k) {
                    FENode& node = mesh.Node(pme->m_node[k]);
                    vmt[k] = node.get_vec3d(m_dofWE[0], m_dofWE[1], m_dofWE[2]);
                    vmp[k] = node.get_vec3d_prev(m_dofWE[0], m_dofWE[1], m_dofWE[2]);
                    emt[k] = node.get(m_dofWE[3]);
                    emp[k] = node.get_prev(m_dofWE[3]);
                }
                vec3d v2 = pme->eval(vmt, pt.m_rs[0], pt.m_rs[1])*alpha + pme->eval(vmp, pt.m_rs[0], pt.m_rs[1])*(1-alpha);
                pt.m_vg = v2 - v1;

                double e2 = pme->eval(emt, pt.m_rs[0], pt.m_rs[1])*alpha + pme->eval(emp, pt.m_rs[0], pt.m_rs[1])*(1-alpha);
                pt.m_Jg = e1 - e2;
            }
            else
            {
                // the node is not tied
                pt.m_vg = vec3d(0,0,0);
                pt.m_Jg = 0;
            }
        }
    }
}

//-----------------------------------------------------------------------------

void FETiedFluidInterface::Update()
{
    // project the surfaces onto each other
    // this will update the gap functions as well
    ProjectSurface(m_s1, m_s2);
    if (m_btwo_pass) ProjectSurface(m_s2, m_s1);
    
}

//-----------------------------------------------------------------------------
void FETiedFluidInterface::LoadVector(FEGlobalVector& R, const FETimeInfo& tp)
{
    vector<int> LM1, LM2, LM, en;
    vector<double> fe;
    const int MI = FEElement::MAX_INTPOINTS;
    double detJ[MI], w[MI], *H1, H2[MI];
    const int MN = FEElement::MAX_NODES;
    vec3d f1[MN], f2[MN];
    double w1[MN], w2[MN];

    // loop over the nr of passes
    int npass = (m_btwo_pass?2:1);
    for (int np=0; np<npass; ++np)
    {
        // get primary and secondary surfaces
        FETiedFluidSurface& s1 = (np == 0? m_s1 : m_s2);
        FETiedFluidSurface& s2 = (np == 0? m_s2 : m_s1);
        
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
				FETiedFluidSurface::Data& pt = static_cast<FETiedFluidSurface::Data&>(*se1.GetMaterialPoint(j));

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
                    int ndpn = 4;
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
                    
                    // gap functions
                    vec3d dg = pt.m_vg;
                    double dJ = pt.m_Jg;
                    
                    // lagrange multipliers
                    vec3d Lmv = pt.m_Lmd;
                    double Lmp = pt.m_Lmp;
                    
                    // penalties
                    double epst = m_epst*pt.m_epst;
                    double epsn = m_epsn*pt.m_epsn;

                    // viscous traction
                    vec3d tv = Lmv + dg*epst;
                    pt.m_tv = tv;

                    // normal velocity jump
                    double vn = Lmp + dJ*epsn;
                    pt.m_vn = vn;
                    
                    // calculate the force vector
                    fe.resize(ndof);
                    zero(fe);
                    
                    for (int a=0; a<neln1; ++a) {
                        f1[a] = tv*H1[a];
                        w1[a] = vn*H1[a];
                    }
                    for (int b=0; b<neln2; ++b) {
                        f2[b] = -tv*H2[b];
                        w2[b] = -vn*H2[b];
                    }
                    
                    for (int a=0; a<neln1; ++a)
                    {
                        fe[4*a  ] += f1[a].x*detJ[j]*w[j];
                        fe[4*a+1] += f1[a].y*detJ[j]*w[j];
                        fe[4*a+2] += f1[a].z*detJ[j]*w[j];
                        fe[4*a+3] += w1[a]*detJ[j]*w[j];
                    }
                    for (int b = 0; b<neln2; ++b) {
                        fe[4*(b+neln1)  ] += f2[b].x*detJ[j]*w[j];
                        fe[4*(b+neln1)+1] += f2[b].y*detJ[j]*w[j];
                        fe[4*(b+neln1)+2] += f2[b].z*detJ[j]*w[j];
                        fe[4*(b+neln1)+3] += w2[b]*detJ[j]*w[j];
                    }

                    // assemble the global residual
                    R.Assemble(en, LM, fe);
                }
            }
        }
    }
}

//-----------------------------------------------------------------------------
void FETiedFluidInterface::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp)
{
    vector<int> LM1, LM2, LM, en;
    const int MI = FEElement::MAX_INTPOINTS;
    const int MN = FEElement::MAX_NODES;
    double detJ[MI], w[MI], *H1, H2[MI], pt[MN], dpr[MN], dps[MN];
    FEElementMatrix ke;
    
    double alpha = tp.alphaf;
    
    // do single- or two-pass
    int npass = (m_btwo_pass?2:1);
    for (int np=0; np < npass; ++np)
    {
		// get primary and secondary surfaces
		FETiedFluidSurface& s1 = (np == 0? m_s1 : m_s2);
        FETiedFluidSurface& s2 = (np == 0? m_s2 : m_s1);
        
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
            for (int j=0; j<nint1; ++j)
            {
                FETiedFluidSurface::Data& pt = static_cast<FETiedFluidSurface::Data&>(*se1.GetMaterialPoint(j));

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
                    ndpn = 4;
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
                    
                    // gap functions
                    vec3d dg = pt.m_vg;
                    double dJ = pt.m_Jg;
                    
                    // lagrange multipliers
                    vec3d Lmv = pt.m_Lmd;
                    double Lmp = pt.m_Lmp;
                    
                    // penalties
                    double epst = m_epst*pt.m_epst;
                    double epsn = m_epsn*pt.m_epsn;

                    // viscous traction
                    vec3d tv = Lmv + dg*epst;
                    pt.m_tv = tv;

                    // normal velocity jump
                    double vn = Lmp + dJ*epsn;
                    pt.m_vn = vn;
                                        
                    // create the stiffness matrix
                    ke.resize(ndof, ndof); ke.zero();
                    
                    //------------------------------------
                    
                    for (int a=0; a<neln1; ++a) {
                        for (int c=0; c<neln1; ++c)
                        {
                            mat3dd K11(-epst*H1[a]*H1[c]*detJ[j]*w[j]*alpha);
                            double k11 = epsn*H1[a]*H1[c]*detJ[j]*w[j]*alpha;
                            ke[ndpn*a    ][ndpn*c    ] -= K11.xx();
                            ke[ndpn*a + 1][ndpn*c + 1] -= K11.yy();
                            ke[ndpn*a + 2][ndpn*c + 2] -= K11.zz();
                            ke[ndpn*a + 3][ndpn*c + 3] -= k11;
                        }
                        for (int d=0; d<neln2; ++d)
                        {
                            mat3dd K12(epst*H1[a]*H2[d]*detJ[j]*w[j]*alpha);
                            double k12 = -epsn*H1[a]*H2[d]*detJ[j]*w[j]*alpha;
                            ke[ndpn*a    ][ndpn*(neln1+d)    ] -= K12.xx();
                            ke[ndpn*a + 1][ndpn*(neln1+d) + 1] -= K12.yy();
                            ke[ndpn*a + 2][ndpn*(neln1+d) + 2] -= K12.zz();
                            ke[ndpn*a + 3][ndpn*(neln1+d) + 3] -= k12;
                        }
                    }

                    for (int b=0; b<neln2; ++b) {
                        for (int c=0; c<neln1; ++c)
                        {
                            mat3dd K21(epst*H2[b]*H1[c]*detJ[j]*w[j]*alpha);
                            double k21 = -epsn*H2[b]*H1[c]*detJ[j]*w[j]*alpha;
                            ke[ndpn*(neln1+b)    ][ndpn*c    ] -= K21.xx();
                            ke[ndpn*(neln1+b) + 1][ndpn*c + 1] -= K21.yy();
                            ke[ndpn*(neln1+b) + 2][ndpn*c + 2] -= K21.zz();
                            ke[ndpn*(neln1+b) + 3][ndpn*c + 3] -= k21;
                        }
                        for (int d=0; d<neln2; ++d)
                        {
                            mat3dd K22(-epst*H2[b]*H2[d]*detJ[j]*w[j]*alpha);
                            double k22 = epsn*H2[b]*H2[d]*detJ[j]*w[j]*alpha;
                            ke[ndpn*(neln1+b)    ][ndpn*(neln1+d)    ] -= K22.xx();
                            ke[ndpn*(neln1+b) + 1][ndpn*(neln1+d) + 1] -= K22.yy();
                            ke[ndpn*(neln1+b) + 2][ndpn*(neln1+d) + 2] -= K22.zz();
                            ke[ndpn*(neln1+b) + 3][ndpn*(neln1+d) + 3] -= k22;
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
bool FETiedFluidInterface::Augment(int naug, const FETimeInfo& tp)
{
    // make sure we need to augment
	if (m_laugon != FECore::AUGLAG_METHOD) return true;

    int i;
    vec3d Ln;
    bool bconv = true;
    
    int N1 = m_s1.Elements();
    int N2 = m_s2.Elements();
    
    // --- c a l c u l a t e   i n i t i a l   n o r m s ---
    // a. normal component
    double normL0 = 0, normJ0 = 0;
    for (int i=0; i<N1; ++i)
    {
		FESurfaceElement& s1 = m_s1.Element(i);
        for (int j=0; j<s1.GaussPoints(); ++j)
        {
			FETiedFluidSurface::Data& d1 = static_cast<FETiedFluidSurface::Data&>(*s1.GetMaterialPoint(j));
			normL0 += d1.m_Lmd*d1.m_Lmd;
            normJ0 += d1.m_Lmp*d1.m_Lmp;
        }
    }
    for (int i=0; i<N2; ++i)
    {
		FESurfaceElement& s2 = m_s2.Element(i);
        for (int j=0; j<s2.GaussPoints(); ++j)
        {
			FETiedFluidSurface::Data& d2 = static_cast<FETiedFluidSurface::Data&>(*s2.GetMaterialPoint(j));
			normL0 += d2.m_Lmd*d2.m_Lmd;
            normJ0 += d2.m_Lmp*d2.m_Lmp;
        }
    }
    
    // b. gap component
    // (is calculated during update)
    double maxgap = 0;
    double maxJg = 0;
    
    // update Lagrange multipliers
    double normL1 = 0, normJ1 = 0, eps, epsn;
    for (i=0; i<N1; ++i)
    {
		FESurfaceElement& s1 = m_s1.Element(i);
		for (int j = 0; j<s1.GaussPoints(); ++j)
		{
			FETiedFluidSurface::Data& d1 = static_cast<FETiedFluidSurface::Data&>(*s1.GetMaterialPoint(j));

            if (d1.m_pme) {
                // update Lagrange multipliers on primary surface
                eps = m_epst*d1.m_epst;
                d1.m_Lmd = d1.m_Lmd + d1.m_vg*eps;
                maxgap = max(maxgap,sqrt(d1.m_vg*d1.m_vg));
                normL1 += d1.m_Lmd*d1.m_Lmd;
                
                epsn = m_epsn*d1.m_epsn;
                d1.m_Lmp = d1.m_Lmp + epsn*d1.m_Jg;
                maxJg = max(maxJg,fabs(d1.m_Jg));
                normJ1 += d1.m_Lmp*d1.m_Lmp;
            }
        }
    }
    
    for (i=0; i<N2; ++i)
    {
		FESurfaceElement& s2 = m_s2.Element(i);
		for (int j = 0; j<s2.GaussPoints(); ++j)
		{
			FETiedFluidSurface::Data& d2 = static_cast<FETiedFluidSurface::Data&>(*s2.GetMaterialPoint(j));

            if (d2.m_pme) {
                // update Lagrange multipliers on secondary surface
                eps = m_epst*d2.m_epst;
                d2.m_Lmd = d2.m_Lmd + d2.m_vg*eps;
                maxgap = max(maxgap,sqrt(d2.m_vg*d2.m_vg));
                normL1 += d2.m_Lmd*d2.m_Lmd;
                
                epsn = m_epsn*d2.m_epsn;
                d2.m_Lmp = d2.m_Lmp + epsn*d2.m_Jg;
                double maxJg = max(maxJg,fabs(d2.m_Jg));
                normJ1 += d2.m_Lmp*d2.m_Lmp;
            }
        }
    }
    
    // calculate relative norms
    double lnorm = (normL1 != 0 ? fabs((normL1 - normL0) / normL1) : fabs(normL1 - normL0));
    double pnorm = (normJ1 != 0 ? fabs((normJ1 - normJ0) / normJ1) : fabs(normJ1 - normJ0));
    
    // check convergence
    if ((m_gtol > 0) && (maxgap > m_gtol)) bconv = false;
    if ((m_etol > 0) && (maxJg > m_etol)) bconv = false;
    
    if ((m_atol > 0) && (lnorm > m_atol)) bconv = false;
    if ((m_atol > 0) && (pnorm > m_atol)) bconv = false;
    
    if (naug < m_naugmin ) bconv = false;
    if (naug >= m_naugmax) bconv = true;
    
    feLog(" tied fluid interface # %d\n", GetID());
    feLog("                        CURRENT        REQUIRED\n");
    feLog("    V multiplier : %15le", lnorm); if (m_atol > 0) feLog("%15le\n", m_atol); else feLog("       ***\n");
    feLog("    P multiplier        : %15le", pnorm); if (m_atol > 0) feLog("%15le\n", m_atol); else feLog("       ***\n");
    
    feLog("    maximum velocity gap  : %15le", maxgap);
    if (m_gtol > 0) feLog("%15le\n", m_gtol); else feLog("       ***\n");
    feLog("    maximum pressure gap : %15le", maxJg);
    if (m_etol > 0) feLog("%15le\n", m_etol); else feLog("       ***\n");

    return bconv;
}

//-----------------------------------------------------------------------------
void FETiedFluidInterface::Serialize(DumpStream &ar)
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
    ar & m_pfluid;
    ar & m_dofWE;
}
