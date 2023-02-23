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
#include "FESoluteNaturalFlux.h"
#include "FEMultiphasic.h"
#include <FECore/FEModel.h>
#include <FECore/FEAnalysis.h>
#include <FECore/FESolidDomain.h>

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FESoluteNaturalFlux, FESurfaceLoad)
    ADD_PARAMETER(m_bshellb, "shell_bottom");
    ADD_PARAMETER(m_isol   , "solute_id")->setEnums("$(solutes)");
    ADD_PARAMETER(m_bup    , "update");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! constructor
FESoluteNaturalFlux::FESoluteNaturalFlux(FEModel* pfem) : FESurfaceLoad(pfem), m_dofC(pfem), m_dofU(pfem), m_dofP(pfem)
{
    m_bshellb = false;
    m_isol = -1;
    m_bup = false;
}
    
//-----------------------------------------------------------------------------
//! allocate storage
void FESoluteNaturalFlux::SetSurface(FESurface* ps)
{
    FESurfaceLoad::SetSurface(ps);
}

//-----------------------------------------------------------------------------
//! serialization
void FESoluteNaturalFlux::Serialize(DumpStream& ar)
{
    FESurfaceLoad::Serialize(ar);

    if (ar.IsShallow() == false)
    {
        ar & m_dofC & m_dofU & m_dofP;
    }
}

//-----------------------------------------------------------------------------
bool FESoluteNaturalFlux::Init()
{
    if (m_isol <= 0) return false;

    // set up the dof lists
    FEModel* fem = GetFEModel();
    m_dofC.Clear();
    m_dofU.Clear();
    m_dofP.Clear();
    if (m_bshellb == false)
    {
        m_dofC.AddDof(fem->GetDOFIndex("concentration", m_isol - 1));

        m_dofU.AddDof(fem->GetDOFIndex("x"));
        m_dofU.AddDof(fem->GetDOFIndex("y"));
        m_dofU.AddDof(fem->GetDOFIndex("z"));
        
        m_dofP.AddDof(fem->GetDOFIndex("p"));
    }
    else
    {
        m_dofC.AddDof(fem->GetDOFIndex("shell concentration", m_isol - 1));

        m_dofU.AddDof(fem->GetDOFIndex("sx"));
        m_dofU.AddDof(fem->GetDOFIndex("sy"));
        m_dofU.AddDof(fem->GetDOFIndex("sz"));

        m_dofP.AddDof(fem->GetDOFIndex("q"));
    }
    m_dof.AddDofs(m_dofU);
    m_dof.AddDofs(m_dofP);
    m_dof.AddDofs(m_dofC);
    return FESurfaceLoad::Init();
}

//-----------------------------------------------------------------------------
void FESoluteNaturalFlux::Update()
{
    if (m_bup) {
        for (int is=0; is<m_psurf->Elements(); ++is)
        {
            // get surface element
            FESurfaceElement& el = m_psurf->Element(is);
            // get underlying solid element
            FESolidElement* pe = dynamic_cast<FESolidElement*>(el.m_elem[0]);
            if (pe == nullptr) break;
            // get element data
            int neln = pe->Nodes();
            int nint = pe->GaussPoints();
            
            FEMaterial* pm = GetFEModel()->GetMaterial(pe->GetMatID());
            // get the local solute id
            FESoluteInterface* psi = dynamic_cast<FESoluteInterface*>(pm);
            if (psi == nullptr) break;
            int sid = psi->FindLocalSoluteID(m_isol);
            if (sid == -1) break;
            
            // identify nodes on the surface
            vector<bool> nsrf(neln,false);
            for (int j=0; j<neln; ++j) {
                for (int k=0; k < el.Nodes(); ++k) {
                    if (el.m_node[k] == pe->m_node[j]) nsrf[j] = true;
                }
            }
            
            // get average effective concentration of nodes not on surface
            double cavg = 0;
            int m = 0;
            for (int i=0; i<neln; ++i) {
                if (!nsrf[i]) {
                    int n = pe->m_node[i];
                    FENode& node = GetMesh().Node(n);
                    int dof = m_dofC[m_isol-1];
                    if (dof != -1) {
                        cavg += node.get(dof);
                        ++m;
                    }
                }
            }
            // assign this average value to surface nodes as initial guess
            if (m) {
                cavg /= m;
                for (int i=0; i<neln; ++i) {
                    if (nsrf[i]) {
                        int n = pe->m_node[i];
                        FENode& node = GetMesh().Node(n);
                        int dof = m_dofC[m_isol-1];
                        if (dof != -1) node.set(dof, cavg);
                    }
                }
            }
        }
    }
}

//-----------------------------------------------------------------------------
void FESoluteNaturalFlux::LoadVector(FEGlobalVector& R)
{
    double dt = CurrentTimeIncrement();

    // element force vector
    vector<double> fe;
    vector<int> lm;
    // jacobian matrix, inverse jacobian matrix and determinants
    double Ji[3][3], detJt;
    const double* Gr, *Gs, *Gt, *H;

    for (int is=0; is<m_psurf->Elements(); ++is)
    {
        // get surface element
        FESurfaceElement& el = m_psurf->Element(is);
        // get surface normal
        vec3d nu(0,0,0);
        for (int n=0; n<el.GaussPoints(); ++n) {
            FESurfaceMaterialPoint* pt = dynamic_cast<FESurfaceMaterialPoint*>(el.GetMaterialPoint(n));
            nu += pt->dxr ^ pt->dxs;
        }
        nu.unit();
        // get underlying solid element
        FESolidElement* pe = dynamic_cast<FESolidElement*>(el.m_elem[0]);
        if (pe == nullptr) break;
        // determine the solid domain to which this solid element belongs
        FESolidDomain* sdom = dynamic_cast<FESolidDomain*>(pe->GetMeshPartition());
        // get element data
        int nint = pe->GaussPoints();
        int neln = pe->Nodes();
        double* gw = pe->GaussWeights();
        
        FEMaterial* pm = GetFEModel()->GetMaterial(pe->GetMatID());
        // get the local solute id
        FESoluteInterface* psi = dynamic_cast<FESoluteInterface*>(pm);
        if (psi == nullptr) break;
        int sid = psi->FindLocalSoluteID(m_isol);
        if (sid == -1) break;
        
        // get the element force vector and initialize it to zero
        fe.assign(neln, 0);    // 1 concentration dof per node
        lm.resize(neln);
        // unpack lm and get nodal effective solute concentrations
        vector<double> ce(neln,0);
        for (int i=0; i<neln; ++i) {
            int n = pe->m_node[i];
            FENode& node = GetMesh().Node(n);
            vector<int>& id = node.m_ID;
            int dof = m_dofC[m_isol-1];
            if (dof != -1) {
                lm[i] = id[dof];
                ce[i] = node.get(dof);
            }
        }

        // for each integration point in the solid element
        for (int n=0; n<nint; ++n) {
            FEMaterialPoint& pt = *pe->GetMaterialPoint(n);
            FEBiphasicMaterialPoint& pb = *(pt.ExtractData<FEBiphasicMaterialPoint>());
            FESolutesMaterialPoint& ps = *(pt.ExtractData<FESolutesMaterialPoint>());

            // calculate the jacobian
            detJt = sdom->invjact(*pe, Ji, n);
            detJt *= gw[n]*dt;
            // get shape functions and their derivatives
            H = pe->H(n);
            Gr = pe->Gr(n);
            Gs = pe->Gs(n);
            Gt = pe->Gt(n);
            
            // get contravariant basis vectors
            vec3d gcntv[3];
            sdom->ContraBaseVectors(*pe, n, gcntv);
            
            // evaluate gradient of shape function and gradient of effective concentration
            // (using ps.m_gradc[n] doesn't work, because it doesn't get updated until convergence)
            vector<vec3d> gradN(neln);
            vec3d gradc(0,0,0);
            for (int i=0; i<neln; ++i) {
                gradN[i] = gcntv[0]*Gr[i] + gcntv[1]*Gs[i] + gcntv[2]*Gt[i];
                gradc += gradN[i]*ce[i];
            }

            for (int i=0; i<neln; ++i)
                fe[i] -= H[i]*(gradc*nu)*detJt;
        }
        
        R.Assemble(pe->m_node, lm, fe);
    }

    // Now do the surface implementation, the normal way
    FESoluteNaturalFlux* flux = this;
    m_psurf->LoadVector(R, m_dofC, false, [=](FESurfaceMaterialPoint& mp, const FESurfaceDofShape& dof_a, std::vector<double>& fa) {

        // get surface element
        FESurfaceElement& el = *mp.SurfaceElement();
        // get underlying solid element
        FEElement* pe = el.m_elem[0];
        FEMaterial* pm = GetFEModel()->GetMaterial(pe->GetMatID());
        // get the local solute id
        FESoluteInterface* psi = dynamic_cast<FESoluteInterface*>(pm);
        if (psi == nullptr) {
            fa[0] = 0;
            return;
        }
        int sid = psi->FindLocalSoluteID(flux->m_isol);
        
        // get element-averaged fluid flux and actual solute concentration
        vec3d w(0,0,0);
        double c = 0;
        int nint = pe->GaussPoints();
        for (int n=0; n<nint; ++n) {
            FEMaterialPoint& pt = *pe->GetMaterialPoint(n);
            FEBiphasicMaterialPoint& pb = *(pt.ExtractData<FEBiphasicMaterialPoint>());
            FESolutesMaterialPoint& ps = *(pt.ExtractData<FESolutesMaterialPoint>());
            w += pb.m_w;
            c += ps.m_ca[sid];
        }
        w /= nint;
        c /= nint;

        // evaluate desired natural solute flux
        vec3d dxt = mp.dxr ^ mp.dxs;
        double jn = c*(w*dxt);
        if (flux->m_bshellb) jn = -jn;

        // molar flow rate
        double f = jn* dt;

        double H_i = dof_a.shape;
        fa[0] = H_i * f;
    });
}

//-----------------------------------------------------------------------------
void FESoluteNaturalFlux::StiffnessMatrix(FELinearSystem& LS)
{
    // time increment
    double dt = CurrentTimeIncrement();

    int ndpn = 4;   // 3 displacement dofs + 1 concentration dof
    // element stiffness matrix
    vector<int> lm;
    // jacobian matrix, inverse jacobian matrix and determinants
    double Ji[3][3], detJt;
    const double* Gr, *Gs, *Gt, *H;
    
    for (int is=0; is<m_psurf->Elements(); ++is)
    {
        // get surface element
        FESurfaceElement& el = m_psurf->Element(is);
        // get surface normal
        vec3d nu(0,0,0);
        for (int n=0; n<el.GaussPoints(); ++n) {
            FESurfaceMaterialPoint* pt = dynamic_cast<FESurfaceMaterialPoint*>(el.GetMaterialPoint(n));
            nu += pt->dxr ^ pt->dxs;
        }
        nu.unit();
        // get underlying solid element
        FESolidElement* pe = dynamic_cast<FESolidElement*>(el.m_elem[0]);
        if (pe == nullptr) break;
        // determine the solid domain to which this solid element belongs
        FESolidDomain* sdom = dynamic_cast<FESolidDomain*>(pe->GetMeshPartition());
        // get element data
        int nint = pe->GaussPoints();
        int neln = pe->Nodes();
        double* gw = pe->GaussWeights();
        int ndof = neln*ndpn;
        FEElementMatrix ke(*pe);
        ke.resize(ndof, ndof);

        FEMaterial* pm = GetFEModel()->GetMaterial(pe->GetMatID());
        // get the local solute id
        FESoluteInterface* psi = dynamic_cast<FESoluteInterface*>(pm);
        if (psi == nullptr) break;
        int sid = psi->FindLocalSoluteID(m_isol);
        if (sid == -1) break;
        
        // initialize stiffness matrix it to zero
        ke.zero();
        lm.resize(ndof);
        // unpack lm and get nodal effective solute concentrations
        vector<double> ce(neln,0);
        for (int i=0; i<neln; ++i) {
            int n = pe->m_node[i];
            FENode& node = GetMesh().Node(n);
            vector<int>& id = node.m_ID;
            lm[ndpn*i  ] = id[m_dofU[0]];
            lm[ndpn*i+1] = id[m_dofU[1]];
            lm[ndpn*i+2] = id[m_dofU[2]];
            int dof = m_dofC[m_isol-1];
            if (dof != -1) {
                lm[ndpn*i+3] = id[dof];
                ce[i] = node.get(dof);
            }
        }
        ke.SetIndices(lm);
        
        // for each integration point in the solid element
        for (int n=0; n<nint; ++n) {
            FEMaterialPoint& pt = *pe->GetMaterialPoint(n);
            FEBiphasicMaterialPoint& pb = *(pt.ExtractData<FEBiphasicMaterialPoint>());
            FESolutesMaterialPoint& ps = *(pt.ExtractData<FESolutesMaterialPoint>());
            
            // calculate the jacobian
            detJt = sdom->invjact(*pe, Ji, n);
            detJt *= gw[n]*dt;
            // get shape functions and their derivatives
            H = pe->H(n);
            Gr = pe->Gr(n);
            Gs = pe->Gs(n);
            Gt = pe->Gt(n);
            
            // get contravariant basis vectors
            vec3d gcntv[3];
            sdom->ContraBaseVectors(*pe, n, gcntv);
            
            // evaluate gradient of shape function and gradient of effective concentration
            // don't use ps.m_gradc[n] as it doesn't get updated until next convergence
            vector<vec3d> gradN(neln);
            vec3d gradc(0,0,0);
            for (int i=0; i<neln; ++i) {
                gradN[i] = gcntv[0]*Gr[i] + gcntv[1]*Gs[i] + gcntv[2]*Gt[i];
                gradc += gradN[i]*ce[i];
            }
            for (int i=0, in = 0; i<neln; ++i, in += ndpn) {
                for (int j=0, jn = 0; j<neln; ++j, jn += ndpn) {
                    vec3d kcu = (gradN[j]*(nu*gradc) - gradc*(gradN[j]*nu))*H[i];
                    double kcc = H[i]*(gradN[j]*nu);
                    
                    ke[in+3][jn  ] += kcu.x*detJt;
                    ke[in+3][jn+1] += kcu.y*detJt;
                    ke[in+3][jn+2] += kcu.z*detJt;
                    ke[in+3][jn+3] += kcc*detJt;
                }
            }
        }
        
        LS.Assemble(ke);
    }

    // Now do the surface implementation, the normal way
    // evaluate the stiffness contribution
    FESoluteNaturalFlux* flux = this;
    m_psurf->LoadStiffness(LS, m_dofC, m_dofU, [=](FESurfaceMaterialPoint& mp, const FESurfaceDofShape& dof_a, const FESurfaceDofShape& dof_b, matrix& Kab) {

        // get surface element
        FESurfaceElement& el = *mp.SurfaceElement();
        // get underlying solid element
        FEElement* pe = el.m_elem[0];
        FEMaterial* pm = GetFEModel()->GetMaterial(pe->GetMatID());
        // get the local solute id
        FESoluteInterface* psi = dynamic_cast<FESoluteInterface*>(pm);
        if (psi == nullptr) return;
        int sid = psi->FindLocalSoluteID(flux->m_isol);
        
        // get element-averaged fluid flux and actual solute concentration
        vec3d w(0,0,0);
        double c = 0;
        int nint = pe->GaussPoints();
        for (int n=0; n<nint; ++n) {
            FEMaterialPoint& pt = *pe->GetMaterialPoint(n);
            FEBiphasicMaterialPoint& pb = *(pt.ExtractData<FEBiphasicMaterialPoint>());
            FESolutesMaterialPoint& ps = *(pt.ExtractData<FESolutesMaterialPoint>());
            w += pb.m_w;
            c += ps.m_ca[sid];
        }
        w /= nint;
        c /= nint;
        
        // shape functions and derivatives
        double H_i  = dof_a.shape;
        double Gr_j = dof_b.shape_deriv_r;
        double Gs_j = dof_b.shape_deriv_s;

        // calculate surface normal
        vec3d dxt = mp.dxr ^ mp.dxs;
        vec3d nu = dxt.normalized();
        double jn = c*(w*nu);
        if (flux->m_bshellb) jn = -jn;

        // calculate stiffness component
        vec3d t1 = nu*jn;
        vec3d t2 = mp.dxs*Gr_j - mp.dxr*Gs_j;
        vec3d kab = (t1 ^ t2)*(H_i)*dt;

        Kab[0][0] = kab.x;
        Kab[0][1] = kab.y;
        Kab[0][2] = kab.z;
    });
}
