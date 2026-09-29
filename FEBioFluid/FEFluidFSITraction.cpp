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
#include "FEFluidFSITraction.h"
#include <FECore/FEModel.h>
#include <FECore/log.h>
#include <FECore/FELinearSystem.h>
#include "FEFluid.h"
#include "FEFluidFSI.h"
#include "FEBioFSI.h"

//-----------------------------------------------------------------------------
// Parameter block for pressure loads
BEGIN_FECORE_CLASS(FEFluidFSITraction, FESurfaceLoad)
    ADD_PARAMETER(m_bshellb , "shell_bottom");
    ADD_PARAMETER(m_btied   , "use_tied_elastic_interface");
END_FECORE_CLASS()

//-----------------------------------------------------------------------------
//! constructor
FEFluidFSITraction::FEFluidFSITraction(FEModel* pfem) : FESurfaceLoad(pfem), m_dofU(pfem), m_dofSU(pfem), m_dofW(pfem)
{
	m_bshellb = false;
    m_btied = false;
    m_tei = nullptr;
    m_psolid = nullptr;

    // get the degrees of freedom
	// TODO: Can this be done in Init, since  there is no error checking
	if (pfem)
	{
		m_dofU.AddVariable(FEBioFSI::GetVariableName(FEBioFSI::DISPLACEMENT));
		m_dofSU.AddVariable(FEBioFSI::GetVariableName(FEBioFSI::SHELL_DISPLACEMENT));
		m_dofW.AddVariable(FEBioFSI::GetVariableName(FEBioFSI::RELATIVE_FLUID_VELOCITY));
		m_dofEF = GetDOFIndex(FEBioFSI::GetVariableName(FEBioFSI::FLUID_DILATATION), 0);

		m_dof.Clear();
		m_dof.AddDofs(m_dofU);
		m_dof.AddDofs(m_dofSU);
		m_dof.AddDofs(m_dofW);
		m_dof.AddDof(m_dofEF);
	}
}

//-----------------------------------------------------------------------------
//! initialize
bool FEFluidFSITraction::Init()
{
    // for now, let's pick the first available tied-elastic interface
    // TODO: Allow user to select the tied-elastic interface
    if (m_btied) {
        FEModel& fem = *GetFEModel();
        // pick first tied-elastic interface
        if (fem.SurfacePairConstraints() > 0)
        {
            // loop over all contact interfaces
            for (int i = 0; i<fem.SurfacePairConstraints(); ++i)
            {
                FEContactInterface* pci = dynamic_cast<FEContactInterface*>(fem.SurfacePairConstraint(i));
                FETiedElasticInterface* pbw = dynamic_cast<FETiedElasticInterface*>(pci);
                if (pbw) {
                    m_tei = pbw;
                    break;
                }
            }
        }
        // war used that there is no tied-elastic interface in this model
        if (m_tei == nullptr) {
            feLogError("This model must include a tied-elastic interface!");
            return false;
        }
        // initialize this interface
        if (m_tei->Init() == false) return false;

        // determine which of the two surfaces in this tied-elastic interface is the fluid-FSI domain
        FESurface& s1 = *(m_tei->GetPrimarySurface());
        FESurface& s2 = *(m_tei->GetSecondarySurface());
        //! Pick the first face associated with either of these surfaces
        FESurfaceElement& el1 = s1.Element(0);
        FESurfaceElement& el2 = s2.Element(0);
        // extract the first solid element on each surface
        FEElement* sel1 = el1.m_elem[0].pe;
        FEElement* sel2 = el2.m_elem[0].pe;
        // get the material of the first solid element in s1 and associated FluidFSI
        FEMaterial* pm1 = GetFEModel()->GetMaterial(sel1->GetMatID());
        FEFluidFSI* pfsi1 = dynamic_cast<FEFluidFSI*>(pm1);
        // get the material of the first solid element in s2 and associated FluidFSI
        FEMaterial* pm2 = GetFEModel()->GetMaterial(sel2->GetMatID());
        FEFluidFSI* pfsi2 = dynamic_cast<FEFluidFSI*>(pm2);
        // check which of these two is a FluidFSI material
        if (pfsi1) {
            SetSurface(m_tei->GetPrimarySurface());
            m_psolid = m_tei->GetSecondarySurface();
        }
        else if (pfsi2) {
            SetSurface(m_tei->GetSecondarySurface());
            m_psolid = m_tei->GetPrimarySurface();
        }
        else {
            feLogError("Neither of the surfaces in the tied-elastic interface belongs to a fluid-FSI domain!");
            return false;
        }
    }

    // TODO: Deal with the case when the surface is a shell domain separating two FSI domains
    // that use different fluid bulk moduli
    // (for now, users have to define two FEFluidFSITraction loads, one on front shell
    // face and the other on back shell face)
    
    FESurface& surf = GetSurface();
    surf.SetShellBottom(m_bshellb);
    surf.SetInterfaceStatus(true);
    if (FESurfaceLoad::Init() == false) return false;
    
    return true;
}

void FEFluidFSITraction::Activate()
{
    // proceed with usual initialization
    FESurface& surf = GetSurface();
    
    // get the list of fluid-FSI elements connected to this interface
    FEModel* fem = GetFEModel();
    int NF = surf.Elements();
    m_elem.resize(NF);
    m_s.resize(NF, 1);
    for (int j = 0; j < NF; ++j)
    {
        bool bself = false;
        FESurfaceElement& el = surf.Element(j);
        // extract the first of two elements on this interface
        m_elem[j] = el.m_elem[0].pe;
        if (el.m_elem[1].pe == nullptr) bself = true;
        // get its material and check if FluidFSI
        FEMaterial* pm = fem->GetMaterial(m_elem[j]->GetMatID());
        FEFluidFSI* pfsi = dynamic_cast<FEFluidFSI*>(pm);
        if (pfsi) {
            double s = m_psurf->FacePointing(el, *m_elem[j]);
            m_s[j] = bself ? -s : s;
            assert(m_s[j]);
        }
        else if (!bself) {
            // extract the second of two elements on this interface
            m_elem[j] = el.m_elem[1].pe;
            pm = fem->GetMaterial(m_elem[j]->GetMatID());
            pfsi = dynamic_cast<FEFluidFSI*>(pm);
            assert(pfsi);
            m_s[j] = m_psurf->FacePointing(el, *m_elem[j]);
            assert(m_s[j]);
        }
        else
            assert(false);
    }
}

//-----------------------------------------------------------------------------
double FEFluidFSITraction::GetFluidDilatation(FESurfaceMaterialPoint& mp, double alpha)
{
	double ef = 0;
	FESurfaceElement& el = *mp.SurfaceElement();
	double* H = el.H(mp.m_index);
	int neln = el.Nodes();
	for (int j = 0; j < neln; ++j) {
		FENode& node = m_psurf->Node(el.m_lnode[j]);
		double ej = node.get(m_dofEF)*alpha + node.get_prev(m_dofEF)*(1.0 - alpha);
		ef += ej*H[j];
	}
	return ef;
}

//-----------------------------------------------------------------------------
mat3ds FEFluidFSITraction::GetFluidStress(FESurfaceMaterialPoint& pt)
{
	FEModel* fem = GetFEModel();
	FESurfaceElement& face = *pt.SurfaceElement();
	int iel = face.m_lid;

	// Get the fluid stress from the fluid-FSI element
	mat3ds sv(mat3dd(0));
	FEElement* pe = m_elem[iel];
	int nint = pe->GaussPoints();
	FEFluidFSI* pfsi = dynamic_cast<FEFluidFSI*>(fem->GetMaterial(pe->GetMatID()));
	for (int n = 0; n<nint; ++n)
	{
		FEMaterialPoint& mp = *pe->GetMaterialPoint(n);
		sv += pfsi->Fluid()->GetViscous()->Stress(mp);
	}
	sv /= nint;
	return sv;
}

//-----------------------------------------------------------------------------
void FEFluidFSITraction::LoadVector(FEGlobalVector& R)
{
	const FETimeInfo& tp = GetTimeInfo();

    if (m_btied == false) {
        // If surface is bottom of shell, we should take shell displacement dofs (i.e. m_dofSU).
        FEDofList dof = m_bshellb ? m_dofSU : m_dofU;
        m_psurf->LoadVector(R, dof, false, [&](FESurfaceMaterialPoint& mp, const FESurfaceDofShape& dof_a, vector<double>& fa) {
            
            // get the surface element
            FESurfaceElement& el = *mp.SurfaceElement();
            int iel = el.m_lid;
            FEFluidFSI* pfsi = dynamic_cast<FEFluidFSI*>(GetFEModel()->GetMaterial(m_elem[iel]->GetMatID()));
            
            // nodal coordinates
            vec3d rt[FEElement::MAX_NODES];
            m_psurf->GetNodalCoordinates(el, tp.alphaf, rt);
            
            // evaluate covariant basis vectors at integration point
            vec3d gr = el.eval_deriv1(rt, mp.m_index)*m_s[iel];
            vec3d gs = el.eval_deriv2(rt, mp.m_index);
            vec3d gt = gr ^ gs;
            
            // Get the fluid viscous stress at integration point
            // necessarily using the attached solid element
            mat3ds sv = GetFluidStress(mp);
            
            // fluid dilatation at integration point
            // only from surface element
            double ef = GetFluidDilatation(mp, tp.alphaf);
            double p = pfsi->Fluid()->Pressure(ef);
            
            // evaluate traction
            vec3d f = gt*p - sv*gt;
            
            double H = dof_a.shape;
            fa[0] = H * f.x;
            fa[1] = H * f.y;
            fa[2] = H * f.z;
        });
    }
    else {
        FETimeInfo& tp = GetFEModel()->GetTime();
        for (int i=0; i<m_psurf->Elements(); ++i) {
            // get the surface element
            FESurfaceElement& elfsi = m_psurf->Element(i);
            int ielfsi = elfsi.m_lid;
            FEFluidFSI* pfsi = dynamic_cast<FEFluidFSI*>(GetFEModel()->GetMaterial(m_elem[ielfsi]->GetMatID()));

            // find number of integration points on fluid-FSI surface
            int nifsi = elfsi.GaussPoints();
            // loop over integration points of fluid-FSI surface
            for (int k=0; k<nifsi; ++k) {
                // get the fluid-FSI surface material point
                FESurfaceMaterialPoint& mp = static_cast<FESurfaceMaterialPoint&>(*elfsi.GetMaterialPoint(k));
                
                // extract fluid-FSI data
                FETiedElasticSurface::Data& data = static_cast<FETiedElasticSurface::Data&>(*elfsi.GetMaterialPoint(k));
                // get surface element of solid surface
                FESurfaceElement* pe = data.m_pme;
                // if it exists, calculate the fluid force acting at the ray intersection points
                if (pe != 0)
                {
                    FESurfaceElement& me = *pe;
                    vector<double> fe;
                    vector<int> lm;
                    int ndpn = 3;
                    int ndof = ndpn*me.Nodes();
                    fe.assign(ndof, 0);
                    lm.assign(ndof,0);
                    // evaluate lm
                    for (int i=0; i<pe->Nodes(); ++i) {
                        FENode& node = GetFEModel()->GetMesh().Node(me.m_node[i]);
                        vector<int>& id = node.m_ID;
                        // displacement dofs
                        for (int j=0; j<ndpn; ++j) lm[ndpn*i+j] = id[m_dofU[j]];
                    }
                    
                    // get nodal coordinates on solid face
                    vec3d rt[FEElement::MAX_NODES];
                    m_psolid->GetNodalCoordinates(me, tp.alphaf, rt);
                    
                    // evaluate covariant basis vectors at the ray intersection with the solid face
                    vec3d gr = me.eval_deriv1(rt, data.m_rs[0], data.m_rs[1]);
                    vec3d gs = me.eval_deriv2(rt, data.m_rs[0], data.m_rs[1]);
                    vec3d gt = gr ^ gs;
                    
                    // Get the fluid viscous stress at fluid-FSI integration point
                    mat3ds sv = GetFluidStress(mp);
                    
                    // fluid dilatation at fluid-FSI integration point
                    double ef = GetFluidDilatation(mp, tp.alphaf);
                    double p = pfsi->Fluid()->Pressure(ef);
                    
                    // evaluate traction on the solid due to the fluid in the fluid-FSI domain
                    vec3d f = -gt*p + sv*gt;
                    
                    // get the shape functions on the solid face
                    double *H;
                    me.shape_fnc(H,data.m_rs[0], data.m_rs[1]);
                    
                    for (int n=0; n<me.Nodes(); ++n) {
                        fe[ndpn*n  ] -= H[n] * f.x;
                        fe[ndpn*n+1] -= H[n] * f.y;
                        fe[ndpn*n+2] -= H[n] * f.z;
                    }

                    // assemble this element load vector into the global load vector
                    R.Assemble(me.m_node, lm, fe);
                }
            }
        }
    }
}

//-----------------------------------------------------------------------------
void FEFluidFSITraction::StiffnessMatrix(FELinearSystem& LS)
{
	FEModel* fem = GetFEModel();
    const FETimeInfo& tp = GetTimeInfo();
    double dt = tp.timeIncrement;
    double alpha = tp.alphaf;
    double a = tp.gamma / (tp.beta*dt);

        FESurface* ps = &GetSurface();
        
        
        // build dof list
        // TODO: If surface is bottom of shell, we should take shell displacement dofs (i.e. m_dofSU).
        FEDofList dofs(fem);
        if (!m_bshellb) dofs.AddDofs(m_dofU); else dofs.AddDofs(m_dofSU);
        dofs.AddDofs(m_dofW);
        dofs.AddDof(m_dofEF);
        
    // evaluate stiffness
    m_psurf->LoadStiffness(LS, dofs, dofs, [&](FESurfaceMaterialPoint& mp, const FESurfaceDofShape& dof_a, const FESurfaceDofShape& dof_b, matrix& Kab) {
        
        FESurfaceElement& el = *mp.SurfaceElement();
        int iel = el.m_lid;
        int neln = el.Nodes();
        
        vector<vec3d> gradN(neln);
        
        // nodal coordinates
        vec3d rt[FEElement::MAX_NODES];
        ps->GetNodalCoordinates(el, tp.alphaf, rt);
        
        // Get the fluid stress and its tangents from the fluid-FSI element
        mat3ds sv(mat3dd(0)), svJ(mat3dd(0));
        tens4ds cv; cv.zero();
        mat3d Ls; Ls.zero();
        FEElement* pe = m_elem[iel];
        int pint = pe->GaussPoints();
        FEFluidFSI* pfsi = dynamic_cast<FEFluidFSI*>(fem->GetMaterial(pe->GetMatID()));
        for (int n = 0; n<pint; ++n)
        {
            FEMaterialPoint& mp = *pe->GetMaterialPoint(n);
            FEElasticMaterialPoint& ep = *(mp.ExtractData<FEElasticMaterialPoint>());
            sv += pfsi->Fluid()->GetViscous()->Stress(mp);
            svJ += pfsi->Fluid()->GetViscous()->Tangent_Strain(mp);
            cv += pfsi->Fluid()->Tangent_RateOfDeformation(mp);
            Ls += ep.m_L;
        }
        sv /= pint;
        svJ /= pint;
        cv /= pint;
        Ls /= pint;
        mat3d M = mat3dd(a) - Ls;
        
        double* N  = el.H (mp.m_index);
        double* Gr = el.Gr(mp.m_index);
        double* Gs = el.Gs(mp.m_index);
        
        // evaluate fluid dilatation
        double ef = GetFluidDilatation(mp, tp.alphaf);
        
        // covariant basis vectors
        vec3d gr = el.eval_deriv1(rt, mp.m_index)*m_s[iel];
        vec3d gs = el.eval_deriv2(rt, mp.m_index);
        vec3d gt = gr ^ gs;
        
        // evaluate fluid pressure
        double p = pfsi->Fluid()->Pressure(ef);
        
        vec3d f = gt*pfsi->Fluid()->GetElastic()->Tangent_Strain(ef,0);
        
        vec3d gcnt[2], gcntp[2];
        ps->ContraBaseVectors(el, mp.m_index, gcnt);
        ps->ContraBaseVectorsP(el, mp.m_index, gcntp);
        for (int i = 0; i<neln; ++i)
            gradN[i] = (gcnt[0] * alpha + gcntp[0] * (1 - alpha))*(Gr[i]*m_s[iel]) +
            (gcnt[1] * alpha + gcntp[1] * (1 - alpha))*Gs[i];
        
        // calculate stiffness component
        int i = dof_a.index;
        int j = dof_b.index;
        vec3d v = gr*Gs[j] - gs*Gr[j];
        mat3d A; A.skew(v);
        mat3d Kv = vdotTdotv(gt, cv, gradN[j]);
        
        mat3d Kuu = (sv*A + Kv*M)*N[i] - A*(N[i] * p); Kuu *= -alpha;
        mat3d Kuw = Kv*N[i]; Kuw *= -alpha;
        vec3d kuJ = svJ*gt*(N[i] * N[j]) - f*(N[i] * N[j]); kuJ *= -alpha;
        
        Kab.zero();
        Kab.sub(0, 0, Kuu);
        Kab.sub(0, 3, Kuw);
        
        Kab[0][6] -= kuJ.x;
        Kab[1][6] -= kuJ.y;
        Kab[2][6] -= kuJ.z;
    });
    
    if (m_btied) {
        // In this supplemental section we evaluate the stiffness matrix acting on the solid
        // due to the traction imposed by the fluid-FSI domain.
        
        // Loop over all surface elements of the fluid-FSI surface
        FEElementMatrix ke;
        vector<int> en;
        for (int i=0; i<m_psurf->Elements(); ++i) {
            // get the surface element
            FESurfaceElement& elfsi = m_psurf->Element(i);
            int ielfsi = elfsi.m_lid;
            // get the fluid-FSI material
            FEFluidFSI* pfsi = dynamic_cast<FEFluidFSI*>(GetFEModel()->GetMaterial(m_elem[ielfsi]->GetMatID()));
            
            // find number of integration points on the fluid-FSI face
            int nifsi = elfsi.GaussPoints();
            
            // loop over integration points of the fluid-FSI face
            for (int k=0; k<nifsi; ++k) {
                // get the fluid-FSI surface material point
                FESurfaceMaterialPoint& mp = static_cast<FESurfaceMaterialPoint&>(*elfsi.GetMaterialPoint(k));
                
                // extract fluid-FSI data
                FETiedElasticSurface::Data& data = static_cast<FETiedElasticSurface::Data&>(*elfsi.GetMaterialPoint(k));
                // get surface element of solid surface
                FESurfaceElement* pe = data.m_pme;
                // if it exists, calculate the solid-solid stiffnes matrix acting at the ray intersection points
                if (pe != 0)
                {
                    FESurfaceElement& me = *pe;
                    vector<double> fe;
                    vector<int> lm;
                    int neln = me.Nodes();
                    int ndpn = 3;
                    int ndof = ndpn*neln;
                    ke.resize(ndof, ndof); ke.zero();
                    // build the en vector
                    en.resize(neln);
                    for (k=0; k<neln; ++k) en[k] = me.m_node[k];
                    lm.assign(ndof,0);
                    // evaluate lm
                    for (int i=0; i<neln; ++i) {
                        FENode& node = GetFEModel()->GetMesh().Node(me.m_node[i]);
                        vector<int>& id = node.m_ID;
                        // displacement dofs
                        for (int j=0; j<ndpn; ++j) lm[ndpn*i+j] = id[m_dofU[j]];
                    }
                    
                    vector<vec3d> gradN(neln);
                    
                    // evaluate fluid dilatation at the fluid-FSI integration point
                    double ef = GetFluidDilatation(mp, tp.alphaf);
                    // evaluate the corresponding fluid pressure
                    double p = pfsi->Fluid()->Pressure(ef);
                    
                    // get the nodal coordinates of the solid face
                    vec3d rt[FEElement::MAX_NODES];
                    m_psolid->GetNodalCoordinates(me, tp.alphaf, rt);
                    
                    // evaluate covariant basis vectors at the ray intersection with the solid face
                    vec3d gr = me.eval_deriv1(rt, data.m_rs[0], data.m_rs[1]);
                    vec3d gs = me.eval_deriv2(rt, data.m_rs[0], data.m_rs[1]);
                    vec3d gt = gr ^ gs;
                    
                    // Get the fluid stress and its tangents from integratino points of the fluid-FSI element
                    mat3ds sv(mat3dd(0));
                    tens4ds cv; cv.zero();
                    mat3d Ls; Ls.zero();
                    FEElement* pe = m_elem[ielfsi];
                    int pint = pe->GaussPoints();
                    for (int n = 0; n<pint; ++n)
                    {
                        FEMaterialPoint& mp = *pe->GetMaterialPoint(n);
                        FEElasticMaterialPoint& ep = *(mp.ExtractData<FEElasticMaterialPoint>());
                        sv += pfsi->Fluid()->GetViscous()->Stress(mp);
                        cv += pfsi->Fluid()->Tangent_RateOfDeformation(mp);
                        Ls += ep.m_L;
                    }
                    sv /= pint;
                    cv /= pint;
                    Ls /= pint;
                    mat3d M = mat3dd(a) - Ls;
                    
                    // get the the shape functions and its derivatives on the solid face (at ray intersection point)
                    double N[FEElement::MAX_INTPOINTS];
                    me.shape_fnc(N,data.m_rs[0], data.m_rs[1]);
                    double Gr[FEElement::MAX_INTPOINTS], Gs[FEElement::MAX_INTPOINTS];
                    me.shape_deriv(Gr,Gs,data.m_rs[0], data.m_rs[1]);
                    vec3d gcnt[2], gcntp[2];
                    m_psolid->ContraBaseVectors(me, data.m_rs[0], data.m_rs[1], gcnt);
                    m_psolid->ContraBaseVectorsP(me, data.m_rs[0], data.m_rs[1], gcntp);
                    for (int i = 0; i<neln; ++i)
                        gradN[i] = (gcnt[0] * alpha + gcntp[0] * (1 - alpha))*Gr[i] +
                        (gcnt[1] * alpha + gcntp[1] * (1 - alpha))*Gs[i];
                    
                    // evaluate the solid-solid stiffness matrix at ray intersection points
                    for (int m=0; m<neln; ++m) {
                        for (int n=0; n<neln; ++n) {
                            vec3d v = gr*Gs[n] - gs*Gr[n];
                            mat3d A; A.skew(v);
                            mat3d Kv = vdotTdotv(gt, cv, gradN[n]);
                            mat3d Kuu = (sv*A + Kv*M)*N[i] - A*(N[i] * p); Kuu *= -alpha;
                            ke[ndpn*m  ][ndpn*n  ] += Kuu(0,0); ke[ndpn*m  ][ndpn*n+1] += Kuu(0,1); ke[ndpn*m  ][ndpn*n+2] += Kuu(0,2);
                            ke[ndpn*m+1][ndpn*n  ] += Kuu(1,0); ke[ndpn*m  ][ndpn*n+1] += Kuu(1,1); ke[ndpn*m  ][ndpn*n+2] += Kuu(1,2);
                            ke[ndpn*m+2][ndpn*n  ] += Kuu(2,0); ke[ndpn*m  ][ndpn*n+1] += Kuu(2,1); ke[ndpn*m  ][ndpn*n+2] += Kuu(2,2);
                        }
                    }
                    
                    // assemble this element stiffness matrix into the global load stiffness matrix
                    ke.SetNodes(en);
                    ke.SetIndices(lm);
                    LS.Assemble(ke);
                }
            }
        }
    }
}

//-----------------------------------------------------------------------------
void FEFluidFSITraction::Serialize(DumpStream& ar)
{
    FESurfaceLoad::Serialize(ar);

	if (ar.IsShallow() == false)
	{
        ar & m_s;
		if (ar.IsSaving())
		{
			int NE = (int)m_elem.size();
			ar << NE;
			for (int i = 0; i < NE; ++i)
			{
				FEElement* pe = m_elem[i];
				int nid = (pe ? pe->GetID() : -1);
				ar << nid;
			}
		}
		else
		{
			FEMesh& mesh = ar.GetFEModel().GetMesh();
			int NE, nid;
			ar >> NE;
			m_elem.resize(NE, nullptr);
			for (int i = 0; i < NE; ++i)
			{
				ar >> nid;
				if (nid != -1)
				{
					FEElement* pe = mesh.FindElementFromID(nid);
					assert(pe);
					m_elem[i] = pe;
				}
			}
		}
	}
}
