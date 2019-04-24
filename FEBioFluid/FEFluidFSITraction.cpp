/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, Columbia University, and others.

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
#include "FECore/FEModel.h"
#include "FEFluid.h"
#include "FEFluidFSI.h"
#include "FEBioFSI.h"

//-----------------------------------------------------------------------------
//! constructor
FEFluidFSITraction::FEFluidFSITraction(FEModel* pfem) : FESurfaceLoad(pfem), m_dofU(pfem), m_dofSU(pfem), m_dofW(pfem)
{
    // get the degrees of freedom
	m_dofU.AddVariable(FEBioFSI::GetVariableName(FEBioFSI::DISPLACEMENT));
	m_dofSU.AddVariable(FEBioFSI::GetVariableName(FEBioFSI::SHELL_DISPLACEMENT));
	m_dofW.AddVariable(FEBioFSI::GetVariableName(FEBioFSI::RELATIVE_FLUID_VELOCITY));
    m_dofEF = pfem->GetDOFIndex(FEBioFSI::GetVariableName(FEBioFSI::FLUID_DILATATION), 0);
}

//-----------------------------------------------------------------------------
//! initialize
bool FEFluidFSITraction::Init()
{
	FESurface& surf = GetSurface();
	surf.SetInterfaceStatus(true);
	if (FESurfaceLoad::Init() == false) return false;

	// get the list of fluid-FSI elements connected to this interface
	FEModel* fem = GetFEModel();
	FEMesh* mesh = surf.GetMesh();
	int NF = surf.Elements();
	m_elem.resize(NF);
	m_K.resize(NF, 0);
	m_s.resize(NF, 1);
	m_bself.resize(NF, false);
	for (int j = 0; j<NF; ++j)
	{
		FESurfaceElement& el = surf.Element(j);
		// extract the first of two elements on this interface
		m_elem[j] = el.m_elem[0];
		if (el.m_elem[1] == nullptr) m_bself[j] = true;
		// get its material and check if FluidFSI
		FEMaterial* pm = fem->GetMaterial(m_elem[j]->GetMatID());
		FEFluidFSI* pfsi = dynamic_cast<FEFluidFSI*>(pm);
		if (pfsi) {
			m_K[j] = pfsi->Fluid()->m_k;
		}
		else if (!m_bself[j]) {
			// extract the second of two elements on this interface
			m_elem[j] = el.m_elem[1];
			pm = fem->GetMaterial(m_elem[j]->GetMatID());
			pfsi = dynamic_cast<FEFluidFSI*>(pm);
			if (pfsi == nullptr) return false;
			m_s[j] = -1;
			m_K[j] = pfsi->Fluid()->m_k;
		}
		else
			return false;
	}

    // TODO: Deal with the case when the surface is a shell domain separating two FSI domains
    // that use different fluid bulk moduli
    
    return true;
}

//-----------------------------------------------------------------------------
void FEFluidFSITraction::UnpackLM(FEElement& el, vector<int>& lm)
{
    FEMesh& mesh = *GetSurface().GetMesh();
    FESurfaceElement& fel = dynamic_cast<FESurfaceElement&>(el);
    FEElement* pe = fel.m_elem[0];
    // get the material
    FEMaterial* pm = GetFEModel()->GetMaterial(pe->GetMatID());
    FEFluidFSI* fsi = dynamic_cast<FEFluidFSI*> (pm);
    if (fsi == nullptr) pe = fel.m_elem[1];
    int N = el.Nodes();
    lm.resize(N*7);
    for (int i=0; i<N; ++i)
    {
        int n = el.m_node[i];
        FENode& node = mesh.Node(n);
        vector<int>& id = node.m_ID;
        
        lm[7*i  ] = id[m_dofU[0]];
        lm[7*i+1] = id[m_dofU[1]];
        lm[7*i+2] = id[m_dofU[2]];
        lm[7*i+3] = id[m_dofW[0]];
        lm[7*i+4] = id[m_dofW[1]];
        lm[7*i+5] = id[m_dofW[2]];
        lm[7*i+6] = id[m_dofEF];
    }

    // substitute interface dofs for solid-shell interfaces
    FESolidElement& sel = static_cast<FESolidElement&>(*pe);
    for (int i = 0; i<sel.m_bitfc.size(); ++i)
    {
        if (sel.m_bitfc[i]) {
            FENode& node = mesh.Node(sel.m_node[i]);
            vector<int>& id = node.m_ID;
            int j = el.FindNode(node.GetID()-1);
            
            // first the displacement dofs
            lm[7*j  ] = id[m_dofSU[0]];
            lm[7*j+1] = id[m_dofSU[1]];
            lm[7*j+2] = id[m_dofSU[2]];
        }
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
void FEFluidFSITraction::Residual(FEGlobalVector& R, const FETimeInfo& tp)
{
	FEModel* fem = GetFEModel();
	FEMesh* mesh = m_psurf->GetMesh();

	// TODO: If surface is bottom of shell, we should take shell displacement dofs (i.e. m_dofSU).
	m_psurf->LoadVector(R, m_dofU, false, [&](FESurfaceMaterialPoint& mp, int node_a, vector<double>& fa) {

		// get the surface element
		FESurfaceElement& el = *mp.SurfaceElement();
		int iel = el.m_lid;

		// nodal coordinates
		vec3d rt[FEElement::MAX_NODES];
		m_psurf->GetNodalCoordinates(el, tp.alpha, rt);

		// evaluate covariant basis vectors at integration point
		vec3d gr = el.eval_deriv1(rt, mp.m_index);
		vec3d gs = el.eval_deriv2(rt, mp.m_index);
		vec3d gt = gr ^ gs;

		// Get the fluid stress at integration point
		mat3ds sv = GetFluidStress(mp);

		// fluid dilatation at integration point
		double ef = GetFluidDilatation(mp, tp.alphaf);

		// evaluate traction
		vec3d f = (sv*gt + gt*(m_K[iel] * ef))*(-m_s[iel]);

		double* N = el.H(mp.m_index);
		fa[0] = N[node_a] * f.x;
		fa[1] = N[node_a] * f.y;
		fa[2] = N[node_a] * f.z;
	});
}

//-----------------------------------------------------------------------------
void FEFluidFSITraction::StiffnessMatrix(FESolver* psolver, const FETimeInfo& tp)
{
	FEModel* fem = GetFEModel();
	FESurface* ps = &GetSurface();

	// build dof list
	// TODO: If surface is bottom of shell, we should take shell displacement dofs (i.e. m_dofSU).
	FEDofList dofs(fem);
	dofs.AddDofs(m_dofU);
	dofs.AddDofs(m_dofW);
	dofs.AddDof(m_dofEF);

	// evaluate stiffness
	m_psurf->LoadStiffness(psolver, dofs, dofs, [&](FESurfaceMaterialPoint& mp, int node_a, int node_b, matrix& Kab) {

		FESurfaceElement& el = *mp.SurfaceElement();
		int iel = el.m_lid;
		int neln = el.Nodes();

		double dt = tp.timeIncrement;
		double alpha = tp.alpha;
		double a = tp.gamma / (tp.beta*dt);

		vector<vec3d> gradN(neln);

		// nodal coordinates
		vec3d rt[FEElement::MAX_NODES];
		ps->GetNodalCoordinates(el, tp.alpha, rt);

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
		vec3d gr = el.eval_deriv1(rt, mp.m_index);
		vec3d gs = el.eval_deriv2(rt, mp.m_index);
		vec3d gt = gr ^ gs;

		// evaluate fluid pressure
		double p = m_K[iel] * ef * m_s[iel];

		vec3d f = gt*(-m_K[iel] * m_s[iel]);

		vec3d gcnt[2], gcntp[2];
		ps->ContraBaseVectors(el, mp.m_index, gcnt);
		ps->ContraBaseVectorsP(el, mp.m_index, gcntp);
		for (int i = 0; i<neln; ++i)
			gradN[i] = (gcnt[0] * alpha + gcntp[0] * (1 - alpha))*Gr[i] +
			(gcnt[1] * alpha + gcntp[1] * (1 - alpha))*Gs[i];

		// calculate stiffness component
		int i = node_a;
		int j = node_b;
		vec3d v = gr*Gs[j] - gs*Gr[j];
		mat3d A; A.skew(v);
		mat3d Kv = vdotTdotv(gt, cv, gradN[j]);

		mat3d Kuu = (sv*A + Kv*M)*(-N[i] *  m_s[iel]) - A*(N[i] * p); Kuu *= alpha;
		mat3d Kuw = Kv*(-N[i] * m_s[iel]); Kuw *= alpha;
		vec3d kuJ = svJ*gt*(-N[i] * N[j] * m_s[iel]) + f*(N[i] * N[j]); kuJ *= alpha;

		Kab.zero();
		Kab[0][0] -= Kuu(0, 0); Kab[0][1] -= Kuu(0, 1); Kab[0][2] -= Kuu(0, 2);
		Kab[1][0] -= Kuu(1, 0); Kab[1][1] -= Kuu(1, 1); Kab[1][2] -= Kuu(1, 2);
		Kab[2][0] -= Kuu(2, 0); Kab[2][1] -= Kuu(2, 1); Kab[2][2] -= Kuu(2, 2);

		Kab[0][3] -= Kuw(0, 0); Kab[0][4] -= Kuw(0, 1); Kab[0][5] -= Kuw(0, 2);
		Kab[1][3] -= Kuw(1, 0); Kab[1][4] -= Kuw(1, 1); Kab[1][5] -= Kuw(1, 2);
		Kab[2][3] -= Kuw(2, 0); Kab[2][4] -= Kuw(2, 1); Kab[2][5] -= Kuw(2, 2);

		Kab[0][6] -= kuJ.x;
		Kab[1][6] -= kuJ.y;
		Kab[2][6] -= kuJ.z;
	});
}

//-----------------------------------------------------------------------------
void FEFluidFSITraction::Serialize(DumpStream& ar)
{
    FESurfaceLoad::Serialize(ar);

	ar & m_K;
	ar & m_s;
	ar & m_bself;

	if (ar.IsShallow() == false)
	{
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
