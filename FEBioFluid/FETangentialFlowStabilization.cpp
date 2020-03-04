/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, The Trustees of Columbia University in 
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
#include "FETangentialFlowStabilization.h"
#include "FEFluidMaterial.h"
#include "FEBioFluid.h"
#include <FECore/FELinearSystem.h>

//-----------------------------------------------------------------------------
// Parameter block for pressure loads
BEGIN_FECORE_CLASS(FETangentialFlowStabilization, FESurfaceLoad)
    ADD_PARAMETER(m_beta, "beta");
END_FECORE_CLASS()

//-----------------------------------------------------------------------------
//! constructor
FETangentialFlowStabilization::FETangentialFlowStabilization(FEModel* pfem) : FESurfaceLoad(pfem), m_dofU(pfem), m_dofW(pfem)
{
    m_beta = 1.0;
    m_rho = 1.0;
}

//-----------------------------------------------------------------------------
//! allocate storage
void FETangentialFlowStabilization::SetSurface(FESurface* ps)
{
    FESurfaceLoad::SetSurface(ps);
}

//-----------------------------------------------------------------------------
//! initialize
bool FETangentialFlowStabilization::Init()
{
	// get the degrees of freedom
	m_dofU.Clear();
	if (m_dofU.AddVariable(FEBioFluid::GetVariableName(FEBioFluid::DISPLACEMENT)) == false) return false;

	m_dofW.Clear();
	if (m_dofW.AddVariable(FEBioFluid::GetVariableName(FEBioFluid::RELATIVE_FLUID_VELOCITY)) == false) return false;

	m_dof = m_dofW;

    if (FESurfaceLoad::Init() == false) return false;
    
    // get fluid density from first surface element
    // assuming the entire surface bounds the same fluid
    FESurfaceElement& el = m_psurf->Element(0);
    FEMesh* mesh = m_psurf->GetMesh();
    FEElement* pe = el.m_elem[0];
    if (pe == nullptr) return false;

	// get the material
    FEMaterial* pm = GetFEModel()->GetMaterial(pe->GetMatID());
    FEFluidMaterial* fluid = pm->ExtractProperty<FEFluidMaterial>();
	if (fluid == nullptr) return false;

	// get the density and bulk modulus
    m_rho = fluid->ReferentialDensity();
    
    return true;
}

//-----------------------------------------------------------------------------
void FETangentialFlowStabilization::Serialize(DumpStream& ar)
{
    FESurfaceLoad::Serialize(ar);
	ar & m_rho;
}

//-----------------------------------------------------------------------------
vec3d FETangentialFlowStabilization::FluidVelocity(FESurfaceMaterialPoint& mp, double alpha)
{
	FESurfaceElement& el = *mp.SurfaceElement();
	vec3d vt[FEElement::MAX_NODES];
	int neln = el.Nodes();
	for (int j = 0; j<neln; ++j) {
		FENode& node = m_psurf->Node(el.m_lnode[j]);
		vt[j] = node.get_vec3d(m_dofW[0], m_dofW[1], m_dofW[2])*alpha + node.get_vec3d_prev(m_dofW[0], m_dofW[1], m_dofW[2])*(1 - alpha);
	}
	return el.eval(vt, mp.m_index);
}

//-----------------------------------------------------------------------------
void FETangentialFlowStabilization::LoadVector(FEGlobalVector& R, const FETimeInfo& tp)
{
	m_psurf->LoadVector(R, m_dofW, false, [=](FESurfaceMaterialPoint& mp, const FESurfaceDofShape& dof_a, vector<double>& fa) {

		vec3d n = mp.dxr ^ mp.dxs;
		double da = n.unit();

		// fluid velocity
		vec3d v = FluidVelocity(mp, tp.alpha);

		// tangential traction = -beta*density*|tangential velocity|*(tangential velocity)
		mat3dd I(1.0);
		vec3d vtau = (I - dyad(n))*v;
		double vmag = vtau.norm();

		// force vector (change sign for inflow vs outflow)
		vec3d f = vtau*(-m_beta*m_rho*vmag*da);

		double H = dof_a.shape;
		fa[0] = H * f.x;
		fa[1] = H * f.y;
		fa[2] = H * f.z;
	});
}

//-----------------------------------------------------------------------------
void FETangentialFlowStabilization::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp)
{
	FEDofList dofs(GetFEModel());
	dofs.AddDofs(m_dofU);
	dofs.AddDofs(m_dofW);
	m_psurf->LoadStiffness(LS, dofs, dofs, [=](FESurfaceMaterialPoint& mp, const FESurfaceDofShape& dof_a, const FESurfaceDofShape& dof_b, matrix& Kab) {
    
		FESurfaceElement& el = *mp.SurfaceElement();
		double alpha = tp.alpha;

		// fluid velocity
		vec3d v = FluidVelocity(mp, alpha);

		// tangent vectors
		vec3d rt[FEElement::MAX_NODES];
		m_psurf->GetNodalCoordinates(el, alpha, rt);
		vec3d dxr = el.eval_deriv1(rt, mp.m_index);
		vec3d dxs = el.eval_deriv2(rt, mp.m_index);
        
        vec3d n = dxr ^ dxs;
        double da = n.unit();
        
		mat3dd I(1.0);
		vec3d vtau = (I - dyad(n))*v;
        double vmag = vtau.unit();
        mat3d K = (I - dyad(n) + dyad(vtau))*(-m_beta*m_rho*vmag*da);
        // force vector (change sign for inflow vs outflow)
        vec3d ttt = vtau*(-m_beta*m_rho*vmag);

		// shape functions and derivatives
		double H_i  = dof_a.shape;
		double Gr_i = dof_a.shape_deriv_r;
		double Gs_i = dof_a.shape_deriv_s;

		double H_j  = dof_b.shape;
		double Gr_j = dof_b.shape_deriv_r;
		double Gs_j = dof_b.shape_deriv_s;

        // calculate stiffness component
		mat3d Kww = K*(H_i*H_j*alpha);
		vec3d g = (dxr*Gs_j - dxs*Gr_j)*(H_i*alpha);
		mat3d Kwu; Kwu.skew(g);
		Kwu = (ttt & n)*Kwu;
		Kab.zero();

		// dw/du
		Kab.sub(3, 0, Kwu);

		// dw/dw
		Kab.sub(3, 3, Kww);
	});
}
