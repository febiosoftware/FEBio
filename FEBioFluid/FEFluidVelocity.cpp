/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2020 University of Utah, The Trustees of Columbia University in 
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
#include "FEFluidVelocity.h"
#include <FECore/FEModel.h>
#include <FECore/log.h>
#include "FEBioFluid.h"

//=============================================================================
BEGIN_FECORE_CLASS(FEFluidVelocity, FESurfaceLoad)
	ADD_PARAMETER(m_scale, "scale");
	ADD_PARAMETER(m_VC   , "velocity"   );
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! constructor
FEFluidVelocity::FEFluidVelocity(FEModel* pfem) : FESurfaceLoad(pfem), m_VC(FE_VEC3D), m_dofW(pfem), m_dofEF(pfem)
{
    m_scale = 1.0;
}

//-----------------------------------------------------------------------------
//! allocate storage
void FEFluidVelocity::SetSurface(FESurface* ps)
{
    FESurfaceLoad::SetSurface(ps);
    m_VC.Create(ps->GetFacetSet());
}

//-----------------------------------------------------------------------------
vec3d FEFluidVelocity::FluidVelocity(FESurfaceMaterialPoint& mp)
{
	FESurfaceElement& el = *mp.SurfaceElement();
	double* N = mp.m_shape;
	vec3d v(0, 0, 0);
	int neln = el.Nodes();
	for (int i = 0; i<neln; ++i)
	{
		v += m_VN[el.m_lnode[i]] * N[i];
	}
	return v;
}

//-----------------------------------------------------------------------------
//! Calculate the residual for the prescribed normal component of velocity
void FEFluidVelocity::LoadVector(FEGlobalVector& R, const FETimeInfo& tp)
{
	m_psurf->LoadVector(R, m_dofEF, true, [=](FESurfaceMaterialPoint& mp, const FESurfaceDofShape& dof_a, vector<double>& fa) {

		// calculate the tangent vectors
		vec3d v = FluidVelocity(mp);

		vec3d nu = mp.dxr ^ mp.dxs;
		double da = nu.unit();
		double vn = (v*nu)*m_scale;

		double H = dof_a.shape;
		fa[0] = H * vn * da;
	});
}

//-----------------------------------------------------------------------------
//! initialize
bool FEFluidVelocity::Init()
{
	m_dofW.Clear();
	if (m_dofW.AddVariable(FEBioFluid::GetVariableName(FEBioFluid::RELATIVE_FLUID_VELOCITY)) == false) return false;
	m_dofEF.Clear();
	if (m_dofEF.AddVariable(FEBioFluid::GetVariableName(FEBioFluid::FLUID_DILATATION)) == false) return false;

    FESurfaceLoad::Init();
    
    // evaluate nodal velocities from boundary cards
	FESurface* ps = &GetSurface();
	m_VN.resize(ps->Nodes(),vec3d(0,0,0));
    vector<int> nf(ps->Nodes(),0);

	for (int iel=0; iel<ps->Elements(); ++iel)
    {
        FESurfaceElement& el = m_psurf->Element(iel);
        
        // nr of element nodes
        int neln = el.Nodes();
        
        for (int i=0; i<neln; ++i) 
		{
            m_VN[el.m_lnode[i]] += m_VC.value<vec3d>(iel, i);
            ++nf[el.m_lnode[i]];
        }
    }
    
    for (int i=0; i<ps->Nodes(); ++i) m_VN[i] /= nf[i];
    
    return true;
}

//-----------------------------------------------------------------------------
//! Activate the degrees of freedom for this BC
void FEFluidVelocity::Activate()
{
    FESurface* ps = &GetSurface();
    
    for (int i=0; i<ps->Nodes(); ++i)
    {
        FENode& node = ps->Node(i);
        // mark node as having prescribed DOF
        node.set_bc(m_dofW[0], DOF_PRESCRIBED);
        node.set_bc(m_dofW[1], DOF_PRESCRIBED);
        node.set_bc(m_dofW[2], DOF_PRESCRIBED);
    }
}

//-----------------------------------------------------------------------------
//! Evaluate and prescribe the nodal velocities
void FEFluidVelocity::Update()
{
    // prescribe this velocity at the nodes
    FESurface* ps = &GetSurface();
    
    for (int i=0; i<ps->Nodes(); ++i)
    {
        // evaluate the nodal velocity
        vec3d v = m_VN[i]*m_scale;
        FENode& node = ps->Node(i);
        if (node.m_ID[m_dofW[0]] < -1) node.set(m_dofW[0], v.x);
        if (node.m_ID[m_dofW[1]] < -1) node.set(m_dofW[1], v.y);
        if (node.m_ID[m_dofW[2]] < -1) node.set(m_dofW[2], v.z);
    }
}

//-----------------------------------------------------------------------------
//! serialization
void FEFluidVelocity::Serialize(DumpStream& ar)
{
	FESurfaceLoad::Serialize(ar);
	if (ar.IsShallow()) return;
	ar & m_bpv & m_dofW & m_dofEF;
	ar & m_VN;
}
