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
#include "FETangentialDamping.h"
#include <FECore/FEMesh.h>
#include "FEBioFluid.h"
#include <FECore/FELinearSystem.h>

//-----------------------------------------------------------------------------
// Parameter block for pressure loads
BEGIN_FECORE_CLASS(FETangentialDamping, FESurfaceLoad)
	ADD_PARAMETER(m_eps, "penalty");
END_FECORE_CLASS()

//-----------------------------------------------------------------------------
//! constructor
FETangentialDamping::FETangentialDamping(FEModel* pfem) : FESurfaceLoad(pfem), m_dofW(pfem)
{
    m_eps = 0.0;
    
    // get the degrees of freedom
	m_dofW.AddVariable(FEBioFluid::GetVariableName(FEBioFluid::RELATIVE_FLUID_VELOCITY));
}

//-----------------------------------------------------------------------------
//! allocate storage
void FETangentialDamping::SetSurface(FESurface* ps)
{
    FESurfaceLoad::SetSurface(ps);
}

//-----------------------------------------------------------------------------
vec3d FETangentialDamping::FluidVelocity(FESurfaceMaterialPoint& mp, double alpha)
{
	vec3d vt[FEElement::MAX_NODES];
	FESurfaceElement& el = *mp.SurfaceElement();
	int neln = el.Nodes();
	for (int j = 0; j<neln; ++j) {
		FENode& node = m_psurf->Node(el.m_lnode[j]);
		vt[j] = node.get_vec3d(m_dofW[0], m_dofW[1], m_dofW[2])*alpha + node.get_vec3d_prev(m_dofW[0], m_dofW[1], m_dofW[2])*(1. - alpha);
	}
	return el.eval(vt, mp.m_index);
}

//-----------------------------------------------------------------------------
void FETangentialDamping::Residual(FEGlobalVector& R, const FETimeInfo& tp)
{
	m_psurf->LoadVector(R, m_dofW, false, [=](FESurfaceMaterialPoint& mp, int node_a, vector<double>& fa) {

		// fluid velocity
		vec3d v = FluidVelocity(mp, tp.alpha);

		vec3d n = mp.dxr ^ mp.dxs;
		double da = n.unit();

		// force vector
		mat3dd I(1);
		vec3d f = (I - dyad(n))*v*(-m_eps*da);

		double* N = mp.m_shape;
		fa[0] = N[node_a] * f.x;
		fa[1] = N[node_a] * f.y;
		fa[2] = N[node_a] * f.z;
	});
}

//-----------------------------------------------------------------------------
void FETangentialDamping::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp)
{
	m_psurf->LoadStiffness(LS, m_dofW, m_dofW, [=](FESurfaceMaterialPoint& mp, int node_a, int node_b, matrix& Kab) {
   
        vec3d n = mp.dxr ^ mp.dxs;
        double da = n.unit();
        
		mat3dd I(1);
		mat3ds K = (I - dyad(n))*(-m_eps*da);

		// shape functions
		double* N = mp.m_shape;
        
        // calculate stiffness component
		mat3ds kab = K*(N[node_a]*N[node_b]);
		Kab.set(0, 0, kab);
	});
}
