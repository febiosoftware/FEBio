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
#include "FEFluidTractionLoad.h"
#include <FECore/FESurface.h>
#include <FECore/FEFacetSet.h>
#include <FECore/FEMesh.h>
#include "FEBioFluid.h"

//=============================================================================
BEGIN_FECORE_CLASS(FEFluidTractionLoad, FESurfaceLoad)
	ADD_PARAMETER(m_scale, "scale"   );
	ADD_PARAMETER(m_TC   , "traction");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! constructor
FEFluidTractionLoad::FEFluidTractionLoad(FEModel* pfem) : FESurfaceLoad(pfem), m_TC(FE_VEC3D), m_dofW(pfem)
{
	m_scale = 1.0;
}

//-----------------------------------------------------------------------------
//! allocate storage
void FEFluidTractionLoad::SetSurface(FESurface* ps)
{
	FESurfaceLoad::SetSurface(ps);
	m_TC.Create(ps); 
}

//-----------------------------------------------------------------------------
bool FEFluidTractionLoad::Init()
{
	m_dofW.Clear();
	if (m_dofW.AddVariable(FEBioFluid::GetVariableName(FEBioFluid::RELATIVE_FLUID_VELOCITY)) == false) return false;
	return FESurfaceLoad::Init();
}

//-----------------------------------------------------------------------------
vec3d FEFluidTractionLoad::TractionLoad(FESurfaceMaterialPoint& mp)
{
	FESurfaceElement& el = *mp.SurfaceElement();
	int iel = el.m_lid;
	int neln = el.Nodes();
	vec3d tn[FEElement::MAX_NODES];
	for (int i = 0; i<neln; ++i)
	{
		tn[i] = m_TC.value<vec3d>(iel, i)*m_scale;
	}
	return el.eval(tn, mp.m_index);
}

//-----------------------------------------------------------------------------
//! Calculate the residual for the traction load
void FEFluidTractionLoad::Residual(FEGlobalVector& R, const FETimeInfo& tp)
{
	m_psurf->LoadVector(R, m_dofW, true, [=](FESurfaceMaterialPoint& mp, int node_a, vector<double>& fa) {

		// fluid traction
		vec3d t = TractionLoad(mp);
		vec3d f = t*((mp.dxr ^ mp.dxs).norm());

		double* N = mp.m_shape;
		fa[0] = N[node_a] * f.x;
		fa[1] = N[node_a] * f.y;
		fa[2] = N[node_a] * f.z;
	});
}

//! calculate traction stiffness (there is none)
void FEFluidTractionLoad::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp)
{

}
