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
#include "FEFluidNormalTraction.h"
#include "FEBioFluid.h"
#include <FECore/FESurface.h>
#include <FECore/FEFacetSet.h>

//=============================================================================
BEGIN_FECORE_CLASS(FEFluidNormalTraction, FESurfaceLoad)
	ADD_PARAMETER(m_traction, "traction");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! constructor
FEFluidNormalTraction::FEFluidNormalTraction(FEModel* pfem) : FESurfaceLoad(pfem), m_dofW(pfem)
{
	m_traction = 1.0;
}

//-----------------------------------------------------------------------------
//! allocate storage
void FEFluidNormalTraction::SetSurface(FESurface* ps)
{
    FESurfaceLoad::SetSurface(ps);
	m_traction.SetItemList(ps->GetFacetSet());
}

//-----------------------------------------------------------------------------
//! initialization
bool FEFluidNormalTraction::Init()
{
	m_dofW.Clear();
	if (m_dofW.AddVariable(FEBioFluid::GetVariableName(FEBioFluid::RELATIVE_FLUID_VELOCITY))) return false;
	return FESurfaceLoad::Init();
}

//-----------------------------------------------------------------------------
void FEFluidNormalTraction::Serialize(DumpStream& ar)
{
    FESurfaceLoad::Serialize(ar);
    if (ar.IsShallow()) return;
    ar & m_dofW;
}

//-----------------------------------------------------------------------------
//! Calculate the residual for the traction load
void FEFluidNormalTraction::LoadVector(FEGlobalVector& R)
{
	const FETimeInfo& tp = GetTimeInfo();

	// evaluate integral over surface
	m_psurf->LoadVector(R, m_dofW, false, [&](FESurfaceMaterialPoint& mp, const FESurfaceDofShape &dof_a, vector<double>& fa) {

		FESurfaceElement& el = *mp.SurfaceElement();

		vec3d rt[FEElement::MAX_NODES];
		m_psurf->GetNodalCoordinates(el, tp.alphaf, rt);

		// calculate the tangent vectors
		vec3d dxr = el.eval_deriv1(rt, mp.m_index);
		vec3d dxs = el.eval_deriv2(rt, mp.m_index);
		vec3d normal = dxr ^ dxs;

		double tn = m_traction(mp);
		vec3d f = normal*tn;

		double H = dof_a.shape;
		fa[0] = H * f.x;
		fa[1] = H * f.y;
		fa[2] = H * f.z;
	});
}
