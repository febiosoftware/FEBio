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
FEFluidTractionLoad::FEFluidTractionLoad(FEModel* pfem) : FESurfaceLoad(pfem)
{
	m_scale = 1.0;
	m_TC = vec3d(0,0,0);
}

//-----------------------------------------------------------------------------
//! allocate storage
void FEFluidTractionLoad::SetSurface(FESurface* ps)
{
	FESurfaceLoad::SetSurface(ps);
	m_TC.SetItemList(ps->GetFacetSet()); 
}

//-----------------------------------------------------------------------------
bool FEFluidTractionLoad::Init()
{
	m_dof.Clear();
	if (m_dof.AddVariable(FEBioFluid::GetVariableName(FEBioFluid::RELATIVE_FLUID_VELOCITY)) == false) return false;
	return FESurfaceLoad::Init();
}

//-----------------------------------------------------------------------------
//! Calculate the residual for the traction load
void FEFluidTractionLoad::LoadVector(FEGlobalVector& R)
{
	m_psurf->LoadVector(R, m_dof, true, [&](FESurfaceMaterialPoint& mp, const FESurfaceDofShape& dof_a, vector<double>& fa) {

		// fluid traction
		vec3d t = m_TC(mp)*m_scale;
		vec3d f = t*((mp.dxr ^ mp.dxs).norm());

		double H = dof_a.shape;
		fa[0] = H * f.x;
		fa[1] = H * f.y;
		fa[2] = H * f.z;
	});
}

//-----------------------------------------------------------------------------
//! calculate traction stiffness (there is none)
void FEFluidTractionLoad::StiffnessMatrix(FELinearSystem& LS)
{

}
