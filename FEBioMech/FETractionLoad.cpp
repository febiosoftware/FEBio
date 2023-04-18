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
#include "FETractionLoad.h"
#include "FEBioMech.h"
#include <FECore/FEFacetSet.h>

//=============================================================================
BEGIN_FECORE_CLASS(FETractionLoad, FESurfaceLoad)
	ADD_PARAMETER(m_scale   , "scale")->SetFlags(FE_PARAM_ADDLC | FE_PARAM_VOLATILE);
	ADD_PARAMETER(m_traction, "traction")->setUnits(UNIT_PRESSURE);
	ADD_PARAMETER(m_bshellb , "shell_bottom");
	ADD_PARAMETER(m_blinear, "linear");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! constructor
FETractionLoad::FETractionLoad(FEModel* pfem) : FESurfaceLoad(pfem)
{
	m_scale = 1.0;
	m_traction = vec3d(0, 0, 0);
	m_bshellb = false;
	m_blinear = false;
}

//-----------------------------------------------------------------------------
//! allocate storage
void FETractionLoad::SetSurface(FESurface* ps)
{
	FESurfaceLoad::SetSurface(ps);
	m_traction.SetItemList(ps->GetFacetSet());
}

//-----------------------------------------------------------------------------
// initialization
bool FETractionLoad::Init()
{
	// get the degrees of freedom
	m_dof.Clear();
	if (m_bshellb == false)
	{
		m_dof.AddVariable(FEBioMech::GetVariableName(FEBioMech::DISPLACEMENT));
	}
	else
	{
		m_dof.AddVariable(FEBioMech::GetVariableName(FEBioMech::SHELL_DISPLACEMENT));
	}
	if (m_dof.IsEmpty()) return false;

	return FESurfaceLoad::Init();
}

//-----------------------------------------------------------------------------
void FETractionLoad::LoadVector(FEGlobalVector& R)
{
	FESurface& surf = GetSurface();
	surf.SetShellBottom(m_bshellb);

	// evaluate the integral
	FETractionLoad* load = this;
	surf.LoadVector(R, m_dof, m_blinear, [=](FESurfaceMaterialPoint& pt, const FESurfaceDofShape& dof_a, std::vector<double>& val) {

		// evaluate traction at this material point
		vec3d t = m_traction(pt)*m_scale;
		if (load->m_bshellb) t = -t;

		double J = (pt.dxr ^ pt.dxs).norm();

		double H_u = dof_a.shape;

		val[0] = H_u * t.x*J;
		val[1] = H_u * t.y*J;
		val[2] = H_u * t.z*J;
	});
}

//-----------------------------------------------------------------------------
void FETractionLoad::StiffnessMatrix(FELinearSystem& LS)
{
	// Nothing to do here.
	// TODO: I think if the linear flag is false, I do need to evaluate a stiffness.
}
