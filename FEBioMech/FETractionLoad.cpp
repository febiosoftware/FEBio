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
#include "FETractionLoad.h"
#include "FEBioMech.h"
#include <FECore/FEFacetSet.h>

//=============================================================================
BEGIN_FECORE_CLASS(FETractionLoad, FESurfaceLoad)
	ADD_PARAMETER(m_scale   , "scale");
	ADD_PARAMETER(m_traction, "traction");
	ADD_PARAMETER(m_bshellb , "shell_bottom");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! constructor
FETractionLoad::FETractionLoad(FEModel* pfem) : FESurfaceLoad(pfem), m_dofList(pfem)
{
	m_scale = 1.0;
	m_traction = vec3d(0, 0, 0);
	m_bshellb = false;
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
	m_dofList.Clear();
	if (m_bshellb == false)
	{
		m_dofList.AddVariable(FEBioMech::GetVariableName(FEBioMech::DISPLACEMENT));
	}
	else
	{
		m_dofList.AddVariable(FEBioMech::GetVariableName(FEBioMech::SHELL_DISPLACEMENT));
	}
	if (m_dofList.IsEmpty()) return false;

	return FESurfaceLoad::Init();
}

//-----------------------------------------------------------------------------
void FETractionLoad::LoadVector(FEGlobalVector& R, const FETimeInfo& tp)
{
	FESurface& surf = GetSurface();
	surf.SetShellBottom(m_bshellb);

	// evaluate the integral
	FETractionLoad* load = this;
	surf.LoadVector(R, m_dofList, true, [=](FESurfaceMaterialPoint& pt, int node_a, std::vector<double>& val) {

		// evaluate traction at this material point
		vec3d t = m_traction(pt)*m_scale;
		if (load->m_bshellb) t = -t;

		double J = (pt.dxr ^ pt.dxs).norm();

		FESurfaceElement& el = *pt.SurfaceElement();
		double* H = el.H(pt.m_index);

		val[0] = H[node_a] * t.x*J;
		val[1] = H[node_a] * t.y*J;
		val[2] = H[node_a] * t.z*J;
	});
}

//-----------------------------------------------------------------------------
void FETractionLoad::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp)
{
	// Nothing to do here.
}
