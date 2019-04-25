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
#include "FEPressureLoad.h"
#include "FEBioMech.h"
#include <FECore/FEFacetSet.h>

//-----------------------------------------------------------------------------
// Parameter block for pressure loads
BEGIN_FECORE_CLASS(FEPressureLoad, FESurfaceLoad)
	ADD_PARAMETER(m_pressure, "pressure");
	ADD_PARAMETER(m_bsymm   , "symmetric_stiffness");
	ADD_PARAMETER(m_blinear , "linear");
	ADD_PARAMETER(m_bshellb , "shell_bottom");
END_FECORE_CLASS()

//-----------------------------------------------------------------------------
//! constructor
FEPressureLoad::FEPressureLoad(FEModel* pfem) : FESurfaceLoad(pfem), m_dofList(pfem)
{ 
	m_pressure = 0.0;
	m_bsymm = true;
	m_bshellb = false;
}

//-----------------------------------------------------------------------------
//! allocate storage
void FEPressureLoad::SetSurface(FESurface* ps)
{
	FESurfaceLoad::SetSurface(ps);
	m_pressure.SetItemList(ps->GetFacetSet());
}

//-----------------------------------------------------------------------------
bool FEPressureLoad::Init()
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
void FEPressureLoad::Residual(FEGlobalVector& R, const FETimeInfo& tp)
{
	FESurface& surf = GetSurface();
	surf.SetShellBottom(m_bshellb);

	// evaluate the integral
	FEPressureLoad* load = this;
	surf.LoadVector(R, m_dofList, m_blinear, [=](FESurfaceMaterialPoint& pt, int node_a, std::vector<double>& val) {
		
		// evaluate pressure at this material point
		double P = -m_pressure(pt);
		if (load->m_bshellb) P = -P;

		double J = (pt.dxr ^ pt.dxs).norm();

		// force vector
		vec3d N = (pt.dxr ^ pt.dxs); N.unit();
		vec3d t = N*P;

		FESurfaceElement& el = *pt.SurfaceElement();

		double* H = el.H(pt.m_index);

		val[0] = H[node_a]*t.x*J;
		val[1] = H[node_a]*t.y*J;
		val[2] = H[node_a]*t.z*J;
	});
}

//-----------------------------------------------------------------------------
void FEPressureLoad::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp)
{
	// Don't calculate stiffness for a linear load
	if (m_blinear) return;

	FESurface& surf = GetSurface();
	surf.SetShellBottom(m_bshellb);

	// evaluate the integral
	FEPressureLoad* load = this;
	surf.LoadStiffness(LS, m_dofList, m_dofList, [=](FESurfaceMaterialPoint& mp, int node_a, int node_b, matrix& kab) {

		// evaluate pressure at this material point
		double P = -(load->m_pressure(mp));
		if (load->m_bshellb) P = -P;

		FESurfaceElement& el = *mp.SurfaceElement();

		double* N = el.H(mp.m_index);
		double* Gr = el.Gr(mp.m_index);
		double* Gs = el.Gs(mp.m_index);

		int i = node_a;
		int j = node_b;

		vec3d vab(0,0,0);
		if (m_bsymm)
			vab = (mp.dxr*(N[j] * Gs[i] - N[i] * Gs[j]) - mp.dxs*(N[j] * Gr[i] - N[i] * Gr[j])) * 0.5*P;
		else  
			vab = (mp.dxs*Gr[j] - mp.dxr*Gs[j])*(P*N[i]);

		mat3da K(vab);
		kab.set(0, 0, K);
	});
}
