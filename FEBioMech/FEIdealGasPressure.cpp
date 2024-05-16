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
#include "FEIdealGasPressure.h"
#include "FEBioMech.h"
#include <FECore/log.h>

BEGIN_FECORE_CLASS(FEIdealGasPressure, FESurfaceLoad)
	ADD_PARAMETER(m_initialPressure, "P0")->setUnits(UNIT_PRESSURE);
	ADD_PARAMETER(m_bsymm, "symmetric_stiffness");
	ADD_PARAMETER(m_bshellb, "shell_bottom");
END_FECORE_CLASS()

FEIdealGasPressure::FEIdealGasPressure(FEModel* pfem) : FESurfaceLoad(pfem)
{
	m_initialPressure = 0.0;
	m_bsymm = true;
	m_bshellb = false;
	m_currentPressure = 0;
	m_currentVolume = 0;
}

bool FEIdealGasPressure::Init()
{
	FESurface& surf = GetSurface();
	surf.SetShellBottom(m_bshellb);

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

void FEIdealGasPressure::Activate()
{
	m_initialVolume = CalculateSurfaceVolume(GetSurface());
	m_currentVolume = m_initialVolume;
	
	m_currentPressure = m_initialPressure;

	feLogDebug("initial volume: %lg\n", m_initialVolume);
	
	FESurfaceLoad::Activate();
}

void FEIdealGasPressure::Update()
{
	m_currentVolume = CalculateSurfaceVolume(GetSurface());
	m_currentPressure = m_initialPressure * (m_initialVolume / m_currentVolume);

	feLogDebug("volume: %lg, pressure: %lg\n", m_currentVolume, m_currentPressure);

	FESurfaceLoad::Update();
}

void FEIdealGasPressure::LoadVector(FEGlobalVector& R)
{
	// evaluate the integral
	FESurface& surf = GetSurface();
	surf.LoadVector(R, m_dof, false, [&](FESurfaceMaterialPoint& pt, const FESurfaceDofShape& dof_a, std::vector<double>& val) {

		// evaluate pressure at this material point
		double P = m_currentPressure;
		if (m_bshellb) P = -P;

		double J = (pt.dxr ^ pt.dxs).norm();

		// force vector
		vec3d N = (pt.dxr ^ pt.dxs); N.unit();
		vec3d t = N * P;

		double H_u = dof_a.shape;

		val[0] = H_u * t.x * J;
		val[1] = H_u * t.y * J;
		val[2] = H_u * t.z * J;
		});
}

void FEIdealGasPressure::StiffnessMatrix(FELinearSystem& LS)
{
	// evaluate the integral
	FESurface& surf = GetSurface();
	surf.LoadStiffness(LS, m_dof, m_dof, [&](FESurfaceMaterialPoint& mp, const FESurfaceDofShape& dof_a, const FESurfaceDofShape& dof_b, matrix& kab) {

		// evaluate pressure at this material point
		double P = m_currentPressure;
		if (m_bshellb) P = -P;

		double H_i = dof_a.shape;
		double Gr_i = dof_a.shape_deriv_r;
		double Gs_i = dof_a.shape_deriv_s;

		double H_j = dof_b.shape;
		double Gr_j = dof_b.shape_deriv_r;
		double Gs_j = dof_b.shape_deriv_s;

		vec3d vab(0, 0, 0);
		if (m_bsymm)
			vab = (mp.dxr * (H_j * Gs_i - H_i * Gs_j) - mp.dxs * (H_j * Gr_i - H_i * Gr_j)) * 0.5 * P;
		else
			vab = (mp.dxs * Gr_j - mp.dxr * Gs_j) * (P * H_i);

		mat3da K(vab);
		kab.set(0, 0, K);
		});
}
