/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2026 University of Utah, The Trustees of Columbia University in
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
#include "FEBioMech.h"
#include "FEScriptedPressureLoad.h"

BEGIN_FECORE_CLASS(FEScriptedPressureLoad, FESurfaceLoad)
	ADD_PARAMETER(m_scriptName, "script");
END_FECORE_CLASS()

FEScriptedPressureLoad::FEScriptedPressureLoad(FEModel* pfem) : FESurfaceLoad(pfem), FEPhysicsProperty(pfem)
{
}

bool FEScriptedPressureLoad::Init()
{
	m_dof.AddVariable(FEBioMech::GetVariableName(FEBioMech::DISPLACEMENT));

	SetSibling(this);
	AddVariable("pos" , FEValueType::Vec3d);
	AddVariable("time", FEValueType::Double, false);

	if (FESurfaceLoad::Init() == false) return false;
	if (FEPhysicsProperty::Init() == false) return false;
	return true;
}

void FEScriptedPressureLoad::LoadVector(FEGlobalVector& R)
{
	double t = GetTimeInfo().currentTime;

	// evaluate the integral
	FESurface& surf = GetSurface();
	surf.LoadVector(R, m_dof, false, [&](FESurfaceMaterialPoint& pt, const FESurfaceDofShape& dof_a, std::vector<double>& val) {

		std::vector<FEValue> vars(2);
		vars[0] = pt.m_rt;
		vars[1] = t;

		// evaluate pressure at this material point
		double P = -Value(vars).d;

		// force vector
		vec3d N = (pt.dxr ^ pt.dxs);
		vec3d t = N * P;

		double H_u = dof_a.shape;

		val[0] = H_u * t.x;
		val[1] = H_u * t.y;
		val[2] = H_u * t.z;
	});
}

void FEScriptedPressureLoad::StiffnessMatrix(FELinearSystem& LS)
{
	double t = GetTimeInfo().currentTime;

	// evaluate the integral
	FESurface& surf = GetSurface();
	surf.LoadStiffness(LS, m_dof, m_dof, [&](FESurfaceMaterialPoint& mp, const FESurfaceDofShape& dof_a, const FESurfaceDofShape& dof_b, matrix& kab) {

		std::vector<FEValue> vars(2);
		vars[0] = mp.m_rt;
		vars[1] = t;

		// evaluate pressure at this material point
		double P = -Value(vars).d;

		// evaluate pressure gradient at this material point
		vec3d dP = -DerivValue(vars, 0).v3;

		double H_i = dof_a.shape;
		double Gr_i = dof_a.shape_deriv_r;
		double Gs_i = dof_a.shape_deriv_s;

		double H_j = dof_b.shape;
		double Gr_j = dof_b.shape_deriv_r;
		double Gs_j = dof_b.shape_deriv_s;

		vec3d vab(0, 0, 0);
		vab = (mp.dxs * Gr_j - mp.dxr * Gs_j) * (P * H_i);

		mat3da K(vab);
		kab.set(0, 0, K);

		vec3d N = (mp.dxr ^ mp.dxs);
		mat3d Kp = (N & dP);
		kab.add(0, 0, Kp);
	});
}
