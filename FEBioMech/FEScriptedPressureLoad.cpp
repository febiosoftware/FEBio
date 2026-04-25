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

FEScriptedPressureLoad::FEScriptedPressureLoad(FEModel* pfem) : FEScripted<FESurfaceLoad>(pfem)
{
	ScriptContext sc;
	sc.returnType = FEValueType::Double;
	sc.addVariable("pos"   , FEValueType::Vec3d, true);
	sc.addVariable("normal", FEValueType::Vec3d, true);
	sc.addVariable("time"  , FEValueType::Double, false);
	SetScriptContext(sc);
}

bool FEScriptedPressureLoad::Init()
{
	m_dof.AddVariable(FEBioMech::GetVariableName(FEBioMech::DISPLACEMENT));
	return FESurfaceLoad::Init();
}

void FEScriptedPressureLoad::LoadVector(FEGlobalVector& R)
{
	double t = GetTimeInfo().currentTime;

	// evaluate the integral
	FESurface& surf = GetSurface();
	surf.LoadVector(R, m_dof, false, [&](FESurfaceMaterialPoint& pt, const FESurfaceDofShape& dof_a, std::vector<double>& val) {

		vec3d normal = (pt.dxr ^ pt.dxs).normalized();

		std::vector<FEValue> vars(3);
		vars[0] = pt.m_rt;
		vars[1] = normal;
		vars[2] = t;

		// evaluate pressure at this material point
		double P = -Value(pt, vars).d;

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
	surf.LoadStiffness(LS, m_dof, m_dof, [&](FESurfaceMaterialPoint& mp, const FESurfaceDofShape& dof_a, const FESurfaceDofShape& dof_b, matrix& Kab) {

		vec3d normal = (mp.dxr ^ mp.dxs).normalized();

		std::vector<FEValue> vars(3);
		vars[0] = mp.m_rt;
		vars[1] = normal;
		vars[2] = t;

		double H_i = dof_a.shape;
		double H_j = dof_b.shape;

		double Gr_j = dof_b.shape_deriv_r;
		double Gs_j = dof_b.shape_deriv_s;

		mat3da Gr(mp.dxr);
		mat3da Gs(mp.dxs);
		mat3d Grs_j = Gr * Gs_j - Gs * Gr_j;

		// evaluate pressure at this material point
		double P = Value(mp, vars).d;

		mat3d K = Grs_j * (P * H_i);
		Kab.set(0, 0, K);

		if (HasDerivative(0))
		{
			vec3d dPx = DerivValue(mp, vars, 0).v3;
			vec3d N = (mp.dxr ^ mp.dxs);
			K = (N & dPx) * (H_i * H_j);
			Kab.add(0, 0, K);
		}

		if (HasDerivative(1))
		{
			mat3dd I(1.0);
			mat3d nxn = (normal & normal);

			vec3d dPn = DerivValue(mp, vars, 1).v3;
			K = Grs_j * ((normal & dPn)* (I - nxn) * H_i);
			Kab.add(0, 0, K);
		}
	});
}
