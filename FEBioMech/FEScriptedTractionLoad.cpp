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
#include "FEScriptedTractionLoad.h"

BEGIN_FECORE_CLASS(FEScriptedTractionLoad, FESurfaceLoad)
	ADD_PARAMETER(m_scriptName, "script");
END_FECORE_CLASS()

FEScriptedTractionLoad::FEScriptedTractionLoad(FEModel* pfem) : FESurfaceLoad(pfem), FEScriptedBehavior(pfem)
{
}

bool FEScriptedTractionLoad::Init()
{
	m_dof.AddVariable(FEBioMech::GetVariableName(FEBioMech::DISPLACEMENT));

	SetSibling(this);
	AddVariable("pos", FEValueType::Vec3d);
	AddVariable("normal", FEValueType::Vec3d);
	AddVariable("time", FEValueType::Double, false);

	SetProgramReturnType(FEValueType::Vec3d);

	if (FESurfaceLoad::Init() == false) return false;
	if (FEScriptedBehavior::Init() == false) return false;
	return true;
}

void FEScriptedTractionLoad::LoadVector(FEGlobalVector& R)
{
	double t = GetTimeInfo().currentTime;

	// evaluate the integral
	FESurface& surf = GetSurface();
	surf.LoadVector(R, m_dof, false, [&](FESurfaceMaterialPoint& pt, const FESurfaceDofShape& dof_a, std::vector<double>& val) {

		vec3d N = (pt.dxr ^ pt.dxs);
		double J = N.unit();

		std::vector<FEValue> vars(3);
		vars[0] = pt.m_rt;
		vars[1] = N;
		vars[2] = t;

		// evaluate traction at this material point
		vec3d t = Value(pt, vars).v3;

		double H_u = dof_a.shape;

		val[0] = H_u * t.x * J;
		val[1] = H_u * t.y * J;
		val[2] = H_u * t.z * J;
	});
}

void FEScriptedTractionLoad::StiffnessMatrix(FELinearSystem& LS)
{
	double t = GetTimeInfo().currentTime;

	// evaluate the integral
	FESurface& surf = GetSurface();
	surf.LoadStiffness(LS, m_dof, m_dof, [&](FESurfaceMaterialPoint& mp, const FESurfaceDofShape& dof_a, const FESurfaceDofShape& dof_b, matrix& Kab) {

		vec3d N = (mp.dxr ^ mp.dxs);
		double J = N.unit();

		std::vector<FEValue> vars(3);
		vars[0] = mp.m_rt;
		vars[1] = N;
		vars[2] = t;

		double H_i = dof_a.shape;
		double H_j = dof_b.shape;

		double Gr_j = dof_b.shape_deriv_r;
		double Gs_j = dof_b.shape_deriv_s;

		mat3da Gr(mp.dxr);
		mat3da Gs(mp.dxs);
		mat3d Grs_j = Gr * Gs_j - Gs * Gr_j;

		// evaluate traction at this material point
		vec3d t = -Value(mp, vars).v3;
		mat3d K = (t & N)*Grs_j*H_i;
		Kab.set(0, 0, K);

		// evaluate traction gradient w.r.t. position
		if (HasDerivative(0))
		{
			mat3d dtdx = -DerivValue(mp, vars, 0).m3;
			K = dtdx * (H_i * H_j * J);
			Kab.add(0, 0, K);
		}

		// evaluate traction gradient w.r.t. normal
		if (HasDerivative(1))
		{
			mat3d dtdn = -DerivValue(mp, vars, 1).m3;
			mat3dd I(1.0);
			K = dtdn * (I - (N & N)) * Grs_j * H_i;
			Kab.add(0, 0, K);
		}
	});
}
