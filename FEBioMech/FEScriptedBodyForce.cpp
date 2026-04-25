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
#include "stdafx.h"
#include "FEScriptedBodyForce.h"

FEScriptedBodyForce::FEScriptedBodyForce(FEModel* pfem) : FEScripted<FEBodyForce>(pfem)
{
	ScriptContext sc;
	sc.returnType = FEValueType::Vec3d;
	sc.addVariable("pos0", FEValueType::Vec3d, false);
	sc.addVariable("pos" , FEValueType::Vec3d, true);
	sc.addVariable("time", FEValueType::Double, false);
	SetScriptContext(sc);
}

vec3d FEScriptedBodyForce::force(FEMaterialPoint& mp)
{
	std::vector<FEValue> var(3);
	var[0] = mp.m_r0;
	var[1] = mp.m_rt;
	var[2] = GetTimeInfo().currentTime;
	return Value(mp, var).v3;
}

mat3d FEScriptedBodyForce::stiffness(FEMaterialPoint& pt)
{
	if (HasDerivative(1))
	{
		std::vector<FEValue> var(3);
		var[0] = pt.m_r0;
		var[1] = pt.m_rt;
		var[2] = GetTimeInfo().currentTime;
		return -DerivValue(pt, var, 1).m3;
	}

	return mat3d(0.0);
}
