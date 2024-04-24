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
#include "FEMathController.h"
#include "FEModel.h"

BEGIN_FECORE_CLASS(FEMathController, FELoadController)
	ADD_PARAMETER(m_var, "var");
	ADD_PARAMETER(m_math, "math");
END_FECORE_CLASS();

FEMathController::FEMathController(FEModel* fem) : FELoadController(fem)
{
}

bool FEMathController::Init()
{
	FEModel& fem = *GetFEModel();
	m_val.AddVariable("t");
	char sz[64] = { 0 };
	for (int i = 0; i < (int)m_var.size(); ++i)
	{
		ParamString ps(m_var[i].c_str());
		FEParamValue param = fem.GetParameterValue(ps);
		if (param.isValid() == false) return false;
		if (param.type() != FE_PARAM_DOUBLE) return false;

		m_param.push_back(param);

		snprintf(sz, 64, "var_%d", i);
		m_val.AddVariable(sz);
	}

	if (m_val.Create(m_math) == false) return false;

	return FELoadController::Init();
}

double FEMathController::GetValue(double time)
{
	vector<double> p(1 + m_param.size());
	p[0] = time;
	for (int i = 0; i < m_param.size(); ++i) p[1 + i] = m_param[i].value<double>();
	return m_val.value_s(p);
}
