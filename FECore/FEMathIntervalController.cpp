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
#include "FEMathIntervalController.h"
#include "FEModel.h"

BEGIN_FECORE_CLASS(FEMathIntervalController, FELoadController)
	ADD_PARAMETER(m_rng, 2, "interval");
	ADD_PARAMETER(m_leftExtend , "left_extend", 0, "zero\0constant\0repeat\0");
	ADD_PARAMETER(m_rightExtend, "right_extend", 0, "zero\0constant\0repeat\0");
	ADD_PARAMETER(m_var, "var");
	ADD_PARAMETER(m_math, "math");
END_FECORE_CLASS();

FEMathIntervalController::FEMathIntervalController(FEModel* fem) : FELoadController(fem)
{
	m_rng[0] = 0.0;
	m_rng[1] = 1.0;
	m_leftExtend = ExtendMode::CONSTANT;
	m_rightExtend = ExtendMode::CONSTANT;
}

bool FEMathIntervalController::Init()
{
	double Dt = m_rng[1] - m_rng[0];
	if (Dt <= 0.0) return false;

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

		sprintf(sz, "var_%d", i);
		m_val.AddVariable(sz);
	}

	if (m_val.Create(m_math) == false) return false;

	return FELoadController::Init();
}

double FEMathIntervalController::GetValue(double time)
{
	if (time <= m_rng[0])
	{
		switch (m_leftExtend)
		{
		case ExtendMode::ZERO: return 0.0; break;
		case ExtendMode::CONSTANT: time = m_rng[0]; break;
		case ExtendMode::REPEAT:
		{
			double Dt = m_rng[1] - m_rng[0];
			time = m_rng[1] - fmod(m_rng[0] - time, Dt);
		}
		break;
		}
	}
	else if (time >= m_rng[1])
	{
		switch (m_rightExtend)
		{
		case ExtendMode::ZERO: return 0.0; break;
		case ExtendMode::CONSTANT: time = m_rng[1]; break;
		case ExtendMode::REPEAT:
		{
			double Dt = m_rng[1] - m_rng[0];
			time = m_rng[0] + fmod(time - m_rng[0], Dt);
		}
		break;
		}
	}

	vector<double> p(1 + m_param.size());
	p[0] = time;
	for (int i = 0; i < m_param.size(); ++i) p[1 + i] = m_param[i].value<double>();
	return m_val.value_s(p);
}
