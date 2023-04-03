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
#include "FEPIDController.h"
#include "FEModel.h"
#include "log.h"

BEGIN_FECORE_CLASS(FEPIDController, FELoadController)
	ADD_PARAMETER(m_paramName, "var");
	ADD_PARAMETER(m_trg, "target");
	ADD_PARAMETER(m_Kp, "Kp");
	ADD_PARAMETER(m_Ki, "Ki");
	ADD_PARAMETER(m_Kd, "Kd");
END_FECORE_CLASS();

FEPIDController::FEPIDController(FEModel* fem) : FELoadController(fem)
{
	m_prev = 0;
	m_int = 0;
	m_prevTime = 0.0;
}

bool FEPIDController::Init()
{
	FEModel& fem = *GetFEModel();

	ParamString ps(m_paramName.c_str());
	m_param = fem.GetParameterValue(ps);
	if (m_param.isValid() == false) return false;
	if (m_param.type() != FE_PARAM_DOUBLE) return false;

	return FELoadController::Init();
}

double FEPIDController::GetValue(double time)
{
	double val = m_param.value<double>();
	double error = m_trg - val;

	double newVal = m_Kp*error;

	double dt = time - m_prevTime;
	if (dt != 0.0)
	{
		double der = (error - m_prev) / dt;
		m_int += error*dt;
		newVal += m_Kd*der + m_Ki*m_int;
	}

	m_prev = error;
	m_prevTime = time;

	if (GetFEModel()->GetPrintParametersFlag())
	{
		feLog("PID controller %d:\n", GetID());
		feLog("\tparameter = %lg\n", val);
		feLog("\terror     = %lg\n", error);
		feLog("\tvalue     = %lg\n", newVal);
	}

	return newVal;
}
