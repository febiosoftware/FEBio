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

	felog.printf("PID controller %d:\n", GetID());
	felog.printf("\tparameter = %lg\n", val);
	felog.printf("\terror     = %lg\n", error);
	felog.printf("\tvalue     = %lg\n", newVal);

	return newVal;
}
