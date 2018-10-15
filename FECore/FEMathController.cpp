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

		sprintf(sz, "var_%d", i);
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
