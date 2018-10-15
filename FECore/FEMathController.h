#pragma once
#include "FELoadController.h"
#include "MathObject.h"
#include <vector>
#include <string>

class FEMathController : public FELoadController
{
public:
	FEMathController(FEModel* fem);

	bool Init() override;

protected:
	double GetValue(double time) override;

private:
	std::vector<std::string>	m_var;
	std::string					m_math;

	MSimpleExpression			m_val;
	std::vector<FEParamValue>	m_param;

	DECLARE_FECORE_CLASS();
};
