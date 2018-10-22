#pragma once
#include "FELoadController.h"
#include "MathObject.h"

class FEPIDController : public FELoadController
{
public:
	FEPIDController(FEModel* fem);

	bool Init() override;

protected:
	double GetValue(double time) override;

private:
	std::string		m_paramName;	// the parameter to target
	double			m_trg;			// the target value

	double	m_Kp;
	double	m_Kd;
	double	m_Ki;

	FEParamValue	m_param;
	double			m_prev;
	double			m_prevTime;
	double			m_int;

	DECLARE_FECORE_CLASS();
};
