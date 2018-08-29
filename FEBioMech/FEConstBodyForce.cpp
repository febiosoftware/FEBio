#include "stdafx.h"
#include "FEConstBodyForce.h"

BEGIN_PARAMETER_LIST(FEConstBodyForce, FEBodyForce);
	ADD_PARAMETER(m_f.x, FE_PARAM_DOUBLE, "x");
	ADD_PARAMETER(m_f.y, FE_PARAM_DOUBLE, "y");
	ADD_PARAMETER(m_f.z, FE_PARAM_DOUBLE, "z");
END_PARAMETER_LIST();
