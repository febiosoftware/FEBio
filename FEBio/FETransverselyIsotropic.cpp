#include "stdafx.h"
#include "FETransverselyIsotropic.h"

//-----------------------------------------------------------------------------
// Material parameters for FETransverselyIsotropic
BEGIN_PARAMETER_LIST(FETransverselyIsotropic, FEUncoupledMaterial)
	ADD_PARAMETER(m_fib.m_c3, FE_PARAM_DOUBLE, "c3");
	ADD_PARAMETER(m_fib.m_c4, FE_PARAM_DOUBLE, "c4");
	ADD_PARAMETER(m_fib.m_c5, FE_PARAM_DOUBLE, "c5");
	ADD_PARAMETER(m_fib.m_lam1, FE_PARAM_DOUBLE, "lam_max");
END_PARAMETER_LIST();
