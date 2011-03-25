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


//-----------------------------------------------------------------------------
//! Serialize data to or from the dump file 
void FETransverselyIsotropic::Serialize(DumpFile &ar)
{
	// serialize the base class parameters
	FEUncoupledMaterial::Serialize(ar);

	if (ar.IsSaving())
	{
		ar << m_fib.m_lcna;
		ar << m_fib.m_ascl;
		ar << m_fib.m_ca0;
		ar << m_fib.m_beta;
		ar << m_fib.m_l0;
		ar << m_fib.m_refl;
	}
	else
	{
		ar >> m_fib.m_lcna;
		ar >> m_fib.m_ascl;
		ar >> m_fib.m_ca0;
		ar >> m_fib.m_beta;
		ar >> m_fib.m_l0;
		ar >> m_fib.m_refl;

//		if (pm->m_fib.m_lcna >= 0) pm->m_fib.m_plc = GetLoadCurve(pm->m_fib.m_lcna);
	}
}
