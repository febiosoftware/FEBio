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
// Data initialization
void FETransverselyIsotropic::Init()
{
	FEUncoupledMaterial::Init();
	m_fib.Init();
}

//-----------------------------------------------------------------------------
// This material has only one property (the fiber material)
int FETransverselyIsotropic::Properties()
{
	return 1;
}

//-----------------------------------------------------------------------------
FEMaterial* FETransverselyIsotropic::GetProperty(int n)
{
	if ((n < 0) || (n > 1)) return 0;
	return &m_fib;
}

//-----------------------------------------------------------------------------
//! Serialize data to or from the dump file 
void FETransverselyIsotropic::Serialize(DumpFile &ar)
{
	// serialize the base class parameters
	FEUncoupledMaterial::Serialize(ar);

	if (ar.IsSaving())
	{
		ar << m_fib.m_ascl;
		ar << m_fib.m_ca0;
		ar << m_fib.m_beta;
		ar << m_fib.m_l0;
		ar << m_fib.m_refl;
	}
	else
	{
		ar >> m_fib.m_ascl;
		ar >> m_fib.m_ca0;
		ar >> m_fib.m_beta;
		ar >> m_fib.m_l0;
		ar >> m_fib.m_refl;
	}
}
