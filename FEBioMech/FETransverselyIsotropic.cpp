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
//! find a material property index ( returns <0 for error)
int FETransverselyIsotropic::FindPropertyIndex(const char* szname)
{
	if (strcmp(szname, "active_contraction") == 0) return 1;
	return -1;
}

//-----------------------------------------------------------------------------
//! set a material property (returns false on error)
bool FETransverselyIsotropic::SetProperty(int i, FEMaterial* pm)
{
	if (i == 1)
	{
		FEActiveFiberContraction* pma = dynamic_cast<FEActiveFiberContraction*>(pm);
		if (pma) { m_fib.SetActiveContraction(pma); return true; }
	}
	return false;
}

//-----------------------------------------------------------------------------
//! Serialize data to or from the dump file 
void FETransverselyIsotropic::Serialize(DumpFile &ar)
{
	// serialize the base class parameters
	FEUncoupledMaterial::Serialize(ar);

	// serialize fiber data
	m_fib.Serialize(ar);
}
