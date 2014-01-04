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
FETransverselyIsotropic::FETransverselyIsotropic(FEModel* pfem) : FEUncoupledMaterial(pfem), m_fib(pfem) 
{
}

//-----------------------------------------------------------------------------
// Data initialization
void FETransverselyIsotropic::Init()
{
	FEUncoupledMaterial::Init();
	m_fib.Init();
}

//-----------------------------------------------------------------------------
// This material has two properties (the fiber material and the active contraction material)
int FETransverselyIsotropic::Properties()
{
	return 2;
}

//-----------------------------------------------------------------------------
FECoreBase* FETransverselyIsotropic::GetProperty(int n)
{
	if (n == 0) return &m_fib;
	if (n == 1) return m_fib.GetActiveContraction();
	return 0;
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
bool FETransverselyIsotropic::SetProperty(int i, FECoreBase* pm)
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
