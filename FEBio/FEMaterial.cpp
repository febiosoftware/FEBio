// FEMaterial.cpp: implementation of the FEMaterial class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include <math.h>
#include <stdarg.h>
#include "FEMaterial.h"

//-----------------------------------------------------------------------------
MaterialError::MaterialError(const char* szfmt, ...)
{
	// get a pointer to the argument list
	va_list	args;

	// print to file
	va_start(args, szfmt);
	vsprintf(m_szerr, szfmt, args);
	va_end(args);
}

//-----------------------------------------------------------------------------
FEParameterList* FEMaterial::GetParameterList()
{
	FEParameterList* pl = new FEParameterList;
	return pl;
}

const char* FEMaterial::GetTypeString() { return "material base"; }

//-----------------------------------------------------------------------------
// Material parameters for FEElasticMaterial
BEGIN_PARAMETER_LIST(FEElasticMaterial, FEMaterial)
	ADD_PARAMETER(m_density, FE_PARAM_DOUBLE, "density");
END_PARAMETER_LIST();

void FEElasticMaterial::Init()
{
	FEMaterial::Init();

	if (m_density <= 0) throw MaterialError("Invalid material density");
}

//-----------------------------------------------------------------------------
// Material parameters for FEIncompressibleMaterial
BEGIN_PARAMETER_LIST(FEIncompressibleMaterial, FEElasticMaterial)
	ADD_PARAMETER(m_K, FE_PARAM_DOUBLE, "k");
	ADD_PARAMETER(m_blaugon, FE_PARAM_BOOL  , "laugon");
	ADD_PARAMETER(m_atol   , FE_PARAM_DOUBLE, "atol"  );
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
// Material parameters for FETransverselyIsotropic
BEGIN_PARAMETER_LIST(FETransverselyIsotropic, FEIncompressibleMaterial)
	ADD_PARAMETER(c3, FE_PARAM_DOUBLE, "c3");
	ADD_PARAMETER(c4, FE_PARAM_DOUBLE, "c4");
	ADD_PARAMETER(c5, FE_PARAM_DOUBLE, "c5");
	ADD_PARAMETER(lam1, FE_PARAM_DOUBLE, "lam_max");
END_PARAMETER_LIST();
