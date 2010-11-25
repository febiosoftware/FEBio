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
FEMaterial::FEMaterial()
{
	static int n = 1;
	m_szname[0] = 0;
	m_nID = n++;
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
// Material parameters for the FENestedMaterial
BEGIN_PARAMETER_LIST(FENestedMaterial, FEMaterial)
	ADD_PARAMETER(m_nBaseMat, FE_PARAM_INT, "solid_id");
END_PARAMETER_LIST();
