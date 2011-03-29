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
	m_nID = -1;
}

//-----------------------------------------------------------------------------
FEParameterList* FEMaterial::GetParameterList()
{
	FEParameterList* pl = new FEParameterList;
	return pl;
}

const char* FEMaterial::GetTypeString() { return "material base"; }

//-----------------------------------------------------------------------------
//! Store the material data to the archive
void FEMaterial::Serialize(DumpFile &ar)
{
	if (ar.IsSaving())
	{
		// store all parameters
		auto_ptr<FEParameterList> pl(GetParameterList());
		int n = pl->Parameters();
		ar << n;
		list<FEParam>::iterator it = pl->first();
		for (int j=0; j<n; ++j, ++it)
		{
			// store the value
			switch (it->m_itype)
			{
			case FE_PARAM_INT    : ar << it->value<int   >(); break;
			case FE_PARAM_BOOL   : ar << it->value<bool  >(); break;
			case FE_PARAM_DOUBLE : ar << it->value<double>(); break;
			case FE_PARAM_DOUBLEV: { for (int k=0; k<it->m_ndim; ++k) ar << it->pvalue<double>()[k]; } break;
			case FE_PARAM_INTV   : { for (int k=0; k<it->m_ndim; ++k) ar << it->pvalue<int   >()[k]; } break;
			default:
				assert(false);
			}

			// store parameter loadcurve data
			ar << it->m_nlc;
		}
	}
	else
	{
		auto_ptr<FEParameterList> pl(GetParameterList());
		int n = 0;
		ar >> n;
		assert(n == pl->Parameters());
		list<FEParam>::iterator it = pl->first();
		for (int j=0; j<n; ++j, ++it)
		{
			// read the value
			switch (it->m_itype)
			{
			case FE_PARAM_INT    : ar >> it->value<int   >(); break;
			case FE_PARAM_BOOL   : ar >> it->value<bool  >(); break;
			case FE_PARAM_DOUBLE : ar >> it->value<double>(); break;
			case FE_PARAM_DOUBLEV: { for (int k=0; k<it->m_ndim; ++k) ar >> it->pvalue<double>()[k]; } break;
			case FE_PARAM_INTV   : { for (int k=0; k<it->m_ndim; ++k) ar >> it->pvalue<int   >()[k]; } break;
			default:
				assert(false);
			}

			// read parameter data
			ar >> it->m_nlc;
		}
	}
}

//-----------------------------------------------------------------------------
// Derivative of stress w.r.t. solute concentration at material point
// Set this to zero by default because elasticity problems do not require it
mat3ds FESolidMaterial::Tangent_Concentration(FEMaterialPoint& pt)
{
	return mat3ds(0,0,0,0,0,0);
}

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

//-----------------------------------------------------------------------------
// Material parameters for FEHydraulicPermeability
BEGIN_PARAMETER_LIST(FEHydraulicPermeability, FEMaterial)
	ADD_PARAMETER(m_phi0, FE_PARAM_DOUBLE, "phi0");
END_PARAMETER_LIST();

void FEHydraulicPermeability::Init()
{
	FEMaterial::Init();
	
	if (!INRANGE(m_phi0, 0.0, 1.0)) throw MaterialError("phi0 must be in the range 0 < phi0 <= 1");
}

//-----------------------------------------------------------------------------
// Derivative of permeability w.r.t. solute concentration at material point
// Set this to zero by default because poroelasticity problems do not require it
mat3ds FEHydraulicPermeability::Tangent_Permeability_Concentration(FEMaterialPoint& pt)
{
	return mat3ds(0,0,0,0,0,0);
}
