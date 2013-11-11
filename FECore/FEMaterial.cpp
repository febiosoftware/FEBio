// FEMaterial.cpp: implementation of the FEMaterial class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FEMaterial.h"
#include <math.h>
#include <stdarg.h>
#include "FEModel.h"

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
	m_pmap = 0;
}

//-----------------------------------------------------------------------------
FEMaterial::~FEMaterial()
{
	if (m_pmap) delete m_pmap; 
}

//-----------------------------------------------------------------------------
//! Sets the name of the material
void FEMaterial::SetName(const char* sz)
{ 
	strcpy(m_szname, sz); 
}

//-----------------------------------------------------------------------------
//! Return the name of the material
const char* FEMaterial::GetName()
{ 
	return m_szname; 
}

//-----------------------------------------------------------------------------
void FEMaterial::SetCoordinateSystemMap(FECoordSysMap* pmap)
{
	m_pmap = pmap;
}

//-----------------------------------------------------------------------------
FECoordSysMap* FEMaterial::GetCoordinateSystemMap()
{
	return m_pmap;
}

//-----------------------------------------------------------------------------
//! This function doesn't initialize anything but should be overridden in 
//! derived classes to initialize and check material parameters
void FEMaterial::Init()
{
}

//-----------------------------------------------------------------------------
//! Store the material data to the archive
void FEMaterial::Serialize(DumpFile &ar)
{
	if (ar.IsSaving())
	{
		ar << m_nID;

		// store all parameters
		FEParameterList& pl = GetParameterList();
		int n = pl.Parameters();
		ar << n;
		list<FEParam>::iterator it = pl.first();
		for (int j=0; j<n; ++j, ++it)
		{
			// store the value
			switch (it->m_itype)
			{
			case FE_PARAM_INT    : ar << it->value<int   >(); break;
			case FE_PARAM_BOOL   : ar << it->value<bool  >(); break;
			case FE_PARAM_DOUBLE : ar << it->value<double>(); break;
			case FE_PARAM_VEC3D  : ar << it->value<vec3d >(); break;
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
		ar >> m_nID;

		FEParameterList& pl = GetParameterList();
		int n = 0;
		ar >> n;
		assert(n == pl.Parameters());
		list<FEParam>::iterator it = pl.first();
		for (int j=0; j<n; ++j, ++it)
		{
			// read the value
			switch (it->m_itype)
			{
			case FE_PARAM_INT    : ar >> it->value<int   >(); break;
			case FE_PARAM_BOOL   : ar >> it->value<bool  >(); break;
			case FE_PARAM_DOUBLE : ar >> it->value<double>(); break;
			case FE_PARAM_VEC3D  : ar >> it->value<vec3d >(); break;
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

//=============================================================================
// FEMultiMaterial
//=============================================================================

//-----------------------------------------------------------------------------
int FEMultiMaterial::FindComponent(const char* sz, int nid)
{
	for (int i=0; i<(int) m_Mat.size(); ++i)
	{
		Property* p = m_Mat[i];
		if ((strcmp(p->GetName(), sz) == 0) && (p->GetID() == nid)) return i;
	}
	return -1;
}
