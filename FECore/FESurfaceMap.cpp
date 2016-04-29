#include "stdafx.h"
#include "FESurfaceMap.h"
#include "FESurface.h"
#include "DumpStream.h"

//-----------------------------------------------------------------------------
FESurfaceMap::FESurfaceMap()
{
}

//-----------------------------------------------------------------------------
bool FESurfaceMap::Create(const FESurface* ps, double def)
{
	int NF = ps->Elements();
	m_val.resize(NF, def);
	return true;
}

//-----------------------------------------------------------------------------
bool FESurfaceMap::SetValue(const FEFacetIndex& n, double v)
{
	m_val[n] = v;
	return true;
}

//-----------------------------------------------------------------------------
void FESurfaceMap::Serialize(DumpStream& ar)
{
	if (ar.IsSaving())
	{
		ar << m_val;
	}
	else
	{
		ar >> m_val;
	}
}
