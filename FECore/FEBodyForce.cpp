#include "stdafx.h"
#include "FEBodyForce.h"

//-----------------------------------------------------------------------------
FEBodyForce::FEBodyForce(FEModel* pfem) : m_pfem(pfem)
{
	s[0] = s[1] = s[2] = 0.0;
	lc[0] = -1; lc[1] = -1; lc[2] = -1;
}

//-----------------------------------------------------------------------------
void FEBodyForce::Serialize(DumpFile& ar)
{
	if (ar.IsSaving())
	{
		ar << s[0] << s[1] << s[2];
		ar << lc[0] << lc[1] << lc[2];
	}
	else
	{
		ar >> s[0] >> s[1] >> s[2];
		ar >> lc[0] >> lc[1] >> lc[2];
	}
}
