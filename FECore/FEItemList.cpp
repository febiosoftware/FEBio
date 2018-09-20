#include "stdafx.h"
#include "FEItemList.h"

void FEItemList::Serialize(DumpStream& ar)
{
	if (ar.IsShallow() == false)
	{
		if (ar.IsSaving()) ar << m_name;
		else ar >> m_name;
	}
}
