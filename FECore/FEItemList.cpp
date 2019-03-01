#include "stdafx.h"
#include "FEItemList.h"
#include "DumpStream.h"

void FEItemList::Serialize(DumpStream& ar)
{
	if (ar.IsShallow()) return;
	ar & m_name;
}
