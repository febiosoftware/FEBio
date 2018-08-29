#include "stdafx.h"
#include "FESoluteInterface.h"
#include "FESolute.h"

//-----------------------------------------------------------------------------
//! Returns the local solute index given the global ID
int FESoluteInterface::FindLocalSoluteID(int nid)
{
	int lsid = -1;
	for (int isol = 0; isol<Solutes(); ++isol)
		if (GetSolute(isol)->GetSoluteID() == nid) {
			lsid = isol;
			break;
		}
	return lsid;
}
