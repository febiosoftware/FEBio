#include "stdafx.h"
#include "FEGlobalData.h"

//-----------------------------------------------------------------------------
FEGlobalData::FEGlobalData(FEModel* fem) : FECoreBase(fem, FEGLOBALDATA_ID)
{
}

//-----------------------------------------------------------------------------
bool FEGlobalData::Init()
{
	return true;
}
