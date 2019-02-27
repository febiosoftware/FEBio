#include "stdafx.h"
#include "FEGlobalData.h"

REGISTER_SUPER_CLASS(FEGlobalData, FEGLOBALDATA_ID);

//-----------------------------------------------------------------------------
FEGlobalData::FEGlobalData(FEModel* fem) : FECoreBase(fem)
{
}

//-----------------------------------------------------------------------------
bool FEGlobalData::Init()
{
	return true;
}
