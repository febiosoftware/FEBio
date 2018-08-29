#include "stdafx.h"
#include "FEGlobalData.h"

//-----------------------------------------------------------------------------
FEGlobalData::FEGlobalData(FEModel* pfem) : FECoreBase(FEGLOBALDATA_ID)
{
	m_pfem = pfem;
}

//-----------------------------------------------------------------------------
bool FEGlobalData::Init()
{
	return true;
}
