#include "stdafx.h"
#include "FEDataMap.h"

//-----------------------------------------------------------------------------
FEDataMap::FEDataMap(int dataSize) : FEDataArray(dataSize) {}

//-----------------------------------------------------------------------------
FEDataMap::FEDataMap(const FEDataMap& map) : FEDataArray(map) {}

//-----------------------------------------------------------------------------
void FEDataMap::SetName(const std::string& name)
{
	m_name = name;
}

//-----------------------------------------------------------------------------
const std::string& FEDataMap::GetName() const { return m_name; }
