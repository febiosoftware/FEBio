#include "stdafx.h"
#include "FEDataMap.h"
#include "fecore_type.h"

//-----------------------------------------------------------------------------
FEDataMap::FEDataMap(FEDataMapType mapType, FEDataType dataType) : FEDataArray(mapType, dataType) {}

//-----------------------------------------------------------------------------
FEDataMap::FEDataMap(const FEDataMap& map) : FEDataArray(map) {}

//-----------------------------------------------------------------------------
void FEDataMap::SetName(const std::string& name)
{
	m_name = name;
}

//-----------------------------------------------------------------------------
const std::string& FEDataMap::GetName() const { return m_name; }
