#pragma once
#include "FEDataArray.h"

class FEMaterialPoint;

//-----------------------------------------------------------------------------
// Base class for all data maps. A data map needs to be able to evaluate data across a domain
// TODO: This is a work in progress. 
// This class was added to create a base for FESurfaceMap and FEDomainMap so that both could be used in 
// FEMappedValue. 
class FEDataMap : public FEDataArray
{
public:
	FEDataMap(int dataSize) : FEDataArray(dataSize) {}
	FEDataMap(const FEDataMap& map) : FEDataArray(map) {}

	// This function needs to be overridden by derived classes
	virtual double value(const FEMaterialPoint& mp) = 0;
};
