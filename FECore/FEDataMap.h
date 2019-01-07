#pragma once
#include "FEDataArray.h"

class FEMaterialPoint;

//-----------------------------------------------------------------------------
// Base class for all data maps. A data map needs to be able to evaluate data across a domain
// TODO: This is a work in progress. 
// This class was added to create a base for FESurfaceMap and FEDomainMap so that both could be used in 
// FEMappedValue. 
class FECORE_API FEDataMap : public FEDataArray
{
public:
	FEDataMap(FEDataType dataType);
	FEDataMap(const FEDataMap& map);

	//! set the name
	void SetName(const std::string& name);

	//! get the name
	const std::string& GetName() const;

public:
	// This function needs to be overridden by derived classes
	virtual double value(const FEMaterialPoint& mp) = 0;
	virtual vec3d valueVec3d(const FEMaterialPoint& mp) = 0;
	virtual mat3d valueMat3d(const FEMaterialPoint& mp) = 0;

protected:
	std::string	m_name;					// name of data map TODO: Move to base class?
};
