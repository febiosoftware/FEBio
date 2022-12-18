/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/



#pragma once
#include "FEDataArray.h"

class FEMaterialPoint;
class FEItemList;

//-----------------------------------------------------------------------------
// Base class for all data maps. A data map needs to be able to evaluate data across a domain
// TODO: This is a work in progress. 
// This class was added to create a base for FESurfaceMap and FEDomainMap so that both could be used in 
// FEMappedValue. 
class FECORE_API FEDataMap : public FEDataArray
{
public:
	FEDataMap(FEDataMapType mapType, FEDataType dataType = FE_INVALID_TYPE);
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
	virtual mat3ds valueMat3ds(const FEMaterialPoint& mp) = 0;

	// return the item list associated with this map
	virtual FEItemList* GetItemList() = 0;

public:
	void Serialize(DumpStream& ar) override;

protected:
	std::string	m_name;					// name of data map TODO: Move to base class?
};
