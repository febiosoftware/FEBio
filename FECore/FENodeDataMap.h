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
#include "FEDataMap.h"

class FENodeSet;

class FECORE_API FENodeDataMap : public FEDataMap
{
public:
	FENodeDataMap();
	FENodeDataMap(FEDataType dataType);

	void Create(const FENodeSet* nodeSet, double val = 0.0);

	const FENodeSet* GetNodeSet() const;

	// return the item list associated with this map
	FEItemList* GetItemList() override;

	void Serialize(DumpStream& ar) override;

public:
	void setValue(int n, double v) override;
	void setValue(int n, const vec2d& v) override;
	void setValue(int n, const vec3d& v) override;
	void setValue(int n, const mat3d& v) override;
	void setValue(int n, const mat3ds& v) override;

	double getValue(int n) const;

	void fillValue(double v) override;
	void fillValue(const vec2d& v) override;
	void fillValue(const vec3d& v) override;
	void fillValue(const mat3d& v) override;
	void fillValue(const mat3ds& v) override;

	double value(const FEMaterialPoint& mp) override;
	vec3d valueVec3d(const FEMaterialPoint& mp) override;
	mat3d valueMat3d(const FEMaterialPoint& mp) override;
	mat3ds valueMat3ds(const FEMaterialPoint& mp) override;

private:
	const FENodeSet*	m_nodeSet;
};
