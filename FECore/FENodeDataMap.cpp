/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, Columbia University, and others.

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

#include "stdafx.h"
#include "FENodeDataMap.h"

FENodeDataMap::FENodeDataMap() : FEDataArray(FE_NODE_DATA_MAP, FE_INVALID_TYPE)
{

}

FENodeDataMap::FENodeDataMap(FEDataType dataType) : FEDataArray(FE_NODE_DATA_MAP, dataType)
{
	
}

void FENodeDataMap::Create(int nsize, double val)
{
	resize(nsize, val);
}

void FENodeDataMap::Add(double val)
{
	push_back(val);
}

double FENodeDataMap::getValue(int n) const
{
	return get<double>(n);
}

void FENodeDataMap::setValue(int n, double v)
{
	set<double>(n, v);
}

void FENodeDataMap::setValue(int n, const vec2d& v)
{
	set<vec2d>(n, v);
}

void FENodeDataMap::setValue(int n, const vec3d& v)
{
	set<vec3d>(n, v);
}

void FENodeDataMap::setValue(int n, const mat3d& v)
{
	set<mat3d>(n, v);
}

void FENodeDataMap::fillValue(double v)
{
	set<double>(v);
}

void FENodeDataMap::fillValue(const vec2d& v)
{
	set<vec2d>(v);
}

void FENodeDataMap::fillValue(const vec3d& v)
{
	set<vec3d>(v);
}

void FENodeDataMap::fillValue(const mat3d& v)
{
	set<mat3d>(v);
}
