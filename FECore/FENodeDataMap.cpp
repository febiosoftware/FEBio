#include "stdafx.h"
#include "FENodeDataMap.h"

FENodeDataMap::FENodeDataMap(FEDataType dataType) : FEDataArray(dataType)
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
