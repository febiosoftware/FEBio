#pragma once
#include "fecore_enum.h"
#include "fecore_api.h"

class vec2d;
class vec3d;
class mat3d;

template <typename T> struct fecoreType {};

template <> struct fecoreType<double>
{
	static FEDataType type() { return FE_DOUBLE; }
	static int size() { return 1; }
};

template <> struct fecoreType<vec2d>
{
	static FEDataType type() { return FE_VEC2D; }
	static int size() { return 2; }
};

template <> struct fecoreType<vec3d>
{
	static FEDataType type() { return FE_VEC3D; }
	static int size() { return 3; }
};

template <> struct fecoreType<mat3d>
{
	static FEDataType type() { return FE_MAT3D; }
	static int size() { return 9; }
};

FECORE_API int fecore_data_size(FEDataType type);
