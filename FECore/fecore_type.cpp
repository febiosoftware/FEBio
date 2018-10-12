#include "stdafx.h"
#include "fecore_type.h"
#include <assert.h>

int fecore_data_size(FEDataType type)
{
	switch (type)
	{
	case FE_DOUBLE: return fecoreType<double>::size(); break;
	case FE_VEC2D : return fecoreType<vec2d >::size(); break;
	case FE_VEC3D : return fecoreType<vec3d >::size(); break;
	case FE_MAT3D : return fecoreType<mat3d >::size(); break;
	default:
		assert(false);
	}

	return 0;
};
