#include "module_vec3.h"

using namespace febcode;

Value febcode::DotVec3(FuncArgs args)
{
	assert(args.count == 2);
	const vec3& a = args.getVec3();
	const vec3& b = args.getVec3();
	return a*b;
}

Value febcode::CrossVec3(FuncArgs args)
{
	assert(args.count == 2);
	const vec3& a = args.getVec3();
	const vec3& b = args.getVec3();
	return a.cross(b);
}

Value febcode::NormalizeVec3(FuncArgs args)
{
	assert(args.count == 1);
	const vec3& a = args.getVec3();
	double x = a.x;
	double y = a.y;
	double z = a.z;

	double D = std::sqrt(x * x + y * y + z * z);
	if (D != 0) {
		x /= D;
		y /= D;
		z /= D;
	}

	return vec3(x,y,z);
}

Value febcode::LengthVec3(FuncArgs args)
{
	assert(args.count == 1);
	const vec3& a = args.getVec3();
	return std::sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
}
