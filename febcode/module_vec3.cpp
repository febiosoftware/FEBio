#include "module_vec3.h"

using namespace febcode;

Value febcode::DotVec3(const std::vector<Value>& args)
{
	const vec3& a = getVec3(args[0]);
	const vec3& b = getVec3(args[1]);
	return a*b;
}

Value febcode::CrossVec3(const std::vector<Value>& args)
{
	const vec3& a = getVec3(args[0]);
	const vec3& b = getVec3(args[1]);
	return a.cross(b);
}

Value febcode::NormalizeVec3(const std::vector<Value>& args)
{
	const vec3& a = getVec3(args[0]);

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
