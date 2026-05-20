#include "module_mat3.h"
using namespace febcode;

Value febcode::TransposeMat3(FuncArgs args)
{
	assert(args.count == 1);
	const mat3& m = args.getMat3();
	return m.transpose();
}

Value febcode::InvertMat3(FuncArgs args)
{
	assert(args.count == 1);
	const mat3& m = args.getMat3();
	return m.inverse();
}

Value febcode::OuterVec3(FuncArgs args)
{
	assert(args.count == 2);
	const vec3& a = args.getVec3();
	const vec3& b = args.getVec3();
	return mat3(
		a.x * b.x, a.x * b.y, a.x * b.z,
		a.y * b.x, a.y * b.y, a.y * b.z,
		a.z * b.x, a.z * b.y, a.z * b.z
	);
}
