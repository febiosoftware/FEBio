#include "module_mat3.h"
using namespace febcode;

Value febcode::TransposeMat3(FuncArgs args)
{
	assert(args.count == 1);
	const mat3& m = args.getMat3();
	return mat3(
		m.m[0][0], m.m[1][0], m.m[2][0],
		m.m[0][1], m.m[1][1], m.m[2][1],
		m.m[0][2], m.m[1][2], m.m[2][2]
	);
}

Value febcode::InvertMat3(FuncArgs args)
{
	assert(args.count == 1);
	const mat3& m = args.getMat3();
	const double (*a)[3] = m.m; // for easier access to elements

	double det = a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) -
				 a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0]) +
				 a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);

	mat3 mInv;
	double (*b)[3] = mInv.m; // for easier access to elements
	b[0][0] = (a[1][1] * a[2][2] - a[1][2] * a[2][1]) / det;
	b[0][1] = (a[0][2] * a[2][1] - a[0][1] * a[2][2]) / det;
	b[0][2] = (a[0][1] * a[1][2] - a[0][2] * a[1][1]) / det;

	b[1][0] = (a[1][2] * a[2][0] - a[1][0] * a[2][2]) / det;
	b[1][1] = (a[0][0] * a[2][2] - a[0][2] * a[2][0]) / det;
	b[1][2] = (a[0][2] * a[1][0] - a[0][0] * a[1][2]) / det;

	b[2][0] = (a[1][0] * a[2][1] - a[1][1] * a[2][0]) / det;
	b[2][1] = (a[0][1] * a[2][0] - a[0][0] * a[2][1]) / det;
	b[2][2] = (a[0][0] * a[1][1] - a[0][1] * a[1][0]) / det;

	return mInv;
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
