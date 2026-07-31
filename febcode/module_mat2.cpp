#include "module_mat2.h"
using namespace febcode;

Value febcode::TransposeMat2(FuncArgs args)
{
	assert(args.count == 1);
	const mat2& m = args.getMat2();
	return mat2(
		m.m[0][0], m.m[1][0],
		m.m[0][1], m.m[1][1]
	);
}

Value febcode::InvertMat2(FuncArgs args)
{
	assert(args.count == 1);
	const mat2& m = args.getMat2();
	double a = m.m[0][0];
	double b = m.m[0][1];
	double c = m.m[1][0];
	double d = m.m[1][1];
	double det = a * d - b * c;
	double invDet = 1.0 / det;
	return mat2(
		d * invDet, -b * invDet,
		-c * invDet, a * invDet
	);
}
