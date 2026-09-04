#include "module_mat2.h"
using namespace febcode;

Value febcode::TransposeMat2(FuncArgs args)
{
	assert(args.count == 1);
	const mat2& m = args.getMat2();
	return m.transpose();
}

Value febcode::InvertMat2(FuncArgs args)
{
	assert(args.count == 1);
	const mat2& m = args.getMat2();
	double a = m(0, 0);
	double b = m(0, 1);
	double c = m(1, 0);
	double d = m(1, 1);
	double det = a * d - b * c;
	double invDet = 1.0 / det;
	return mat2(
		d * invDet, -b * invDet,
		-c * invDet, a * invDet
	);
}

Value febcode::TraceMat2(FuncArgs args)
{
	assert(args.count == 1);
	const mat2& m = args.getMat2();
	return m.trace();
}

Value febcode::DeterminantMat2(FuncArgs args)
{
	assert(args.count == 1);
	const mat2& m = args.getMat2();
	return m.det();
}
