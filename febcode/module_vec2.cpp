#include "module_vec2.h"
#include <assert.h>

using namespace febcode;

Value febcode::DotVec2(FuncArgs args)
{
	assert(args.count == 2);
	vec2 a = args.getVec2();
	vec2 b = args.getVec2();
	return a*b;
}

Value febcode::LengthVec2(FuncArgs args)
{
	assert(args.count == 1);
	const vec2& a = args.getVec2();
	return std::sqrt(a.x*a.x + a.y*a.y);
}

Value febcode::NormalizeVec2(FuncArgs args)
{
	assert(args.count == 1);
	const vec2& a = args.getVec2();

	double x = a.x;
	double y = a.y;

	double D = std::sqrt(x * x + y * y);
	if (D != 0) {
		x /= D;
		y /= D;
	}

	return vec2(x, y);
}

Value febcode::OuterVec2(FuncArgs args)
{
	assert(args.count == 2);
	const vec2& a = args.getVec2();
	const vec2& b = args.getVec2();
	return mat2(
		a.x * b.x, a.x * b.y,
		a.y * b.x, a.y * b.y
	);
}
