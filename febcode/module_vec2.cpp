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

Value febcode::NormalizeVec2(FuncArgs args)
{
	assert(args.count == 1);
	const vec2& a = getVec2(args.getDouble());

	double x = a.x;
	double y = a.y;

	double D = std::sqrt(x * x + y * y);
	if (D != 0) {
		x /= D;
		y /= D;
	}

	return vec2(x, y);
}
