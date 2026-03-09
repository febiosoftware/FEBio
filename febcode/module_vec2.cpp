#include "module_vec2.h"
#include <assert.h>

using namespace febcode;

Value febcode::DotVec2(const Value* args, int argc)
{
	assert(argc == 2);
	const vec2& a = getVec2(args[0]);
	const vec2& b = getVec2(args[1]);
	return a*b;
}

Value febcode::NormalizeVec2(const Value* args, int argc)
{
	assert(argc == 1);
	const vec2& a = getVec2(args[0]);

	double x = a.x;
	double y = a.y;

	double D = std::sqrt(x * x + y * y);
	if (D != 0) {
		x /= D;
		y /= D;
	}

	return vec2(x, y);
}
