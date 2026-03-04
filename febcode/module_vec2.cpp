#include "module_vec2.h"

using namespace febcode;

Value febcode::Vec2Add(const Value& a, const Value& b)
{
	auto lv = copyStruct(getStruct(a));
	auto& va = lv->fields;
	const auto& vb = getStruct(b).fields;

	va[0] = getDouble(va[0]) + getDouble(vb[0]);
	va[1] = getDouble(va[1]) + getDouble(vb[1]);
	return lv;
}

Value febcode::Vec2Sub(const Value& a, const Value& b)
{
	auto lv = copyStruct(getStruct(a));
	auto& va = lv->fields;
	const auto& vb = getStruct(b).fields;

	va[0] = getDouble(va[0]) - getDouble(vb[0]);
	va[1] = getDouble(va[1]) - getDouble(vb[1]);
	return lv;
}

Value febcode::Vec2Mul(const Value& a, const Value& b)
{
	const auto& va = getStruct(a).fields;
	const auto& vb = getStruct(b).fields;

	double l[2] = { getDouble(va[0]), getDouble(va[1]) };
	double r[2] = { getDouble(vb[0]), getDouble(vb[1]) };

	return l[0] * r[0] + l[1] * r[1];
}

// vec2 = vec2 * double
Value febcode::Vec2ScaleRight(const Value& a, const Value& b)
{
	auto ret = copyStruct(getStruct(a));
	auto& va = ret->fields;

	double vb = getDouble(b);

	va[0] = getDouble(va[0]) * vb;
	va[1] = getDouble(va[1]) * vb;

	return ret;
}

// vec2 = double * vec2
Value febcode::Vec2ScaleLeft(const Value& a, const Value& b)
{
	auto ret = copyStruct(getStruct(b));
	auto& vb = ret->fields;

	double va = getDouble(a);

	vb[0] = getDouble(vb[0]) * va;
	vb[1] = getDouble(vb[1]) * va;

	return ret;
}

Value febcode::NormalizeVec2(const std::vector<Value>& args)
{
	auto& a = getStruct(args[0]);

	double x = getDouble(a.fields[0]);
	double y = getDouble(a.fields[1]);

	double D = std::sqrt(x * x + y * y);
	if (D != 0) {
		x /= D;
		y /= D;
	}

	auto ret = copyStruct(getStruct(args[0]));
	ret->fields[0] = x;
	ret->fields[1] = y;

	return ret;
}
