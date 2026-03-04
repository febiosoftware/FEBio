#include "module_vec3.h"

using namespace febcode;

Value febcode::Vec3Add(const Value& a, const Value& b)
{
	auto lv = copyStruct(getStruct(a));
	auto& va = lv->fields;
	const auto& vb = getStruct(b).fields;

	va[0] = getDouble(va[0]) + getDouble(vb[0]);
	va[1] = getDouble(va[1]) + getDouble(vb[1]);
	va[2] = getDouble(va[2]) + getDouble(vb[2]);
	return lv;
}

Value febcode::Vec3Sub(const Value& a, const Value& b)
{
	auto lv = copyStruct(getStruct(a));
	auto& va = lv->fields;
	const auto& vb = getStruct(b).fields;

	va[0] = getDouble(va[0]) - getDouble(vb[0]);
	va[1] = getDouble(va[1]) - getDouble(vb[1]);
	va[2] = getDouble(va[2]) - getDouble(vb[2]);
	return lv;
}

Value febcode::Vec3Mul(const Value& a, const Value& b)
{
	const auto& va = getStruct(a).fields;
	const auto& vb = getStruct(b).fields;

	double l[3] = { getDouble(va[0]), getDouble(va[1]), getDouble(va[2]) };
	double r[3] = { getDouble(vb[0]), getDouble(vb[1]), getDouble(vb[2]) };

	return l[0] * r[0] + l[1] * r[1] + l[2] * r[2];
}

// vec3 = vec3 * double
Value febcode::Vec3ScaleRight(const Value& a, const Value& b)
{
	auto ret = copyStruct(getStruct(a));
	auto& va = ret->fields;

	double vb = getDouble(b);

	va[0] = getDouble(va[0]) * vb;
	va[1] = getDouble(va[1]) * vb;
	va[2] = getDouble(va[2]) * vb;

	return ret;
}

// vec3 = double * vec3
Value febcode::Vec3ScaleLeft(const Value& a, const Value& b)
{
	auto ret = copyStruct(getStruct(b));
	auto& vb = ret->fields;

	double va = getDouble(a);

	vb[0] = getDouble(vb[0]) * va;
	vb[1] = getDouble(vb[1]) * va;
	vb[2] = getDouble(vb[2]) * va;

	return ret;
}

Value febcode::CrossProduct(const std::vector<Value>& args)
{
	auto& a = getStruct(args[0]);
	auto& b = getStruct(args[1]);

	double x = getDouble(a.fields[1]) * getDouble(b.fields[2]) - getDouble(a.fields[2]) * getDouble(b.fields[1]);
	double y = getDouble(a.fields[2]) * getDouble(b.fields[0]) - getDouble(a.fields[0]) * getDouble(b.fields[2]);
	double z = getDouble(a.fields[0]) * getDouble(b.fields[1]) - getDouble(a.fields[1]) * getDouble(b.fields[0]);

	auto ret = copyStruct(getStruct(args[0]));
	ret->fields[0] = x;
	ret->fields[1] = y;
	ret->fields[2] = z;

	return ret;
}

Value febcode::NormalizeVec3(const std::vector<Value>& args)
{
	auto& a = getStruct(args[0]);

	double x = getDouble(a.fields[0]);
	double y = getDouble(a.fields[1]);
	double z = getDouble(a.fields[2]);

	double D = std::sqrt(x * x + y * y + z * z);
	if (D != 0) {
		x /= D;
		y /= D;
		z /= D;
	}

	auto ret = copyStruct(getStruct(args[0]));
	ret->fields[0] = x;
	ret->fields[1] = y;
	ret->fields[2] = z;

	return ret;
}
