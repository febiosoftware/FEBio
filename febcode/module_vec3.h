#pragma once
#include "module.h"

namespace febcode
{
	// vec3 = vec3 + vec3
	Value Vec3Add(const Value& a, const Value& b);

	// vec3 = vec3 - vec3
	Value Vec3Sub(const Value& a, const Value& b);

	// double = vec3 * vec3
	Value Vec3Mul(const Value& a, const Value& b);

	// vec3 = vec3 * double
	Value Vec3ScaleRight(const Value& a, const Value& b);

	// vec3 = double * vec3
	Value Vec3ScaleLeft(const Value& a, const Value& b);

	// vec3 = cross(vec3, vec3)
	Value CrossProduct(const std::vector<Value>& args);

	// normalized = Normalize(vec3)
	Value NormalizeVec3(const std::vector<Value>& args);

	class Vec3Module : public Module
	{
	public:
		void Register(Program& prg) override
		{
			// 1. Register type
			Type vec3Type = prg.types.defineStructType("vec3", {
				{TypeKind::Double, "x"},
				{TypeKind::Double, "y"},
				{TypeKind::Double, "z"}
				});

			// 2. Register operators
			Type doubleType = prg.types.getDoubleType();
			prg.RegisterBinaryOperator(BinaryOp::Plus    , vec3Type, vec3Type, vec3Type, Vec3Add);
			prg.RegisterBinaryOperator(BinaryOp::Minus   , vec3Type, vec3Type, vec3Type, Vec3Sub);
			prg.RegisterBinaryOperator(BinaryOp::Multiply, doubleType, vec3Type, vec3Type, Vec3Mul);
			prg.RegisterBinaryOperator(BinaryOp::Multiply, vec3Type, vec3Type, doubleType, Vec3ScaleRight);
			prg.RegisterBinaryOperator(BinaryOp::Multiply, vec3Type, doubleType, vec3Type, Vec3ScaleLeft);

			// 3. Register native functions
			prg.registerNative("cross", vec3Type, { vec3Type, vec3Type }, CrossProduct);
			prg.registerNative("normalize", vec3Type, { vec3Type }, NormalizeVec3);
		}
	};
}
