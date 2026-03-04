#pragma once
#include "module.h"

namespace febcode
{
	// vec2 = vec2 + vec2
	Value Vec2Add(const Value& a, const Value& b);

	// vec2 = vec2 - vec2
	Value Vec2Sub(const Value& a, const Value& b);

	// double = vec2 * vec2
	Value Vec2Mul(const Value& a, const Value& b);

	// vec2 = vec2 * double
	Value Vec2ScaleRight(const Value& a, const Value& b);

	// vec2 = double * vec2
	Value Vec2ScaleLeft(const Value& a, const Value& b);

	// normalized = Normalize(vec2)
	Value NormalizeVec2(const std::vector<Value>& args);

	class Vec2Module : public Module
	{
	public:
		void Register(Program& prg) override
		{
			// 1. Register type
			Type vec2Type = prg.types.defineStructType("vec2", {
				{TypeKind::Double, "x"},
				{TypeKind::Double, "y"}
				});

			// 2. Register operators
			Type doubleType = prg.types.getDoubleType();
			prg.RegisterBinaryOperator(BinaryOp::Plus    , vec2Type  , vec2Type  , vec2Type  , Vec2Add);
			prg.RegisterBinaryOperator(BinaryOp::Minus   , vec2Type  , vec2Type  , vec2Type  , Vec2Sub);
			prg.RegisterBinaryOperator(BinaryOp::Multiply, doubleType, vec2Type  , vec2Type  , Vec2Mul);
			prg.RegisterBinaryOperator(BinaryOp::Multiply, vec2Type  , vec2Type  , doubleType, Vec2ScaleRight);
			prg.RegisterBinaryOperator(BinaryOp::Multiply, vec2Type  , doubleType, vec2Type  , Vec2ScaleLeft);

			// 3. Register native functions
			prg.registerNative("normalize", vec2Type, { vec2Type }, NormalizeVec2);

		}
	};
}

