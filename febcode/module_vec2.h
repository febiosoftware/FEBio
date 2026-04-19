#pragma once
#include "module.h"

namespace febcode
{
	// double = dot(vec2)
	Value DotVec2(FuncArgs args);

	// double = length(vec2)
	Value LengthVec2(FuncArgs args);

	// normalized = Normalize(vec2)
	Value NormalizeVec2(FuncArgs);

	// outer product = outer(vec2, vec2) -> mat2
	Value OuterVec2(FuncArgs args);

	class Vec2Module : public Module
	{
	public:
		void Register(Program& prg) override
		{
			Type bool_t = prg.types.Bool();
			Type doubleType = prg.types.Double();
			Type vec2Type = prg.types.Vec2();
			Type mat2Type = prg.types.Mat2();

			// binary operators                            LHS       RHS       Result
			prg.binaryOps[BinaryOp::Plus    ].push_back({ vec2Type, vec2Type, vec2Type });
			prg.binaryOps[BinaryOp::Minus   ].push_back({ vec2Type, vec2Type, vec2Type });
			prg.binaryOps[BinaryOp::Multiply].push_back({ vec2Type, doubleType, vec2Type });
			prg.binaryOps[BinaryOp::Multiply].push_back({ doubleType, vec2Type, vec2Type });
			prg.binaryOps[BinaryOp::Multiply].push_back({ vec2Type, vec2Type, doubleType }); // dot product
			prg.binaryOps[BinaryOp::Divide  ].push_back({ vec2Type, doubleType, vec2Type });

			prg.binaryOps[BinaryOp::EqualEqual].push_back({ vec2Type, vec2Type, bool_t });
			prg.binaryOps[BinaryOp::NotEqual  ].push_back({ vec2Type, vec2Type, bool_t });

			prg.registerNative("dot"      , doubleType, { vec2Type, vec2Type }, DotVec2);
			prg.registerNative("length"   , doubleType, { vec2Type }, LengthVec2);
			prg.registerNative("normalize", vec2Type, { vec2Type }, NormalizeVec2);
			prg.registerNative("outer"    , mat2Type, { vec2Type, vec2Type }, OuterVec2);
		}
	};
}
