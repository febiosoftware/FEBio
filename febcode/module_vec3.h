#pragma once
#include "module.h"

namespace febcode
{
	// vec3 = dot(vec3, vec3)
	Value DotVec3(FuncArgs args);

	// vec3 = cross(vec3, vec3)
	Value CrossVec3(FuncArgs args);

	// normalized = Normalize(vec3)
	Value NormalizeVec3(FuncArgs args);

	class Vec3Module : public Module
	{
	public:
		void Register(Program& prg) override
		{
			Type bool_t = prg.types.Bool();
			Type doubleType = prg.types.Double();
			Type vec3Type = prg.types.Vec3();
			Type mat3Type = prg.types.Mat3();

			// binary operators                            LHS       RHS       Result
			prg.binaryOps[BinaryOp::Plus    ].push_back({ vec3Type, vec3Type, vec3Type });
			prg.binaryOps[BinaryOp::Minus   ].push_back({ vec3Type, vec3Type, vec3Type });
			prg.binaryOps[BinaryOp::Multiply].push_back({ vec3Type, doubleType, vec3Type });
			prg.binaryOps[BinaryOp::Multiply].push_back({ doubleType, vec3Type, vec3Type });
			prg.binaryOps[BinaryOp::Multiply].push_back({ vec3Type, vec3Type, doubleType }); // dot product

			prg.binaryOps[BinaryOp::EqualEqual].push_back({ vec3Type, vec3Type, bool_t });
			prg.binaryOps[BinaryOp::NotEqual  ].push_back({ vec3Type, vec3Type, bool_t });

			// Register native functions
			prg.registerNative("dot"      , vec3Type, { vec3Type, vec3Type }, DotVec3);
			prg.registerNative("cross"    , vec3Type, { vec3Type, vec3Type }, CrossVec3);
			prg.registerNative("normalize", vec3Type, { vec3Type }, NormalizeVec3);
		}
	};
}
