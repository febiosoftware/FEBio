#pragma once
#include "module.h"

namespace febcode
{
	Value TransposeMat3(FuncArgs args);

	Value InvertMat3(FuncArgs args);

	// outer product = outer(vec3, vec3) -> mat3
	Value OuterVec3(FuncArgs args);

	class Mat3Module : public Module
	{
	public:
		void Register(Program& prg) override
		{
			Type bool_t = prg.types.Bool();
			Type double_t = prg.types.Double();
			Type vec3_t = prg.types.Vec3();
			Type mat3_t = prg.types.Mat3();

			// binary operators                            LHS     RHS    Result
			prg.binaryOps[BinaryOp::Plus    ].push_back({ mat3_t, mat3_t, mat3_t });
			prg.binaryOps[BinaryOp::Minus   ].push_back({ mat3_t, mat3_t, mat3_t });
			prg.binaryOps[BinaryOp::Multiply].push_back({ mat3_t, mat3_t, mat3_t });
			prg.binaryOps[BinaryOp::Multiply].push_back({ mat3_t, vec3_t, vec3_t });
			prg.binaryOps[BinaryOp::Multiply].push_back({ mat3_t, double_t, mat3_t });
			prg.binaryOps[BinaryOp::Multiply].push_back({ double_t, mat3_t, mat3_t });
			prg.binaryOps[BinaryOp::Divide  ].push_back({ mat3_t, double_t, mat3_t });

			prg.binaryOps[BinaryOp::EqualEqual].push_back({ mat3_t, mat3_t, bool_t });
			prg.binaryOps[BinaryOp::NotEqual  ].push_back({ mat3_t, mat3_t, bool_t });

			prg.registerNative("transpose", mat3_t, { mat3_t }, TransposeMat3);
			prg.registerNative("inverse", mat3_t, { mat3_t }, InvertMat3);
			prg.registerNative("outer", mat3_t, { vec3_t, vec3_t }, OuterVec3);
		}
	};
}

