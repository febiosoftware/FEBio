#pragma once
#include "module.h"

namespace febcode
{
	Value TransposeMat2(FuncArgs args);

	Value InvertMat2(FuncArgs args);

	class Mat2Module : public Module
	{
	public:
		void Register(Program& prg) override
		{
			Type bool_t = prg.types.Bool();
			Type double_t = prg.types.Double();
			Type vec2_t = prg.types.Vec2();
			Type mat2_t = prg.types.Mat2();

			// binary operators                            LHS     RHS    Result
			prg.binaryOps[BinaryOp::Plus    ].push_back({ mat2_t, mat2_t, mat2_t });
			prg.binaryOps[BinaryOp::Minus   ].push_back({ mat2_t, mat2_t, mat2_t });
			prg.binaryOps[BinaryOp::Multiply].push_back({ mat2_t, mat2_t, mat2_t });
			prg.binaryOps[BinaryOp::Multiply].push_back({ mat2_t, vec2_t, vec2_t });
			prg.binaryOps[BinaryOp::Multiply].push_back({ mat2_t, double_t, mat2_t });
			prg.binaryOps[BinaryOp::Multiply].push_back({ double_t, mat2_t, mat2_t });
			prg.binaryOps[BinaryOp::Divide  ].push_back({ mat2_t, double_t, mat2_t });

			prg.binaryOps[BinaryOp::EqualEqual].push_back({ mat2_t, mat2_t, bool_t });
			prg.binaryOps[BinaryOp::NotEqual  ].push_back({ mat2_t, mat2_t, bool_t });

			prg.registerNative("transpose", mat2_t, { mat2_t }, TransposeMat2);
			prg.registerNative("inverse", mat2_t, { mat2_t }, InvertMat2);
		}
	};
}
