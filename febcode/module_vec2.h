#pragma once
#include "module.h"

namespace febcode
{
	// double = dot(vec2)
	Value DotVec2(FuncArgs args);

	// normalized = Normalize(vec2)
	Value NormalizeVec2(FuncArgs);

	class Vec2Module : public Module
	{
	public:
		void Register(Program& prg) override
		{
			Type vec2Type = prg.types.getVec2Type();

			prg.registerNative("dot"      , vec2Type, { vec2Type }, DotVec2);
			prg.registerNative("normalize", vec2Type, { vec2Type }, NormalizeVec2);
		}
	};
}
