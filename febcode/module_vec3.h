#pragma once
#include "module.h"

namespace febcode
{
	// vec3 = dot(vec3, vec3)
	Value DotVec3(const Value* args, int argc);

	// vec3 = cross(vec3, vec3)
	Value CrossVec3(const Value* args, int argc);

	// normalized = Normalize(vec3)
	Value NormalizeVec3(const Value* args, int argc);

	class Vec3Module : public Module
	{
	public:
		void Register(Program& prg) override
		{
			Type vec3Type = prg.types.getVec3Type();

			// Register native functions
			prg.registerNative("dot"      , vec3Type, { vec3Type, vec3Type }, DotVec3);
			prg.registerNative("cross"    , vec3Type, { vec3Type, vec3Type }, CrossVec3);
			prg.registerNative("normalize", vec3Type, { vec3Type }, NormalizeVec3);
		}
	};
}
