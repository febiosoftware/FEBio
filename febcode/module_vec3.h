#pragma once
#include "module.h"

namespace febcode
{
	// vec3 = dot(vec3, vec3)
	Value DotVec3(const std::vector<Value>& args);

	// vec3 = cross(vec3, vec3)
	Value CrossVec3(const std::vector<Value>& args);

	// normalized = Normalize(vec3)
	Value NormalizeVec3(const std::vector<Value>& args);

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
