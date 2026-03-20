#pragma once
#include "module.h"

namespace febcode {

	class MathModule : public Module
	{
	public:
		void Register(Program& prg) override
		{
			// Math functions
			prg.registerNative("abs"  , std::abs);
			prg.registerNative("acos" , std::acos);
			prg.registerNative("acosh", std::acosh);
			prg.registerNative("asin" , std::asin);
			prg.registerNative("asinh", std::asinh);
			prg.registerNative("atan" , std::atan);
			prg.registerNative("atanh", std::atanh);
			prg.registerNative("cos"  , std::cos);
			prg.registerNative("cosh" , std::cosh);
			prg.registerNative("exp"  , std::exp);
			prg.registerNative("log"  , std::log);
			prg.registerNative("log10", std::log10);
			prg.registerNative("sin"  , std::sin);
			prg.registerNative("sinh" , std::sinh);
			prg.registerNative("sqrt" , std::sqrt);
			prg.registerNative("tan"  , std::tan);
			prg.registerNative("tanh" , std::tanh);
		}
	};
}