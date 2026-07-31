#pragma once
#include "program.h"

namespace febcode 
{
	class Module
	{
	public:
		virtual void Register(Program& prg) = 0;
		virtual ~Module() = default;
	};
}
