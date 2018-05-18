#pragma once
#include "fecore_api.h"

//-----------------------------------------------------------------------------
//! The FECore namespace encapsulates all classes that belong to the FECore library
namespace FECore
{
	// retrieve version numbers
	FECORE_API void get_version(int& version, int& subversion);

	// retrieve version number string
	FECORE_API const char* get_version_string();

	// initialize the module
	FECORE_API void InitModule();

} // namespace FECore
