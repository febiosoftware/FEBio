#pragma once
#include "fecore_export.h"

//-----------------------------------------------------------------------------
//! The FECore namespace encapsulates all classes that belong to the FECore library
namespace FECore
{
	// retrieve version numbers
	FECOREDLL_EXPORT void get_version(int& version, int& subversion);

	// retrieve version number string
	FECOREDLL_EXPORT const char* get_version_string();

	// initialize the module
	FECOREDLL_EXPORT void InitModule();

} // namespace FECore
