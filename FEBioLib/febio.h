#pragma once
#include <FECore/fecore_export.h>
#include <cstddef>

//-----------------------------------------------------------------------------
// Defines the FEBio namespace
namespace febio 
{
	// Initialize all the FEBio modules
	FECOREDLL_EXPORT void InitLibrary();

	// read the configuration file
	FECOREDLL_EXPORT bool Configure(const char* szfile);

	// load a plugin
	FECOREDLL_EXPORT bool ImportPlugin(const char* szfile);

	// call this to clean up all FEBio data
	FECOREDLL_EXPORT void FinishLibrary();

	// helper function for retrieving the executable's path
	FECOREDLL_EXPORT int get_app_path(char *pname, size_t pathsize);

	// print hello message
	FECOREDLL_EXPORT int Hello();
}
