#pragma once
#include <FECore/fecore_export.h>

//-----------------------------------------------------------------------------
// Defines the FEBio namespace
namespace febio 
{
	// Initialize all the FEBio modules
	FECORE_EXPORT void InitLibrary();

	// read the configuration file
	FECORE_EXPORT bool Configure(const char* szfile);

	// load a plugin
	FECORE_EXPORT bool ImportPlugin(const char* szfile);

	// call this to clean up all FEBio data
	FECORE_EXPORT void FinishLibrary();

	// helper function for retrieving the executable's path
	FECORE_EXPORT int get_app_path(char *pname, size_t pathsize);

	// print hello message
	FECORE_EXPORT int Hello();
}
