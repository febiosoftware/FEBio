#pragma once
#include <cstddef>
#include <FECore/FEBox.h>
#include <FECore/quatd.h>
#include "febiolib_api.h"

//-----------------------------------------------------------------------------
// Defines the FEBio namespace
namespace febio 
{
	// Initialize all the FEBio modules
	FEBIOLIB_API void InitLibrary();

	// read the configuration file
	FEBIOLIB_API bool Configure(const char* szfile);

	// load a plugin
	FEBIOLIB_API bool ImportPlugin(const char* szfile);

	// call this to clean up all FEBio data
	FEBIOLIB_API void FinishLibrary();

	// helper function for retrieving the executable's path
	FEBIOLIB_API int get_app_path(char *pname, size_t pathsize);

	// print hello message
	FEBIOLIB_API int Hello();

	// set the number of OMP threads
	FEBIOLIB_API void SetOMPThreads(int n);
}

//-----------------------------------------------------------------------------
// These are some dummy classes to force export of the base class
// TODO: Probably should find a better solution. Maybe create a DLL for the FECore library
class _dummy_FEBox : public FEBox {};
class _dummy_quatd : public quatd {};
