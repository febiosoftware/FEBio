#pragma once

//-----------------------------------------------------------------------------
// Defines the FEBio namespace
namespace febio 
{
	// Initialize all the FEBio modules
	void InitLibrary();

	// read the configuration file
	bool Configure(const char* szfile);

	// load a plugin
	bool ImportPlugin(const char* szfile);

	// call this to clean up all FEBio data
	void FinishLibrary();
}
