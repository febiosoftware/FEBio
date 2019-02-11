#pragma once
#include <cstddef>
#include "febiolib_api.h"
#include "FEBioModel.h"
#include <ostream>

class CompactMatrix;

//-----------------------------------------------------------------------------
// Defines the FEBio namespace
namespace febio 
{
	// get the kernel
	FEBIOLIB_API FECoreKernel* GetFECoreKernel();

	// Initialize all the FEBio modules
	FEBIOLIB_API void InitLibrary();

	// read the configuration file
	FEBIOLIB_API bool Configure(const char* szfile);

	// load a plugin
	FEBIOLIB_API bool ImportPlugin(const char* szfile);

	// load all the plugins in a folder
	FEBIOLIB_API void ImportPluginFolder(const char* szfolder);

	// call this to clean up all FEBio data
	FEBIOLIB_API void FinishLibrary();

	// helper function for retrieving the executable's path
	FEBIOLIB_API int get_app_path(char *pname, size_t pathsize);

	// print hello message
	FEBIOLIB_API int Hello();

	// set the number of OMP threads
	FEBIOLIB_API void SetOMPThreads(int n);

	// run an FEBioModel
	FEBIOLIB_API bool SolveModel(FEBioModel& fem, const char* sztask = nullptr, const char* szctrl = nullptr);

	// write a matrix to file
	FEBIOLIB_API bool write_hb(CompactMatrix& K, const char* szfile);

	// print matrix sparsity pattern to svn file
	FEBIOLIB_API void print_svg(CompactMatrix* m, std::ostream &out, int i0 = 0, int j0 = 0, int i1 = -1, int j1 = -1);

	// write a vector to file
	FEBIOLIB_API bool write_vector(const vector<double>& a, const char* szfile);
}
