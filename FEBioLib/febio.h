/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/



#pragma once
#include <cstddef>
#include "febiolib_api.h"
#include "FEBioModel.h"
#include "FEBioConfig.h"
#include "cmdoptions.h"
#include <ostream>

class CompactMatrix;
class LogStream;

//-----------------------------------------------------------------------------
// Defines the FEBio namespace
namespace febio 
{
	// get the kernel
	FEBIOLIB_API FECoreKernel* GetFECoreKernel();

	// Initialize all the FEBio modules
	FEBIOLIB_API void InitLibrary();

	// read the configuration file
	FEBIOLIB_API bool Configure(const char* szfile, FEBioConfig& config);

	// load a plugin
	FEBIOLIB_API bool ImportPlugin(const char* szfile);

	// load all the plugins in a folder
	FEBIOLIB_API void ImportPluginFolder(const char* szfolder);

	// get the name of the plugin from its allocator Id
	FEBIOLIB_API const char* GetPluginName(int allocId);

	// call this to clean up all FEBio data
	FEBIOLIB_API void FinishLibrary();

	// helper function for retrieving the executable's path
	FEBIOLIB_API int get_app_path(char *pname, size_t pathsize);

	// print hello message
	FEBIOLIB_API int Hello(LogStream& log);

	// set the number of OMP threads
	FEBIOLIB_API void SetOMPThreads(int n);

	// run an FEBioModel
	FEBIOLIB_API bool SolveModel(FEBioModel& fem, const char* sztask = nullptr, const char* szctrl = nullptr);

	// run an FEBioModel
	FEBIOLIB_API int RunModel(FEBioModel& fem, CMDOPTIONS* ops);

	// write a matrix to file
	FEBIOLIB_API bool write_hb(CompactMatrix& K, const char* szfile, int mode = 0);

	// print matrix sparsity pattern to svn file
	FEBIOLIB_API void print_svg(CompactMatrix* m, std::ostream &out, int i0 = 0, int j0 = 0, int i1 = -1, int j1 = -1);

	// write a vector to file
	FEBIOLIB_API bool write_vector(const vector<double>& a, const char* szfile, int mode = 0);

	// run a material test
	FEBIOLIB_API bool RunMaterialTest(FEMaterial* mat, double simtime, int steps, double strain, const char* sztest, std::vector<pair<double, double> >& out);
}
