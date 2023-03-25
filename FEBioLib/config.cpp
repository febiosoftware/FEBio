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



#include "stdafx.h"
#include "febio.h"
#include <XML/XMLReader.h>
#include <FEBioXML/xmltool.h>
#include <FECore/FEModel.h>
#include <FECore/FECoreTask.h>
#include <FECore/FEMaterial.h>
#include <NumCore/MatrixTools.h>
#include <FECore/LinearSolver.h>
#include <FEBioTest/FEMaterialTest.h>
#include "plugin.h"
#include <map>
#include <iostream>

#ifdef WIN32
// TODO: This is deprecated and <filesystem> should be used instead when switching to C++17. At that point, also remove this define.
// #include <filesystem>
#define _SILENCE_EXPERIMENTAL_FILESYSTEM_DEPRECATION_WARNING
#include <experimental/filesystem>
#endif

#ifndef WIN32
#include <dlfcn.h>
#endif


namespace febio {

	//-----------------------------------------------------------------------------
	bool parse_tags(XMLTag& tag);
	bool parse_default_linear_solver(XMLTag& tag);
	bool parse_import(XMLTag& tag);
	bool parse_import_folder(XMLTag& tag);
	bool parse_set(XMLTag& tag);
	bool parse_output_negative_jacobians(XMLTag& tag);

	// create a map for the variables (defined with set)
	static std::map<string, string> vars;
	static bool boutput = true;

	//-----------------------------------------------------------------------------
	// configure FEBio
	bool Configure(const char* szfile, FEBioConfig& config)
	{
		vars.clear();

		config.Defaults();
		boutput = (config.m_noutput != 0);

		// open the configuration file
		XMLReader xml;
		if (xml.Open(szfile) == false)
		{
			fprintf(stderr, "FATAL ERROR: Failed reading FEBio configuration file %s.", szfile);
			return false;
		}

		// unload all plugins	
		FEBioPluginManager& pm = *FEBioPluginManager::GetInstance();
		pm.UnloadAllPlugins();

		// loop over all child tags
		try
		{
			// Find the root element
			XMLTag tag;
			if (xml.FindTag("febio_config", tag) == false) return false;
			const char* szversion = tag.AttributeValue("version");
			if (strcmp(szversion, "3.0") == 0)
			{
				if (!tag.isleaf())
				{
					// Read version 1.0
					++tag;
					do
					{
						// parse the tags
						if (tag == "if_debug")
						{
#ifdef _DEBUG
							++tag;
							if (parse_tags(tag) == false) return false;
							++tag;
#else
							tag.skip();
#endif // DEBUG
						}
						else if (tag == "if_release")
						{
#ifndef _DEBUG
							++tag;
							if (parse_tags(tag) == false) return false;
							++tag;
#else
							tag.skip();
#endif // !_DEBUG
						}
						else if (tag == "print_model_params")
						{
							tag.value(config.m_printParams);
						}
						else if (tag == "show_warnings_and_errors")
						{
							tag.value(config.m_bshowErrors);
						}
						else
						{
							if (parse_tags(tag) == false) return false;
						}

						++tag;
					} 
					while (!tag.isend());
				}
			}
			else
			{
				fprintf(stderr, "FATAL ERROR: Invalid version for FEBio configuration file.");
				return false;
			}
		}
		catch (XMLReader::Error& e)
		{
			fprintf(stderr, "FATAL ERROR: %s\n", e.what());
			return false;
		}
		catch (...)
		{
			fprintf(stderr, "FATAL ERROR: unrecoverable error (line %d)", xml.GetCurrentLine());
			return false;
		}

		xml.Close();

		return true;
	}

	//-----------------------------------------------------------------------------
	bool parse_tags(XMLTag& tag)
	{
		if (tag == "set")
		{
			if (parse_set(tag) == false) return false;
		}
		else if (tag == "default_linear_solver")
		{
			if (parse_default_linear_solver(tag) == false) return false;
		}
		else if (tag == "import")
		{
			if (parse_import(tag) == false) return false;
		}
		else if (tag == "import_folder")
		{
			if (parse_import_folder(tag) == false) return false;
		}
		else if (tag == "output_negative_jacobians")
		{
			if (parse_output_negative_jacobians(tag) == false) return false;
		}
		else throw XMLReader::InvalidTag(tag);

		return true;
	}

	//-----------------------------------------------------------------------------
	bool parse_set(XMLTag& tag)
	{
		char szname[256] = { 0 };
		strcpy(szname, tag.AttributeValue("name"));
		string key(szname);
		string val(tag.szvalue());
		vars[key] = val;
		return true;
	}

	//-----------------------------------------------------------------------------
	bool parse_output_negative_jacobians(XMLTag& tag)
	{
		int n;
		tag.value(n);
		NegativeJacobian::m_boutput = (n != 0);
		return true;
	}

	//-----------------------------------------------------------------------------
	bool parse_default_linear_solver(XMLTag& tag)
	{
		const char* szt = tag.AttributeValue("type");

		// read the solver parameters
		FEClassDescriptor* cd = fexml::readParameterList(tag);
		if (cd == nullptr)
		{
			delete cd;
			return false;
		}
		else
		{
			// set this as the default solver
			FECoreKernel& fecore = FECoreKernel::GetInstance();
			fecore.SetDefaultSolver(cd);
			if (boutput) fprintf(stderr, "Default linear solver: %s\n", fecore.GetLinearSolverType());
		}

		return true;
	}
	
	//-----------------------------------------------------------------------------
	bool process_aliases(char* szout, const char* szbuf)
	{
		bool bok = true;
		char sztmp[64] = { 0 };
		const char* ch = szbuf;
		char* s = szout;
		while (*ch)
		{
			if (*ch == '$')
			{
				++ch;
				if (*ch++ == '(')
				{
					const char* ch2 = strchr(ch, ')');
					if (ch2)
					{
						int l = (int)(ch2 - ch);
						strncpy(sztmp, ch, l);
						sztmp[l] = 0;
						string key(sztmp);
						ch = ch2 + 1;
						std::map<string, string>::iterator it = vars.find(key);
						if (it != vars.end())
						{
							string v = it->second;
							const char* sz = v.c_str();
							while (*sz) *s++ = *sz++;
						}
						else { bok = false; break; }
					}
					else { bok = false; break; }

				}
				else { bok = false; break; }
			}
			else *s++ = *ch++;
		}
		return bok;
	}

	//-----------------------------------------------------------------------------
	bool parse_import(XMLTag& tag)
	{
		// get the file name
		const char* szfile = tag.szvalue();

		// process any aliases
		char szbuf[1024] = { 0 };
		bool bok = process_aliases(szbuf, szfile);

		// load the plugin
		if (bok) febio::ImportPlugin(szbuf);

		return bok;
	}

	//-----------------------------------------------------------------------------
	bool parse_import_folder(XMLTag& tag)
	{
		// get the folder name
		const char* szfolder = tag.szvalue();

		// process any aliases
		char szbuf[1024] = { 0 };
		bool bok = process_aliases(szbuf, szfolder);

		// load the plugin
		if (bok) febio::ImportPluginFolder(szbuf);

		return bok;
	}

	//-----------------------------------------------------------------------------
	const char* GetFileTitle(const char* szfile)
	{
		const char* ch = strrchr(szfile, '\\');
		if (ch == 0) {
			ch = strrchr(szfile, '/');
			if (ch == 0) ch = szfile; else ch++;
		}
		else ch++;
		return ch;
	}

	//-----------------------------------------------------------------------------
	bool ImportPlugin(const char* szfile)
	{
		const char* sztitle = GetFileTitle(szfile);
		FEBioPluginManager* pPM = FEBioPluginManager::GetInstance();
		PLUGIN_INFO info;
		int nerr = pPM->LoadPlugin(szfile, info);
		switch (nerr)
		{
		case 0: if (boutput) fprintf(stderr, "Success loading plugin %s (version %d.%d.%d)\n", sztitle, info.major, info.minor, info.patch); return true; break;
		case 1:
			fprintf(stderr, "Failed loading plugin %s\n Reason: Failed to load the file.\n\n", szfile);
#ifndef WIN32
			fprintf(stderr, "dlopen failed: %s\n\n", dlerror());
#endif
			break;
		case 2: fprintf(stderr, "Failed loading plugin %s\n Reason: Required plugin function PluginNumClasses not found.\n\n", szfile); break;
		case 3: fprintf(stderr, "Failed loading plugin %s\n Reason: Required plugin function PluginGetFactory not found.\n\n", szfile); break;
		case 4: fprintf(stderr, "Failed loading plugin %s\n Reason: Invalid number of classes returned by PluginNumClasses.\n\n", szfile); break;
		case 5: fprintf(stderr, "Failed loading plugin %s\n Reason: Required plugin function GetSDKVersion not found.\n\n", szfile); break;
		case 6: fprintf(stderr, "Failed loading plugin %s\n Reason: Invalid SDK version.\n\n", szfile); break;
		case 7: fprintf(stderr, "Failed loading plugin %s\n Reason: Plugin is already loaded.\n\n", szfile); break;
		default:
			fprintf(stderr, "Failed loading plugin %s\n Reason: unspecified.\n\n", szfile); break;
		}

		return false;
	}

	//-----------------------------------------------------------------------------

#ifdef WIN32
	namespace fs = std::experimental::filesystem;

	void ImportPluginFolder(const char* szfolder)
	{
		std::string path = szfolder;

		// get the default (system-dependant) extension
		std::wstring defExt = L".dll";
		//	std::wstring defExt = L".dylib";
		//	std::wstring defExt = L".so";

		size_t extLength = defExt.length();

		// loop over all the items in a directory
		for (auto & p : fs::directory_iterator(path))
		{
			// only get regular files with the default extension 
			if (p.status().type() == fs::file_type::regular)
			{
				std::wstring fileName = p.path();
				size_t l = fileName.length();
				if (l > extLength) {
					std::wstring ext = fileName.substr(l - extLength, extLength);
					if (ext == defExt)
					{
						// we can only deal with strings for now, so convert
						std::string s = p.path().string<char>();

						// try to load the plugin
						ImportPlugin(s.c_str());
					}
				}
			}
		}
	}
#else
void ImportPluginFolder(const char* szfolder)
{
}
#endif

//-----------------------------------------------------------------------------
FEBIOLIB_API const char* GetPluginName(int allocId)
{
	FEBioPluginManager& pm = *FEBioPluginManager::GetInstance();

	for (int i = 0; i < pm.Plugins(); ++i)
	{
		const FEBioPlugin& pi = pm.GetPlugin(i);
		if (pi.GetAllocatorID() == allocId)
		{
			return pi.GetName();
		}
	}
	return nullptr;
}

//-----------------------------------------------------------------------------
// run an FEBioModel
FEBIOLIB_API bool SolveModel(FEBioModel& fem, const char* sztask, const char* szctrl)
{
	// Make sure we have a task
	if (sztask == nullptr) sztask = "solve";

	// find a task
	FECoreTask* ptask = fecore_new<FECoreTask>(sztask, &fem);
	if (ptask == 0)
	{
		fprintf(stderr, "Don't know how to do task: %s\n", sztask);
		return false;
	}

	// initialize the task
	if (ptask->Init(szctrl) == false)
	{
		fprintf(stderr, "Failed initializing the task: %s\n", sztask);
		return false;
	}

	// run the task
	bool bret = true;
	try {
		bret = ptask->Run();
	}
	catch (std::exception e)
	{
		fprintf(stderr, "\nException detected: %s\n\n", e.what());
		bret = false;
	}

	return bret;
}

//-----------------------------------------------------------------------------
// run an FEBioModel
FEBIOLIB_API int RunModel(FEBioModel& fem, CMDOPTIONS* ops)
{
	// set options that were passed on the command line
	if (ops)
	{
		fem.SetDebugLevel(ops->ndebug);
		fem.SetDumpLevel(ops->dumpLevel);

		// set the output filenames
		fem.SetLogFilename(ops->szlog);
		fem.SetPlotFilename(ops->szplt);
		fem.SetDumpFilename(ops->szdmp);
	}

	// read the input file if specified
	int nret = 0;
	if (ops && ops->szfile[0])
	{
		// read the input file
		if (fem.Input(ops->szfile) == false) nret = 1;
	}

	// solve the model with the task and control file
	if (nret == 0)
	{
		const char* sztask = (ops && ops->sztask[0] ? ops->sztask : nullptr);
		const char* szctrl = (ops && ops->szctrl[0] ? ops->szctrl : nullptr);
		bool b = febio::SolveModel(fem, sztask, szctrl);
		nret = (b ? 0 : 1);
	}

	return nret;
}

// write a matrix to file
bool write_hb(CompactMatrix& K, const char* szfile, int mode)
{
	return NumCore::write_hb(K, szfile, mode);
}

// print matrix sparsity pattern to svn file
void print_svg(CompactMatrix* m, std::ostream &out, int i0, int j0, int i1, int j1)
{
	NumCore::print_svg(m, out, i0, j0, i1, j1);
}

// write a vector to file
bool write_vector(const vector<double>& a, const char* szfile, int mode)
{
	return NumCore::write_vector(a, szfile, mode);
}

bool RunMaterialTest(FEMaterial* mat, double simtime, int steps, double strain, const char* sztest, std::vector<pair<double, double> >& out)
{
	FEModel fem;

	FEMaterial* matcopy = dynamic_cast<FEMaterial*>(CopyFEBioClass(mat, &fem));
	if (matcopy == nullptr) return false;

	fem.AddMaterial(matcopy);

	FECoreKernel& febio = FECoreKernel::GetInstance();

	FEMaterialTest diag(&fem);
	diag.SetOutputFileName(nullptr);

	FEDiagnosticScenario* s = diag.CreateScenario(sztest);
	s->GetParameterList();
	s->SetParameter<double>("strain", strain);

	FEAnalysis* step = fem.GetStep(0);
	step->m_ntime = steps;
	step->m_dt0 = simtime / steps;
	fem.SetCurrentStepIndex(0);

	if (diag.Init() == false) return false;

	if (fem.Init() == false) return false;

	bool b = diag.Run();

	if (b)
	{
		out = diag.GetOutputData();
	}

	return b;
}

} // namespace febio
