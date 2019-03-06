#include "stdafx.h"
#include "febio.h"
#include <FEBioXML/XMLReader.h>
#include <FEBioXML/xmltool.h>
#include <FECore/FEModel.h>
#include <FECore/FECoreTask.h>
#include <NumCore/MatrixTools.h>
#include "plugin.h"
#include <map>
#include <iostream>

#ifdef WIN32
#include <filesystem>
#endif

#ifndef WIN32
#include <dlfcn.h>
#endif

#ifdef WIN32
extern "C" void __cdecl omp_set_num_threads(int);
#else
extern "C" void omp_set_num_threads(int);
#endif

namespace febio {

	//-----------------------------------------------------------------------------
	bool parse_tags(XMLTag& tag);
	bool parse_linear_solver_props(XMLTag& tag);
	bool parse_default_linear_solver(XMLTag& tag);
	bool parse_import(XMLTag& tag);
	bool parse_import_folder(XMLTag& tag);
	bool parse_set(XMLTag& tag);
	bool parse_omp_num_threads(XMLTag& tag);
	bool parse_output_negative_jacobians(XMLTag& tag);

	// create a map for the variables (defined with set)
	static std::map<string, string> vars;

	//-----------------------------------------------------------------------------
	// configure FEBio
	bool Configure(const char* szfile)
	{
		vars.clear();

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

			if (strcmp(tag.m_att[0].m_szatv, "3.0") == 0)
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
							xml.SkipTag(tag);
#endif // DEBUG
						}
						else if (tag == "if_release")
						{
#ifndef _DEBUG
							++tag;
							if (parse_tags(tag) == false) return false;
							++tag;
#else
							xml.SkipTag(tag);
#endif // !_DEBUG
						}
						else
						{
							if (parse_tags(tag) == false) return false;
						}
					} while (!tag.isend());
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
			fprintf(stderr, "FATAL ERROR: %s (line %d)\n", e.GetErrorString(), xml.GetCurrentLine());
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
		else if (tag == "linear_solver_props")
		{
			if (parse_linear_solver_props(tag) == false) return false;
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
		else if (tag == "omp_num_threads")
		{
			if (parse_omp_num_threads(tag) == false) return false;
		}
		else if (tag == "output_negative_jacobians")
		{
			if (parse_output_negative_jacobians(tag) == false) return false;
		}
		else throw XMLReader::InvalidTag(tag);

		++tag;
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
	bool parse_omp_num_threads(XMLTag& tag)
	{
		int n;
		tag.value(n);
		omp_set_num_threads(n);
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
	bool parse_linear_solver_props(XMLTag& tag)
	{
		FECoreKernel& fecore = FECoreKernel::GetInstance();

		const char* szt = tag.AttributeValue("type");
		if (tag.isleaf() == false)
		{
			FECoreFactory* fac = fecore.FindFactoryClass(FELINEARSOLVER_ID, szt);
			if (fac == 0) throw XMLReader::InvalidAttributeValue(tag, "type", szt);

			// read the solver parameters
			FEParameterList& PL = fac->GetParameterList();
			++tag;
			do
			{
				FEParam* p = PL.FindFromName(tag.m_sztag);
				if (p)
				{
					if (fexml::readParameter(tag, PL) == false)
					{
						fprintf(stderr, "Invalid linear solver parameter value\n"); return false;
					}
				}
				else { fprintf(stderr, "Invalid linear solver parameter\n"); return false; }
				++tag;
			} while (!tag.isend());
		}

		return true;
	}

	//-----------------------------------------------------------------------------
	bool parse_default_linear_solver(XMLTag& tag)
	{
		FECoreKernel& fecore = FECoreKernel::GetInstance();

		const char* szt = tag.AttributeValue("type");
		FECoreFactory* fac = fecore.SetDefaultSolver(szt);
		if (fac == nullptr)
		{
			fprintf(stderr, "Invalid linear solver\n");
			return false;
		}
		else fprintf(stdout, "Default linear solver: %s\n", szt);

		if (tag.isleaf() == false) parse_linear_solver_props(tag);

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
		case 0: fprintf(stderr, "Success loading plugin %s (version %d.%d.%d)\n", sztitle, info.major, info.minor, info.patch); return true; break;
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
void SetOMPThreads(int n)
{
	omp_set_num_threads(n);
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
	try {
		bool bret = ptask->Run();
	}
	catch (std::exception e)
	{
		fprintf(stderr, "\nException detected: %s\n\n", e.what());
		return false;
	}

	return true;
}

// write a matrix to file
bool write_hb(CompactMatrix& K, const char* szfile)
{
	return NumCore::write_hb(K, szfile);
}

// print matrix sparsity pattern to svn file
void print_svg(CompactMatrix* m, std::ostream &out, int i0, int j0, int i1, int j1)
{
	NumCore::print_svg(m, out, i0, j0, i1, j1);
}

// write a vector to file
bool write_vector(const vector<double>& a, const char* szfile)
{
	return NumCore::write_vector(a, szfile);
}

} // namespace febio
