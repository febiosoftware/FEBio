#include "stdafx.h"
#include "FECore/FECore.h"
#include "NumCore/NumCore.h"
#include "FEBioMech/FEBioMech.h"
#include "FEBioMix/FEBioMix.h"
#include "FEBioHeat/FEBioHeat.h"
#include "FEBioOpt/FEBioOpt.h"
#include "FEBioFluid/FEBioFluid.h"
#include <FEBioXML/XMLReader.h>
#include <FECore/FEModel.h>
#include <FECore/log.h>
#include "febio.h"
#include "plugin.h"
#include <map>

#ifndef WIN32
#include <dlfcn.h>
#endif

namespace febio {

//-----------------------------------------------------------------------------
// import all modules
void InitLibrary()
{
	FECore::InitModule();
	NumCore::InitModule();
	FEBioMech::InitModule();
	FEBioMix::InitModule();
	FEBioHeat::InitModule();
	FEBioOpt::InitModule();
	FEBioFluid::InitModule();
}

//-----------------------------------------------------------------------------
// configure FEBio
bool Configure(const char* szfile)
{
	// create a map for the variables (defined with set)
	std::map<string, string> vars;

	// open the configuration file
	XMLReader xml;
	if (xml.Open(szfile) == false)
	{
		fprintf(stderr, "FATAL ERROR: Failed reading FEBio configuration file %s.", szfile);
		return false;
	}

	// loop over all child tags
	try
	{
		// Find the root element
		XMLTag tag;
		if (xml.FindTag("febio_config", tag) == false) return false;

		if (strcmp(tag.m_att[0].m_szatv, "1.0") == 0)
		{
			if (!tag.isleaf())
			{
				// Read version 1.0
				++tag;
				do
				{
					if (tag == "set")
					{
						char szname[256] = {0};
						strcpy(szname, tag.AttributeValue("name"));
						string key(szname);
						string val(tag.szvalue());
						vars[key] = val;
					}
					else if (tag == "linear_solver")
					{
						const char* szt = tag.AttributeValue("type");
						if      (strcmp(szt, "skyline"           ) == 0) FECoreKernel::SetDefaultSolver(SKYLINE_SOLVER   );
						else if (strcmp(szt, "psldlt"            ) == 0) FECoreKernel::SetDefaultSolver(PSLDLT_SOLVER    );
						else if (strcmp(szt, "superlu"           ) == 0) FECoreKernel::SetDefaultSolver(SUPERLU_SOLVER   );
						else if (strcmp(szt, "superlu_mt"        ) == 0) FECoreKernel::SetDefaultSolver(SUPERLU_MT_SOLVER);
						else if (strcmp(szt, "pardiso"           ) == 0) FECoreKernel::SetDefaultSolver(PARDISO_SOLVER   );
						else if (strcmp(szt, "rcicg"             ) == 0) FECoreKernel::SetDefaultSolver(RCICG_SOLVER     );
						else if (strcmp(szt, "fgmres"            ) == 0) FECoreKernel::SetDefaultSolver(FGMRES_SOLVER    );
						else if (strcmp(szt, "fgmres_ilut"       ) == 0) FECoreKernel::SetDefaultSolver(FGMRES_ILUT_SOLVER);
						else if (strcmp(szt, "fgmres_ilu0"       ) == 0) FECoreKernel::SetDefaultSolver(FGMRES_ILU0_SOLVER);
						else if (strcmp(szt, "wsmp"              ) == 0) FECoreKernel::SetDefaultSolver(WSMP_SOLVER      );
					}
					else if (tag == "import")
					{
						const char* szfile = tag.szvalue();

						char szbuf[1024] = {0};
						strcpy(szbuf, szfile);

						bool bok = true;
						char sz[2048] = {0};
						char* ch = szbuf;
						char* s = sz;
						while (*ch)
						{
							if (*ch == '$')
							{
								++ch;
								if (*ch++ == '(')
								{
									char* ch2 = strchr(ch, ')');
									if (ch2)
									{
										*ch2 = 0;
										string key(ch);
										ch = ch2+1;
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

						if (bok) febio::ImportPlugin(sz);
					}
/*					else if (tag == "import_folder")
					{
						const char* szfile = tag.szvalue();
						if (LoadPluginFolder(szfile) == false) throw XMLReader::InvalidTag(tag);
					}
*/					else if (tag == "omp_num_threads")
					{
						int n;
						tag.value(n);
						  //						omp_set_num_threads(n);
					}
					else if (tag == "output_negative_jacobians")
					{
						int n;
						tag.value(n);
						NegativeJacobian::m_boutput = (n != 0);
					}
					else throw XMLReader::InvalidTag(tag);

					// go to the next tag
					++tag;
				}
				while (!tag.isend());
			}
		}
		else
		{
			felog.printbox("FATAL ERROR", "Invalid version for FEBio configuration file.");
			return false;
		}
	}
	catch (XMLReader::Error& e)
	{
		felog.printf("FATAL ERROR: %s (line %d)\n", e.GetErrorString(), xml.GetCurrentLine());
		return false;
	}
	catch (...)
	{
		felog.printf("FATAL ERROR: unrecoverable error (line %d)", xml.GetCurrentLine());
		return false;
	}

	xml.Close();

	return true;
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
		fprintf(stderr, "Failed loading plugin %s\n Reason: Failed to load the file.\n\n", szfile); break;
#ifndef WIN32
		fprintf(stderr, "dlopen failed: %s\n\n", dlerror());
#endif
		break;
	case 2: fprintf(stderr, "Failed loading plugin %s\n Reason: Required plugin function PluginNumClasses not found.\n\n", szfile); break;
	case 3: fprintf(stderr, "Failed loading plugin %s\n Reason: Required plugin function PluginGetFactory not found.\n\n", szfile); break;
	case 4: fprintf(stderr, "Failed loading plugin %s\n Reason: Invalid number of classes returned by PluginNumClasses.\n\n", szfile); break;
	case 5: fprintf(stderr, "Failed loading plugin %s\n Reason: Required plugin function GetSDKVersion not found.\n\n", szfile); break;
	case 6: fprintf(stderr, "Failed loading plugin %s\n Reason: Invalid SDK version.\n\n", szfile); break;
	default:
		fprintf(stderr, "Failed loading plugin %s\n Reason: unspecified.\n\n", szfile); break;
	}

	return false;
}

//-----------------------------------------------------------------------------
void FinishLibrary()
{
	FEBioPluginManager* pPM = FEBioPluginManager::GetInstance();
	pPM->DeleteThis();
}

} // namespace febio
