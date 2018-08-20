#include "stdafx.h"
#include "FECore/FECore.h"
#include "NumCore/NumCore.h"
#include "FEBioMech/FEBioMech.h"
#include "FEBioMix/FEBioMix.h"
#include "FEBioOpt/FEBioOpt.h"
#include "FEBioFluid/FEBioFluid.h"
#include <FEBioFluid/FEBioFSI.h>
#include <FEBioTest/FEBioTest.h>
#include <FEBioXML/XMLReader.h>
#include <FEBioXML/xmltool.h>
#include <FECore/FEModel.h>
#include <FECore/log.h>
#include "febio.h"
#include "plugin.h"
#include <map>

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
// import all modules
void InitLibrary()
{
	FECore::InitModule();
	NumCore::InitModule();
	FEBioMech::InitModule();
	FEBioMix::InitModule();
	FEBioOpt::InitModule();
	FEBioFluid::InitModule();
	FEBioFSI::InitModule();
	FEBioTest::InitModule();
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
						int nsolver = -1;
						if      (strcmp(szt, "skyline"           ) == 0) FECoreKernel::SetDefaultSolver(nsolver = SKYLINE_SOLVER   );
						else if (strcmp(szt, "psldlt"            ) == 0) FECoreKernel::SetDefaultSolver(nsolver = PSLDLT_SOLVER    );
						else if (strcmp(szt, "superlu"           ) == 0) FECoreKernel::SetDefaultSolver(nsolver = SUPERLU_SOLVER   );
						else if (strcmp(szt, "superlu_mt"        ) == 0) FECoreKernel::SetDefaultSolver(nsolver = SUPERLU_MT_SOLVER);
						else if (strcmp(szt, "pardiso"           ) == 0) FECoreKernel::SetDefaultSolver(nsolver = PARDISO_SOLVER   );
						else if (strcmp(szt, "rcicg"             ) == 0) FECoreKernel::SetDefaultSolver(nsolver = RCICG_SOLVER     );
						else if (strcmp(szt, "fgmres"            ) == 0) FECoreKernel::SetDefaultSolver(nsolver = FGMRES_SOLVER    );
						else if (strcmp(szt, "fgmres_ilut"       ) == 0) FECoreKernel::SetDefaultSolver(nsolver = FGMRES_ILUT_SOLVER);
						else if (strcmp(szt, "fgmres_ilu0"       ) == 0) FECoreKernel::SetDefaultSolver(nsolver = FGMRES_ILU0_SOLVER);
						else if (strcmp(szt, "wsmp"              ) == 0) FECoreKernel::SetDefaultSolver(nsolver = WSMP_SOLVER      );
						else if (strcmp(szt, "bipn"              ) == 0) FECoreKernel::SetDefaultSolver(nsolver = BIPN_SOLVER      );
						else if (strcmp(szt, "hypre_gmres"       ) == 0) FECoreKernel::SetDefaultSolver(nsolver = HYPRE_GMRES      );
						else if (strcmp(szt, "stokes"            ) == 0) FECoreKernel::SetDefaultSolver(nsolver = STOKES_SOLVER    );
						else if (strcmp(szt, "cg_stokes"         ) == 0) FECoreKernel::SetDefaultSolver(nsolver = CG_STOKES_SOLVER );
						else if (strcmp(szt, "schur"             ) == 0) FECoreKernel::SetDefaultSolver(nsolver = SCHUR_SOLVER     );
						else { fprintf(stderr, "Invalid linear solver\n"); return false; }

						if (tag.isleaf() == false)
						{
							FECoreKernel& fecore = FECoreKernel::GetInstance();
							FELinearSolverFactory* fac = fecore.FindLinearSolverFactory(nsolver);
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
							}
							while (!tag.isend());
						}
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
						omp_set_num_threads(n);
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
	default:
		fprintf(stderr, "Failed loading plugin %s\n Reason: unspecified.\n\n", szfile); break;
	}

	return false;
}

//-----------------------------------------------------------------------------
void SetOMPThreads(int n)
{
	omp_set_num_threads(n);
}

//-----------------------------------------------------------------------------
void FinishLibrary()
{
	FEBioPluginManager* pPM = FEBioPluginManager::GetInstance();
	pPM->DeleteThis();
}

} // namespace febio
