#include "stdafx.h"
#include "fem.h"
#include "FEBioXML/XMLReader.h"
#include "FEBioMech/FESlidingInterface.h"
#include "FEBioMech/FETiedInterface.h"
#include "FEBioMix/FETiedBiphasicInterface.h"
#include "FEBioMech/FERigidWallInterface.h"
#include "FEBioMech/FEFacet2FacetSliding.h"
#include "FEBioMech/FEFacet2FacetTied.h"
#include "FEBioMix/FESlidingInterface2.h"
#include "FEBioMix/FESlidingInterface3.h"
#include "FEBioMech/FEPeriodicBoundary.h"
#include "FEBioMech/FESurfaceConstraint.h"
#include "FEBioMech/FESlidingInterfaceBW.h"
#include "FEBioMech/FERigidBody.h"
#include "FECore/log.h"
#include "plugin.h"
#include "Interrupt.h"
#include "console.h"

#ifdef WIN32
extern "C" void __cdecl omp_set_num_threads(int);
#else
extern "C" void omp_set_num_threads(int);
#endif

//-----------------------------------------------------------------------------
//! Constructor of the FEM class
//! Initializes default variables
//!
FEM::FEM()
{
	// User may interrupt run
	m_bInterruptable = true;
}

//-----------------------------------------------------------------------------
//! destructor of FEM class.
FEM::~FEM()
{
}

//-----------------------------------------------------------------------------
//! Reads the FEBio configuration file. This file contains some default settings.

bool FEM::Configure(const char *szfile)
{
	// open the configuration file
	XMLReader xml;
	if (xml.Open(szfile) == false)
	{
		fprintf(stderr, "FATAL ERROR: Failed reading FEBio configuration file %s.\n", szfile);
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
					if (tag == "linear_solver")
					{
						const char* szt = tag.AttributeValue("type");
						if      (strcmp(szt, "skyline"           ) == 0) m_nsolver = SKYLINE_SOLVER;
						else if (strcmp(szt, "psldlt"            ) == 0) m_nsolver = PSLDLT_SOLVER;
						else if (strcmp(szt, "superlu"           ) == 0) m_nsolver = SUPERLU_SOLVER;
						else if (strcmp(szt, "superlu_mt"        ) == 0) m_nsolver = SUPERLU_MT_SOLVER;
						else if (strcmp(szt, "pardiso"           ) == 0) m_nsolver = PARDISO_SOLVER;
						else if (strcmp(szt, "rcicg"             ) == 0) m_nsolver = RCICG_SOLVER;
						else if (strcmp(szt, "wsmp"              ) == 0) m_nsolver = WSMP_SOLVER;
					}
					else if (tag == "import")
					{
						const char* szfile = tag.szvalue();
						if (LoadPlugin(szfile) == false) throw XMLReader::InvalidValue(tag);
						printf("Plugin \"%s\" loaded successfully\n", szfile);
					}
					else if (tag == "import_folder")
					{
						const char* szfile = tag.szvalue();
						if (LoadPluginFolder(szfile) == false) throw XMLReader::InvalidTag(tag);
					}
					else if (tag == "omp_num_threads")
					{
						int n;
						tag.value(n);
						omp_set_num_threads(n);
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
			clog.printbox("FATAL ERROR", "Invalid version for FEBio configuration file.");
			return false;
		}
	}
	catch (XMLReader::XMLSyntaxError)
	{
		clog.printf("FATAL ERROR: Syntax error (line %d)\n", xml.GetCurrentLine());
		return false;
	}
	catch (XMLReader::InvalidTag e)
	{
		clog.printf("FATAL ERROR: unrecognized tag \"%s\" (line %d)\n", e.tag.m_sztag, e.tag.m_nstart_line);
		return false;
	}
	catch (XMLReader::InvalidAttributeValue e)
	{
		const char* szt = e.tag.m_sztag;
		const char* sza = e.szatt;
		const char* szv = e.szval;
		int l = e.tag.m_nstart_line;
		clog.printf("FATAL ERROR: unrecognized value \"%s\" for attribute \"%s.%s\" (line %d)\n", szv, szt, sza, l);
		return false;
	}
	catch (XMLReader::InvalidValue e)
	{
		clog.printf("FATAL ERROR: the value for tag \"%s\" is invalid (line %d)\n", e.tag.m_sztag, e.tag.m_nstart_line);
		return false;
	}
	catch (XMLReader::MissingAttribute e)
	{
		clog.printf("FATAL ERROR: Missing attribute \"%s\" of tag \"%s\" (line %d)\n", e.szatt, e.tag.m_sztag, e.tag.m_nstart_line);
		return false;
	}
	catch (XMLReader::UnmatchedEndTag e)
	{
		const char* sz = e.tag.m_szroot[e.tag.m_nlevel];
		clog.printf("FATAL ERROR: Unmatched end tag for \"%s\" (line %d)\n", sz, e.tag.m_nstart_line);
		return false;
	}
	catch (...)
	{
		clog.printf("FATAL ERROR: unrecoverable error (line %d)\n", xml.GetCurrentLine());
		return false;
	}

	xml.Close();

	return true;
}

//-----------------------------------------------------------------------------
//! This function is called in several places to see if the user requested to
//! pause the run. If so, the FEBio prompt appears and users can enter a command.
void FEM::CheckInterruption()
{
	if (m_bInterruptable)
	{
		Interruption itr;
		if (itr.m_bsig)
		{
			itr.m_bsig = false;
			itr.interrupt();
		}
	}
}
