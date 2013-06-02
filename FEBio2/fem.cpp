#include "stdafx.h"
#include "fem.h"
#include "FEBioXML/XMLReader.h"
#include "FEBioLib/FESlidingInterface.h"
#include "FEBioLib/FETiedInterface.h"
#include "FEBioLib/FETiedBiphasicInterface.h"
#include "FEBioLib/FERigidWallInterface.h"
#include "FEBioLib/FEFacet2FacetSliding.h"
#include "FEBioLib/FEFacet2FacetTied.h"
#include "FEBioLib/FESlidingInterface2.h"
#include "FEBioLib/FESlidingInterface3.h"
#include "FEBioLib/FEPeriodicBoundary.h"
#include "FEBioLib/FESurfaceConstraint.h"
#include "FEBioLib/FESlidingInterfaceBW.h"
#include "FECore/FERigidBody.h"
#include "FECore/log.h"
#include "plugin.h"
#include "Interrupt.h"
#include "console.h"

extern "C" void __cdecl omp_set_num_threads(int);

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
static FEM* pfem_copy = 0;

//-----------------------------------------------------------------------------
void FEM::PushState()
{
	if (pfem_copy == 0) pfem_copy = new FEM;
	pfem_copy->ShallowCopy(*this);
}

//-----------------------------------------------------------------------------
void FEM::PopState()
{
	assert(pfem_copy);
	ShallowCopy(*pfem_copy);
	delete pfem_copy;
	pfem_copy = 0;
}

//-----------------------------------------------------------------------------
// This function is used when pushing the FEM state data. Since we don't need
// to copy all the data, this function only copies the data that needs to be 
// restored for a running restart.

// TODO: Shallow copy nonlinear constraints
void FEM::ShallowCopy(FEM& fem)
{
	int i;

	// copy time data
	m_ftime = fem.m_ftime;

	// copy the mesh
	m_mesh = fem.GetMesh();

	// copy rigid body data
	if (m_Obj.empty())
	{
		for (int i=0; i<(int) fem.m_Obj.size();	++i)
		{
			FERigidBody* prb = new FERigidBody(this);
			m_Obj.push_back(prb);
		}
	}
	assert(m_Obj.size() == fem.m_Obj.size());
	for (i=0; i<(int) m_Obj.size(); ++i) m_Obj[i]->ShallowCopy(fem.m_Obj[i]);

	// copy contact data
	if (m_CI.empty())
	{
		FEContactInterface* pci;
		for (int i=0; i<fem.ContactInterfaces(); ++i)
		{
			switch (fem.m_CI[i]->Type())
			{
			case FE_CONTACT_SLIDING      : pci = new FESlidingInterface     (this); break;
			case FE_FACET2FACET_SLIDING  : pci = new FEFacet2FacetSliding   (this); break;
			case FE_CONTACT_TIED         : pci = new FETiedInterface        (this); break;
			case FE_CONTACT_RIGIDWALL    : pci = new FERigidWallInterface   (this); break;
			case FE_CONTACT_SLIDING2     : pci = new FESlidingInterface2    (this); break;
			case FE_PERIODIC_BOUNDARY    : pci = new FEPeriodicBoundary     (this); break;
			case FE_SURFACE_CONSTRAINT   : pci = new FESurfaceConstraint    (this); break;
			case FE_CONTACT_SLIDING3     : pci = new FESlidingInterface3    (this); break;
			case FE_CONTACT_TIED_BIPHASIC: pci = new FETiedBiphasicInterface(this); break;
			case FE_CONTACT_SLIDINGBW    : pci = new FESlidingInterfaceBW   (this); break;
			case FE_FACET2FACET_TIED     : pci = new FEFacet2FacetTied      (this); break;
			default:
				assert(false);
			}
			m_CI.push_back(pci);
		}
	}
	assert(ContactInterfaces() == fem.ContactInterfaces());
	for (i=0; i<ContactInterfaces(); ++i) m_CI[i]->ShallowCopy(*fem.m_CI[i]);
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
