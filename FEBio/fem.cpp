#include "stdafx.h"
#include "fem.h"
#include <string.h>
#include "Archive.h"
#include "XMLReader.h"
#include "FEFacet2FacetSliding.h"
#include "FESlidingInterface2.h"

//-----------------------------------------------------------------------------
//! Constructor of the FEM class
//! Initializes default variables
//!
FEM::FEM()
{
	// --- Analysis Data ---
	// we create at least one step
	m_pStep = new FEAnalysis(*this);
	m_Step.add(m_pStep);
	m_nStep = 0;
	m_nhex8 = FE_HEX;
	m_bsym_poro = true;	// use symmetric poro implementation

	m_ftime = 0;

	// --- Geometry Data ---
	m_nreq = 0;
	m_nrb = 0;
	m_nrm = 0;
	m_nrj = 0;

	m_bcontact = false;		// assume no contact

	m_bsymm = true;	// assume symmetric stiffness matrix

	m_psurf = new FESurface(&m_mesh);

	// --- Material Data ---
	// (nothing to initialize yet)

	// --- Load Curve Data ---
	// (nothing to initialize yet)

	// --- Boundary Condition Data ---
	for (int i=0; i<3; ++i)
	{
		m_BF[i].lc = -1;
		m_BF[i].s = 0.0;
	}

	// --- Direct Solver Data ---
	// set the skyline solver as default
	m_nsolver = SKYLINE_SOLVER;

	// However, if available use the PSLDLT solver instead
#ifdef PSLDLT
	m_nsolver = PSLDLT_SOLVER;
#endif

	m_neq = 0;
	m_bwopt = 0;

	// --- I/O-Data ---
	strcpy(m_szplot, "n3plot");
	strcpy(m_szlog , "n3log" );
	strcpy(m_szdump, "n3dump");

	m_sztitle[0] = 0;
	m_debug = false;
	m_pcb = 0;
	m_pcbd = 0;
}

///////////////////////////////////////////////////////////////////////////////
// FUNCTION : FEM::~FEM
// destructor of FEM class
//

FEM::~FEM()
{

}

///////////////////////////////////////////////////////////////////////////////
// FUNCTION: FEM::FEM(FEM& fem)
// copy constructor.
// The copy constructor and assignment operator are used for push/pop'ing.
// Note that not all data is copied. We only copy the data that is relevant
// for push/pop'ing
//
// TODO: apparently, the copy constructor never gets called. Maybe we should
// make it private.

FEM::FEM(FEM& fem)
{
	assert(false);
}

///////////////////////////////////////////////////////////////////////////////
// FUNCTION: FEM::operator = (FEM& fem)
// assignment operator
// The copy constructor and assignment operator are used for push/pop'ing.
// Note that not all data is copied. We only copy the data that is relevant
// for push/pop'ing
//

void FEM::operator =(FEM& fem)
{
	int i;

	// keep a pointer to the current analysis step
	// note that we do not keep the entire analysis history
	// since that would be waste of space and time
	// however, that does imply that for *this fem object
	// the m_Step array remains empty!
	m_pStep = fem.m_pStep;
	m_nStep = m_nStep;

	// copy the mesh
	m_mesh = fem.m_mesh;

	m_nrb = fem.m_nrb;
	m_RB = fem.m_RB;

	m_ftime = fem.m_ftime;

	// copy rigid joint data
	if (m_nrj == 0)
	{
		for (i=0; i<fem.m_nrj; ++i) m_RJ.add(new FERigidJoint(this));
		m_nrj = m_RJ.size();
	}
	assert(m_nrj == fem.m_nrj);
	for (i=0; i<m_nrj; ++i) m_RJ[i].ShallowCopy(fem.m_RJ[i]);

	// copy contact data
	if (ContactInterfaces() == 0)
	{
		FEContactInterface* pci;
		for (int i=0; i<fem.ContactInterfaces(); ++i)
		{
			switch (fem.m_CI[i].Type())
			{
			case FE_CONTACT_SLIDING:
				pci = new FESlidingInterface(this);
				break;
			case FE_FACET2FACET_SLIDING:
				pci = new FEFacet2FacetSliding(this);
				break;
			case FE_CONTACT_TIED:
				pci = new FETiedInterface(this);
				break;
			case FE_CONTACT_RIGIDWALL:
				pci = new FERigidWallInterface(this);
				break;
			case FE_CONTACT_SLIDING2:
				pci = new FESlidingInterface2(this);
			default:
				assert(false);
			}

			m_CI.add(pci);
		}
	}
	assert(ContactInterfaces() == fem.ContactInterfaces());
	for (i=0; i<ContactInterfaces(); ++i) m_CI[i].ShallowCopy(fem.m_CI[i]);
}

///////////////////////////////////////////////////////////////////////////////
// FUNCTION: FEM::DoCallback
// Call the callback function if there is one defined
//

void FEM::DoCallback()
{
	if (m_pcb) m_pcb(this, m_pcbd);
}

//-----------------------------------------------------------------------------
//! Return a pointer to the named variable

//! This function returns a pointer to a named variable. Currently, we only
//! support names of the form material_name.parameter name. The material
//! name is a user defined name for a material and the parameter name is the
//! predefined name of the variable
//! \todo perhaps I should use XPath to refer to material parameters ?

double* FEM::FindParameter(const char* szparam)
{
	int i, nmat;

	char szname[256];
	strcpy(szname, szparam);

	// get the material and parameter name
	char* ch = strchr((char*)szname, '.');
	if (ch == 0) return 0;
	*ch = 0;
	const char* szmat = szname;
	const char* szvar = ch+1;

	// find the material with the same name
	FEMaterial* pmat;

	for (i=0; i<Materials(); ++i)
	{
		pmat = GetMaterial(i);
		nmat = i;

		if (strcmp(szmat, pmat->GetName()) == 0)
		{
			break;
		}

		pmat = 0;
	}

	// make sure we found a material with the same name
	if (pmat == 0) return false;

	// if the variable is a vector, then we require an index
	char* szarg = strchr((char*) szvar, '[');
	int index = -1;
	if (szarg)
	{
		*szarg = 0; szarg++;
		const char* ch = strchr(szarg, ']');
		assert(ch);
		index = atoi(szarg);
	}

	// get the material's parameter list
	auto_ptr<FEParameterList> pl(pmat->GetParameterList());

	// find the parameter
	FEParam* pp = pl->Find(szvar);

	if (pp)
	{
		switch (pp->m_itype)
		{
		case FE_PARAM_DOUBLE:
			{
				assert(index<0);
				return &pp->value<double>();
			}
			break;
		case FE_PARAM_DOUBLEV:
			{
				assert((index >= 0) && (index < pp->m_ndim));
				return &(pp->pvalue<double>()[index]);
			}
			break;
		default:
			// other types are not supported yet
			assert(false);
			return 0;
		}
	}


	// the rigid bodies are dealt with differently
	for (i=0; i<m_nrb; ++i)
	{
		FERigidBody& rb = m_RB[i];

		if (rb.m_mat == nmat)
		{
			if (strcmp(szvar, "Fx") == 0) return &rb.m_Fr.x;
			if (strcmp(szvar, "Fy") == 0) return &rb.m_Fr.y;
			if (strcmp(szvar, "Fz") == 0) return &rb.m_Fr.z;
			if (strcmp(szvar, "Mx") == 0) return &rb.m_Mr.x;
			if (strcmp(szvar, "My") == 0) return &rb.m_Mr.y;
			if (strcmp(szvar, "Mz") == 0) return &rb.m_Mr.z;
		}
	}

	// oh, oh, we didn't find it
	return 0;
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

	// get the logfile
	Logfile& log = m_log;

	// loop over all child tags
	try
	{
		// Find the root element
		XMLTag tag;
		if (xml.FindTag("febio_config", tag) == false) return false;

		if (strcmp(tag.m_szatv[0], "1.0") == 0)
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
						else if (strcmp(szt, "wsmp"              ) == 0) m_nsolver = WSMP_SOLVER;
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
			log.printbox("FATAL ERROR", "Invalid version for FEBio configuration file.");
			return false;
		}
	}
	catch (XMLReader::XMLSyntaxError)
	{
		log.printf("FATAL ERROR: Syntax error (line %d)\n", xml.GetCurrentLine());
		return false;
	}
	catch (XMLReader::InvalidTag e)
	{
		log.printf("FATAL ERROR: unrecognized tag \"%s\" (line %d)\n", e.tag.m_sztag, e.tag.m_nstart_line);
		return false;
	}
	catch (XMLReader::InvalidAttributeValue e)
	{
		const char* szt = e.tag.m_sztag;
		const char* sza = e.szatt;
		const char* szv = e.szval;
		int l = e.tag.m_nstart_line;
		log.printf("FATAL ERROR: unrecognized value \"%s\" for attribute \"%s.%s\" (line %d)\n", szv, szt, sza, l);
		return false;
	}
	catch (XMLReader::InvalidValue e)
	{
		log.printf("FATAL ERROR: the value for tag \"%s\" is invalid (line %d)\n", e.tag.m_sztag, e.tag.m_nstart_line);
		return false;
	}
	catch (XMLReader::MissingAttribute e)
	{
		log.printf("FATAL ERROR: Missing attribute \"%s\" of tag \"%s\" (line %d)\n", e.szatt, e.tag.m_sztag, e.tag.m_nstart_line);
		return false;
	}
	catch (XMLReader::UnmatchedEndTag e)
	{
		const char* sz = e.tag.m_szroot[e.tag.m_nlevel];
		log.printf("FATAL ERROR: Unmatched end tag for \"%s\" (line %d)\n", sz, e.tag.m_nstart_line);
		return false;
	}
	catch (...)
	{
		log.printf("FATAL ERROR: unrecoverable error (line %d)\n", xml.GetCurrentLine());
		return false;
	}

	xml.Close();

	return true;
}

//-----------------------------------------------------------------------------

FEBoundaryCondition* FEM::FindBC(int nid)
{
	int i;
	for (i=0; i<m_DC.size(); ++i) if (m_DC[i].GetID() == nid) return &m_DC[i];

	for (i=0; i<m_FC.size(); ++i) if (m_FC[i].GetID() == nid) return &m_FC[i];

	for (i=0; i<m_PC.size(); ++i) if (m_PC[i].GetID() == nid) return &m_PC[i];

	for (i=0; i<m_RDC.size(); ++i) if (m_RDC[i].GetID() == nid) return &m_RDC[i];

	for (i=0; i<m_RFC.size(); ++i) if (m_RFC[i].GetID() == nid) return &m_RFC[i];

	return 0;
}
