// FEBioImport.cpp: implementation of the FEFEBioImport class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FEBioImport.h"
#include "FEBioModuleSection.h"
#include "FEBioControlSection.h"
#include "FEBioGlobalsSection.h"
#include "FEBioMaterialSection.h"
#include "FEBioGeometrySection.h"
#include "FEBioBoundarySection.h"
#include "FEBioLoadsSection.h"
#include "FEBioContactSection.h"
#include "FEBioConstraintsSection.h"
#include "FEBioInitialSection.h"
#include "FEBioLoadDataSection.h"
#include "FEBioOutputSection.h"
#include "FEBioStepSection.h"
#include <FEBioLib/DataStore.h>
#include <FECore/log.h>
#include <string.h>
#include <FECore/Image.h>

//-----------------------------------------------------------------------------
FEModel* FEBioFileSection::GetFEModel() { return m_pim->GetFEModel(); }
FEAnalysisStep* FEBioFileSection::GetStep() { return dynamic_cast<FEAnalysisStep*>(m_pim->GetStep()); }

//-----------------------------------------------------------------------------
FEBioFileSectionMap::~FEBioFileSectionMap()
{
	// clear the map
	FEBioFileSectionMap::iterator is;
	for (is = begin(); is != end(); ++is)
	{
		FEBioFileSection* ps = is->second; delete ps;
	}
	clear();
}

//=============================================================================
//  The FEBioImport class imports an XML formatted FEBio input file.
//  The actual file is parsed using the XMLReader class.
//
bool FEFEBioImport::Load(FEModel& fem, const char* szfile)
{
	// Open the XML file
	XMLReader xml;
	if (xml.Open(szfile) == false) return errf("FATAL ERROR: Failed opening input file %s\n\n", szfile);

	// keep a pointer to the fem object
	m_pfem = &fem;

	// keep a pointer to the mesh
	m_pMesh = &fem.GetMesh();

	// intialize some variables
	m_pStep = 0;	// zero step pointer
	m_nsteps = 0;	// reset step section counter
	m_nstep_type = FE_SOLID;	// default step type in case Module section is not defined
	m_nversion = -1;
	m_szdmp[0] = 0;
	m_szlog[0] = 0;
	m_szplt[0] = 0;

	// default element type
	m_ntet4  = ET_TETG1;
	m_nhex8  = FE_HEX8G8;
	m_nut4   = FE_TET4G1;
	m_ntet10 = FE_TET10G4;
	m_ntri6  = FE_TRI6G3;
	m_ntri3  = FE_TRI3G3;

	// 3-field formulation on by default
	m_b3field = true;

	// Find the root element
	XMLTag tag;
	try
	{
		if (xml.FindTag("febio_spec", tag) == false) return errf("FATAL ERROR: febio_spec tag was not found. This is not a valid input file.\n\n");
	}
	catch (...)
	{
		clog.printf("An error occured while finding the febio_spec tag.\nIs this a valid FEBio input file?\n\n");
		return false;
	}

	// parse the file
	try
	{
		// get the version number
		ParseVersion(tag);

		// FEBio2 only supports file version 1.1 and up
		if ((m_nversion != 0x0101) &&
			(m_nversion != 0x0102) && 
			(m_nversion != 0x0200)) throw InvalidVersion();

		// define the file structure
		FEBioFileSectionMap map;
		map["Module"     ] = new FEBioModuleSection     (this);
		map["Control"    ] = new FEBioControlSection    (this);
		map["Material"   ] = new FEBioMaterialSection   (this);
		map["Geometry"   ] = new FEBioGeometrySection   (this);
		map["Boundary"   ] = new FEBioBoundarySection   (this);
		map["Initial"    ] = new FEBioInitialSection    (this);
		map["LoadData"   ] = new FEBioLoadDataSection   (this);
		map["Globals"    ] = new FEBioGlobalsSection    (this);
		map["Output"     ] = new FEBioOutputSection     (this);
		map["Constraints"] = new FEBioConstraintsSection(this);
		map["Step"       ] = new FEBioStepSection       (this);

		// version 1.2 only
		if (m_nversion >= 0x0102)
		{
			map["Loads"] = new FEBioLoadsSection(this);
		}

		// version 2.0 only!
		if (m_nversion >= 0x0200)
		{
			map["Contact"] = new FEBioContactSection(this);
		}

		// parse the file
		++tag;
		do
		{
			// try to find a section parser
			FEBioFileSectionMap::iterator is = map.find(tag.Name());

			// make sure we found a section reader
			if (is == map.end()) throw XMLReader::InvalidTag(tag);

			// see if the file has the "from" attribute (for version 2.0 and up)
			if (m_nversion >= 0x0200)
			{
				const char* szinc = tag.AttributeValue("from", true);
				if (szinc)
				{
					// make sure this is a leaf
					if (tag.isleaf() == false) return errf("FATAL ERROR: included sections may not have child sections.\n\n");

					// read this section from an included file.
					XMLReader xml2;
					if (xml2.Open(szinc) == false) return errf("FATAL ERROR: failed opening input file %s\n\n", szinc);

					// find the febio_spec tag
					XMLTag tag2;
					if (xml2.FindTag("febio_spec", tag2) == false) return errf("FATAL ERROR: febio_spec tag was not found. This is not a valid input file.\n\n");

					// find the section we are looking for
					if (xml2.FindTag(tag.Name(), tag2) == false) return errf("FATAL ERROR: Couldn't find %s section in file %s.\n\n", tag.Name(), szinc);

					// parse the section
					is->second->Parse(tag2);
				}
				else is->second->Parse(tag);
			}
			else
			{
				// parse the section
				is->second->Parse(tag);
			}

			// go to the next tag
			++tag;
		}
		while (!tag.isend());
	}
	// --- XML Reader Exceptions ---
	catch (XMLReader::XMLSyntaxError)
	{
		clog.printf("FATAL ERROR: Syntax error (line %d)\n", xml.GetCurrentLine());
		return false;
	}
	catch (XMLReader::InvalidAttributeValue e)
	{
		const char* szt = e.tag.m_sztag;
		const char* sza = e.szatt;
		const char* szv = e.szval;
		int l = e.tag.m_nstart_line;
		clog.printf("FATAL ERROR: invalid value \"%s\" for attribute \"%s.%s\" (line %d)\n", szv, szt, sza, l);
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
	// --- FEBio Exceptions ---
	catch (InvalidVersion)
	{
		clog.printbox("FATAL ERROR", "Invalid version for FEBio specification.");
		return false;
	}
	catch (InvalidMaterial e)
	{
		clog.printbox("FATAL ERROR:", "Element %d has an invalid material type.", e.m_nel);
		return false;
	}
	catch (XMLReader::InvalidTag e)
	{
		clog.printf("FATAL ERROR: unrecognized tag \"%s\" (line %d)\n", e.tag.m_sztag, e.tag.m_nstart_line);
		return false;
	}
	catch (InvalidDomainType)	
	{
		clog.printf("Fatal Error: Invalid domain type\n");
		return false;
	}
	catch (FailedCreatingDomain)
	{
		clog.printf("Fatal Error: Failed creating domain\n");
		return false;
	}
	catch (InvalidElementType)
	{
		clog.printf("Fatal Error: Invalid element type\n");
		return false;
	}
	catch (FailedLoadingPlugin e)
	{
		clog.printf("Fatal Error: failed loading plugin %s\n", e.FileName());
		return false;
	}
	catch (UnknownDataField e)
	{
		clog.printf("Fatal Error: \"%s\" is not a valid field variable name (line %d)\n", e.m_szdata, xml.GetCurrentLine()-1);
		return false;
	}
	catch (DuplicateMaterialSection)
	{
		clog.printf("Fatal Error: Material section has already been defined (line %d).\n", xml.GetCurrentLine()-1);
		return false;
	}
	// --- Unknown exceptions ---
	catch (...)
	{
		clog.printf("FATAL ERROR: unrecoverable error (line %d)\n", xml.GetCurrentLine());
		return false;
	}

	// close the XML file
	xml.Close();

	// we're done!
	return true;
}

//-----------------------------------------------------------------------------
FEAnalysis* FEFEBioImport::GetStep()
{
	if (m_pStep == 0)
	{
		m_pStep = CreateNewStep();
		m_pfem->AddStep(m_pStep);	
		if (m_pfem->Steps() == 1) 
		{
			m_pfem->SetCurrentStep(m_pStep);
			m_pfem->m_nStep = 0;
		}
	}
	return m_pStep;
}

//-----------------------------------------------------------------------------
FEAnalysis* FEFEBioImport::CreateNewStep()
{
	FEAnalysis* pstep = 0;
	switch (m_nstep_type)
	{
	case FE_SOLID         : pstep = new FESolidAnalysis         (*m_pfem); break;
	case FE_EXPLICIT_SOLID: pstep = new FEExplicitSolidAnalysis (*m_pfem); break;
	case FE_LINEAR_SOLID  : pstep = new FELinearSolidAnalysis   (*m_pfem); break;
	case FE_BIPHASIC      : pstep = new FEBiphasicAnalysis      (*m_pfem); break;
	case FE_POROSOLUTE    : pstep = new FEBiphasicSoluteAnalysis(*m_pfem); break;
	case FE_HEAT          : pstep = new FEHeatTransferAnalysis  (*m_pfem); break;
	case FE_HEAT_SOLID    : pstep = new FEThermoElasticAnalysis (*m_pfem); break;
	default:
		assert(false);
	}
	return pstep;
}

//-----------------------------------------------------------------------------
//! This function parses the febio_spec tag for the version number
void FEFEBioImport::ParseVersion(XMLTag &tag)
{
	const char* szv = tag.AttributeValue("version");
	assert(szv);
	int n1, n2;
	int nr = sscanf(szv, "%d.%d", &n1, &n2);
	if (nr != 2) throw InvalidVersion();
	if ((n1 < 1) || (n1 > 0xFF)) throw InvalidVersion();
	if ((n2 < 0) || (n2 > 0xFF)) throw InvalidVersion();
	m_nversion = (n1 << 8) + n2;
}

//-----------------------------------------------------------------------------
//! This function parese a parameter list
bool FEFEBioImport::ReadParameter(XMLTag& tag, FEParameterList& pl, const char* szparam)
{
	// see if we can find this parameter
	FEParam* pp = pl.Find((szparam == 0 ? tag.Name() : szparam));
	if (pp)
	{
		switch (pp->m_itype)
		{
		case FE_PARAM_DOUBLE : tag.value(pp->value<double>() ); break;
		case FE_PARAM_INT    : tag.value(pp->value<int   >() ); break;
		case FE_PARAM_BOOL   : tag.value(pp->value<bool  >() ); break;
		case FE_PARAM_VEC3D  : tag.value(pp->value<vec3d >() ); break;
		case FE_PARAM_STRING : tag.value(pp->cvalue() ); break;
		case FE_PARAM_INTV   : tag.value(pp->pvalue<int   >(), pp->m_ndim); break;
		case FE_PARAM_DOUBLEV: tag.value(pp->pvalue<double>(), pp->m_ndim); break;
		case FE_PARAM_IMAGE_3D:
			{
				const char* szfile = tag.AttributeValue("file");
				++tag;
				int n[3] = {0};
				do
				{
					if (tag == "size") tag.value(n, 3);
					else throw XMLReader::InvalidTag(tag);
					++tag;
				}
				while (!tag.isend());
				Image& im = pp->value<Image>();
				im.Create(n[0], n[1], n[2]);
				if (im.Load(szfile) == false) throw XMLReader::InvalidValue(tag);
			}
			break;
		default:
			assert(false);
			return false;
		}

		int nattr = tag.m_natt;
		for (int i=0; i<nattr; ++i)
		{
			const char* szat = tag.m_att[i].m_szatt;
			if (strcmp(szat, "lc") == 0)
			{
				int lc = atoi(tag.m_att[i].m_szatv)-1;
				if (lc < 0) throw XMLReader::InvalidAttributeValue(tag, szat, tag.m_att[i].m_szatv);
				pp->m_nlc = lc;
				switch (pp->m_itype)
				{
				case FE_PARAM_DOUBLE: pp->m_scl = pp->value<double>(); break;
				}
			}
			// This is not true. Parameters can have attributes that are used for other purposed. E.g. The local fiber option.
//			else clog.printf("WARNING: attribute \"%s\" of parameter \"%s\" ignored (line %d)\n", szat, tag.Name(), tag.m_ncurrent_line-1);
		}

		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
void FEFEBioImport::ReadList(XMLTag& tag, vector<int>& l)
{
	// make sure the list is empty
	l.clear();

	// get a pointer to the value
	const char* sz = tag.szvalue();

	// parse the string
	const char* ch;
	do
	{
		int n0, n1, nn;
		int nread = sscanf(sz, "%d:%d:%d", &n0, &n1, &nn);
		switch (nread)
		{
		case 1:
			n1 = n0;
			nn = 1;
			break;
		case 2:
			nn = 1;
			break;
		}

		for (int i=n0; i<=n1; i += nn) l.push_back(i);

		ch = strchr(sz, ',');
		if (ch) sz = ch+1;
	}
	while (ch != 0);
}
