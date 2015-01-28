// FEBioImport.cpp: implementation of the FEFEBioImport class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FEBioImport.h"
#include "FEBioParametersSection.h"
#include "FEBioIncludeSection.h"
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
#include "FEBioDiscreteSection.h"
#include "FECore/DataStore.h"
#include "FECore/log.h"
#include "FECore/Image.h"
#include "FECore/FEModel.h"
#include "FECore/FECoreKernel.h"
#include "FECore/DOFS.h"
#include <string.h>

//-----------------------------------------------------------------------------
FEModel* FEBioFileSection::GetFEModel() { return m_pim->GetFEModel(); }
FEAnalysis* FEBioFileSection::GetStep() { return m_pim->GetStep(); }

//-----------------------------------------------------------------------------
FEBioFileSectionMap::~FEBioFileSectionMap()
{
	Clear();
}

//-----------------------------------------------------------------------------
void FEBioFileSectionMap::Clear()
{
	// clear the map
	FEBioFileSectionMap::iterator is;
	for (is = begin(); is != end(); ++is)
	{
		FEBioFileSection* ps = is->second; delete ps;
	}
	clear();
}

//-----------------------------------------------------------------------------
FEFEBioImport::PlotVariable::PlotVariable(const FEFEBioImport::PlotVariable& pv)
{
	strcpy(m_szvar, pv.m_szvar);
	m_item = pv.m_item;
}

FEFEBioImport::PlotVariable::PlotVariable(const char* szvar, vector<int>& item)
{
	strcpy(m_szvar, szvar);
	m_item = item;
}


//-----------------------------------------------------------------------------
void FEFEBioImport::ClearParams()
{
	m_Param.clear();
}

//-----------------------------------------------------------------------------
FEFEBioImport::XMLParam* FEFEBioImport::FindParameter(const char* sz)
{
	for (size_t i=0; i<m_Param.size(); ++i)
	{
		XMLParam& p = m_Param[i];
		if (strcmp(p.m_szname, sz) == 0) return &p;
	}
	return 0;
}

//-----------------------------------------------------------------------------
void FEFEBioImport::AddParameter(const char* szname, const char* szval)
{
	XMLParam p;
	strcpy(p.m_szname, szname);
	strcpy(p.m_szval , szval);
	m_Param.push_back(p);
}

//-----------------------------------------------------------------------------
const char* FEFEBioImport::get_value_string(XMLTag& tag)
{
	const char* sz = tag.szvalue();
	if (sz[0]=='@')
	{
		XMLParam* p = FindParameter(sz+1);
		if (p==0) throw XMLReader::InvalidValue(tag);
		sz = p->m_szval;
	}
	return sz;
}

//-----------------------------------------------------------------------------
void FEFEBioImport::value(XMLTag& tag, int& n)
{
	const char* sz = get_value_string(tag);
	n = atoi(sz);
}

//-----------------------------------------------------------------------------
void FEFEBioImport::value(XMLTag& tag, double& g)
{
	const char* sz = get_value_string(tag);
	g = atof(sz);
}

//-----------------------------------------------------------------------------
void FEFEBioImport::value(XMLTag& tag, bool& b)
{
	const char* sz = get_value_string(tag);
	int n=0; 
	sscanf(sz, "%d", &n); 
	b = (n != 0); 
}

//-----------------------------------------------------------------------------
void FEFEBioImport::value(XMLTag& tag, vec3d& v)
{
	const char* sz = get_value_string(tag);
	int n = sscanf(sz, "%lg,%lg,%lg", &v.x, &v.y, &v.z);
	if (n != 3) throw XMLReader::XMLSyntaxError();
}

//-----------------------------------------------------------------------------
void FEFEBioImport::value(XMLTag& tag, mat3d& m)
{
	const char* sz = get_value_string(tag);
	double xx, xy, xz, yx, yy, yz, zx, zy, zz;
	int n = sscanf(sz, "%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg", &xx, &xy, &xz, &yx, &yy, &yz, &zx, &zy, &zz);
	if (n != 9) throw XMLReader::XMLSyntaxError();
	m = mat3d(xx, xy, xz, yx, yy, yz, zx, zy, zz);
}

//-----------------------------------------------------------------------------
void FEFEBioImport::value(XMLTag& tag, mat3ds& m)
{
	const char* sz = get_value_string(tag);
	double x, y, z, xy, yz, xz;
	int n = sscanf(sz, "%lg,%lg,%lg,%lg,%lg,%lg", &x, &y, &z, &xy, &yz, &xz);
	if (n != 6) throw XMLReader::XMLSyntaxError();
	m = mat3ds(x, y, z, xy, yz, xz);
}

//-----------------------------------------------------------------------------
void FEFEBioImport::value(XMLTag& tag, char* szstr)
{
	const char* sz = get_value_string(tag);
	strcpy(szstr, sz); 
}

//-----------------------------------------------------------------------------
int FEFEBioImport::value(XMLTag& tag, int* pi, int n)
{
	const char* sz = get_value_string(tag);
	int nr = 0;
	for (int i=0; i<n; ++i)
	{
		const char* sze = strchr(sz, ',');

		pi[i] = atoi(sz);
		nr++;

		if (sze) sz = sze+1;
		else break;
	}
	return nr;
}

//-----------------------------------------------------------------------------
int FEFEBioImport::value(XMLTag& tag, double* pf, int n)
{
	const char* sz = get_value_string(tag);
	int nr = 0;
	for (int i=0; i<n; ++i)
	{
		const char* sze = strchr(sz, ',');

		pf[i] = atof(sz);
		nr++;

		if (sze) sz = sze+1;
		else break;
	}
	return nr;
}

//-----------------------------------------------------------------------------
FEFEBioImport::FEFEBioImport()
{
	// define the file structure
	m_map["Module"     ] = new FEBioModuleSection     (this);
	m_map["Control"    ] = new FEBioControlSection    (this);
	m_map["Material"   ] = new FEBioMaterialSection   (this);
	m_map["Geometry"   ] = new FEBioGeometrySection   (this);
	m_map["Boundary"   ] = new FEBioBoundarySection   (this);
	m_map["Loads"      ] = new FEBioLoadsSection      (this);
	m_map["Initial"    ] = new FEBioInitialSection    (this);
	m_map["LoadData"   ] = new FEBioLoadDataSection   (this);
	m_map["Globals"    ] = new FEBioGlobalsSection    (this);
	m_map["Output"     ] = new FEBioOutputSection     (this);
	m_map["Constraints"] = new FEBioConstraintsSection(this);
	m_map["Step"       ] = new FEBioStepSection       (this);

	// version 2.0 only!
	m_map["Parameters" ] = new FEBioParametersSection (this);
	m_map["Include"    ] = new FEBioIncludeSection    (this);
	m_map["Contact"    ] = new FEBioContactSection(this);
	m_map["Discrete"   ] = new FEBioDiscreteSection(this);
}

//=============================================================================
//  The FEBioImport class imports an XML formatted FEBio input file.
//  The actual file is parsed using the XMLReader class.
//
bool FEFEBioImport::Load(FEModel& fem, const char* szfile)
{
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

	// plot output
	m_szplot_type[0] = 0;
	m_plot.clear();
	m_nplot_compression = 0;

	// default element type
	m_ntet4  = FE_TET4G1;
	m_nhex8  = FE_HEX8G8;
	m_ntet10 = FE_TET10G4;
	m_ntet15 = FE_TET15G15;
	m_ntri6  = FE_TRI6G3;
	m_ntri3  = FE_TRI3G3;
	m_ntri7  = FE_TRI7G7;

	// 3-field formulation on by default
	m_b3field = true;

	// UT4 formulation off by default
	m_but4 = false;

	// Reset degrees of freedom (TODO: Can I do this elsewhere?)
    DOFS& fedofs = *DOFS::GetInstance();
	fedofs.Reset();

	// extract the path
	strcpy(m_szpath, szfile);
	char* ch = strrchr(m_szpath, '\\');
	if (ch==0) ch = strrchr(m_szpath, '/');
	if (ch==0) m_szpath[0] = 0; else *(ch+1)=0;

	// clear the parameters 
	ClearParams();

	// read the file
	return ReadFile(szfile);
}

//-----------------------------------------------------------------------------
bool FEFEBioImport::ReadFile(const char* szfile)
{
	// Open the XML file
	XMLReader xml;
	if (xml.Open(szfile) == false) return errf("FATAL ERROR: Failed opening input file %s\n\n", szfile);

	// Find the root element
	XMLTag tag;
	try
	{
		if (xml.FindTag("febio_spec", tag) == false) return errf("FATAL ERROR: febio_spec tag was not found. This is not a valid input file.\n\n");
	}
	catch (...)
	{
		felog.printf("An error occured while finding the febio_spec tag.\nIs this a valid FEBio input file?\n\n");
		return false;
	}

	// parse the file
	try
	{
		// get the version number
		ParseVersion(tag);

		// FEBio2 only supports file version 1.2 and 2.0
		if ((m_nversion != 0x0102) && 
			(m_nversion != 0x0200)) throw InvalidVersion();

		// parse the file
		++tag;
		do
		{
			// try to find a section parser
			FEBioFileSectionMap::iterator is = m_map.find(tag.Name());

			// make sure we found a section reader
			if (is == m_map.end()) throw XMLReader::InvalidTag(tag);

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
		felog.printf("FATAL ERROR: Syntax error (line %d)\n", xml.GetCurrentLine());
		return false;
	}
	catch (XMLReader::InvalidAttributeValue e)
	{
		const char* szt = e.tag.m_sztag;
		const char* sza = e.szatt;
		const char* szv = e.szval;
		int l = e.tag.m_nstart_line;
		felog.printf("FATAL ERROR: invalid value \"%s\" for attribute \"%s.%s\" (line %d)\n", szv, szt, sza, l);
		return false;
	}
	catch (XMLReader::InvalidValue e)
	{
		felog.printf("FATAL ERROR: the value for tag \"%s\" is invalid (line %d)\n", e.tag.m_sztag, e.tag.m_nstart_line);
		return false;
	}
	catch (XMLReader::MissingAttribute e)
	{
		felog.printf("FATAL ERROR: Missing attribute \"%s\" of tag \"%s\" (line %d)\n", e.szatt, e.tag.m_sztag, e.tag.m_nstart_line);
		return false;
	}
	catch (XMLReader::UnmatchedEndTag e)
	{
		const char* sz = e.tag.m_szroot[e.tag.m_nlevel];
		felog.printf("FATAL ERROR: Unmatched end tag for \"%s\" (line %d)\n", sz, e.tag.m_nstart_line);
		return false;
	}
	// --- FEBio Exceptions ---
	catch (InvalidVersion)
	{
		felog.printbox("FATAL ERROR", "Invalid version for FEBio specification.");
		return false;
	}
	catch (InvalidMaterial e)
	{
		felog.printbox("FATAL ERROR:", "Element %d has an invalid material type.", e.m_nel);
		return false;
	}
	catch (XMLReader::InvalidTag e)
	{
		felog.printf("FATAL ERROR: unrecognized tag \"%s\" (line %d)\n", e.tag.m_sztag, e.tag.m_nstart_line);
		return false;
	}
	catch (InvalidDomainType)	
	{
		felog.printf("Fatal Error: Invalid domain type\n");
		return false;
	}
	catch (InvalidDomainMaterial)
	{
		felog.printf("Fatal Error: Invalid domain material (line %d)\n", xml.GetCurrentLine());
		return false;
	}
	catch (FailedCreatingDomain)
	{
		felog.printf("Fatal Error: Failed creating domain\n");
		return false;
	}
	catch (InvalidElementType)
	{
		felog.printf("Fatal Error: Invalid element type\n");
		return false;
	}
	catch (FailedLoadingPlugin e)
	{
		felog.printf("Fatal Error: failed loading plugin %s\n", e.FileName());
		return false;
	}
	catch (UnknownDataField e)
	{
		felog.printf("Fatal Error: \"%s\" is not a valid field variable name (line %d)\n", e.m_szdata, xml.GetCurrentLine()-1);
		return false;
	}
	catch (DuplicateMaterialSection)
	{
		felog.printf("Fatal Error: Material section has already been defined (line %d).\n", xml.GetCurrentLine()-1);
		return false;
	}
	// --- Unknown exceptions ---
	catch (...)
	{
		felog.printf("FATAL ERROR: unrecoverable error (line %d)\n", xml.GetCurrentLine());
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
	case FE_SOLID         : pstep = fecore_new<FEAnalysis>(FEANALYSIS_ID, "solid"          , m_pfem); break;
	case FE_SOLID2        : pstep = fecore_new<FEAnalysis>(FEANALYSIS_ID, "solid"          , m_pfem); break;
	case FE_CG_SOLID      : pstep = fecore_new<FEAnalysis>(FEANALYSIS_ID, "solid"          , m_pfem); break;
	case FE_EXPLICIT_SOLID: pstep = fecore_new<FEAnalysis>(FEANALYSIS_ID, "explicit-solid" , m_pfem); break;
	case FE_LINEAR_SOLID  : pstep = fecore_new<FEAnalysis>(FEANALYSIS_ID, "linear-solid"   , m_pfem); break;
	case FE_BIPHASIC      : pstep = fecore_new<FEAnalysis>(FEANALYSIS_ID, "biphasic"       , m_pfem); break;
	case FE_POROSOLUTE    : pstep = fecore_new<FEAnalysis>(FEANALYSIS_ID, "biphasic-solute", m_pfem); break;
	case FE_MULTIPHASIC   : pstep = fecore_new<FEAnalysis>(FEANALYSIS_ID, "multiphasic"    , m_pfem); break;
	case FE_HEAT          : pstep = fecore_new<FEAnalysis>(FEANALYSIS_ID, "heat transfer"  , m_pfem); break;
//	case FE_HEAT_SOLID    : pstep = new FEThermoElasticAnalysis (*m_pfem); break;
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
		case FE_PARAM_DOUBLE : value(tag, pp->value<double>()); break;
		case FE_PARAM_INT    : value(tag, pp->value<int   >()); break;
		case FE_PARAM_BOOL   : value(tag, pp->value<bool  >()); break;
		case FE_PARAM_VEC3D  : value(tag, pp->value<vec3d >()); break;
		case FE_PARAM_MAT3D  : value(tag, pp->value<mat3d >()); break;
		case FE_PARAM_MAT3DS : value(tag, pp->value<mat3ds>()); break;
		case FE_PARAM_STRING : value(tag, pp->cvalue()); break;
		case FE_PARAM_INTV   : value(tag, pp->pvalue<int   >(), pp->m_ndim); break;
		case FE_PARAM_DOUBLEV: value(tag, pp->pvalue<double>(), pp->m_ndim); break;
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

				// see if we need to pre-pend a path
				char szin[512];
				strcpy(szin, szfile);
				char* ch = strrchr(szin, '\\');
				if (ch==0) ch = strrchr(szin, '/');
				if (ch==0)
				{
					// pre-pend the name with the input path
					sprintf(szin, "%s%s", m_szpath, szfile);
				}

				// Try to load the image file
				if (im.Load(szin) == false) throw XMLReader::InvalidValue(tag);
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
			if (pl.GetContainer()->SetParameterAttribute(*pp, szat, tag.m_att[i].m_szatv) == false)
			{
				// If we get here, the container did not understand the attribute.
				// If the attribute is a "lc", we interpret it as a load curve
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
/*				else 
				{
					throw XMLReader::InvalidAttributeValue(tag, szat, tag.m_att[i].m_szatv);
				}
*/			}
			// This is not true. Parameters can have attributes that are used for other purposed. E.g. The local fiber option.
//			else felog.printf("WARNING: attribute \"%s\" of parameter \"%s\" ignored (line %d)\n", szat, tag.Name(), tag.m_ncurrent_line-1);
		}

		// give the parameter container a chance to do additional processing
		pl.GetContainer()->SetParameter(*pp);

		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
//! This function parses a parameter list
bool FEFEBioImport::ReadParameter(XMLTag& tag, FECoreBase* pc, const char* szparam)
{
	FEParameterList& pl = pc->GetParameterList();

	// see if we can find this parameter
	FEParam* pp = pl.Find((szparam == 0 ? tag.Name() : szparam));
	if (pp)
	{
		switch (pp->m_itype)
		{
		case FE_PARAM_DOUBLE : value(tag, pp->value<double>() ); break;
		case FE_PARAM_INT    : value(tag, pp->value<int   >() ); break;
		case FE_PARAM_BOOL   : value(tag, pp->value<bool  >() ); break;
		case FE_PARAM_VEC3D  : value(tag, pp->value<vec3d >() ); break;
		case FE_PARAM_STRING : value(tag, pp->cvalue() ); break;
		case FE_PARAM_INTV   : value(tag, pp->pvalue<int   >(), pp->m_ndim); break;
		case FE_PARAM_DOUBLEV: value(tag, pp->pvalue<double>(), pp->m_ndim); break;
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
			if (pl.GetContainer()->SetParameterAttribute(*pp, szat, tag.m_att[i].m_szatv) == false)
			{
				// If we get here, the container did not understand the attribute.
				// If the attribute is a "lc", we interpret it as a load curve
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
/*				else 
				{
					throw XMLReader::InvalidAttributeValue(tag, szat, tag.m_att[i].m_szatv);
				}
*/			}
			// This is not true. Parameters can have attributes that are used for other purposed. E.g. The local fiber option.
//			else felog.printf("WARNING: attribute \"%s\" of parameter \"%s\" ignored (line %d)\n", szat, tag.Name(), tag.m_ncurrent_line-1);
		}

		// give the parameter container a chance to do additional processing
		pl.GetContainer()->SetParameter(*pp);

		return true;
	}
	else
	{
		// if we get here, the parameter is not found.
		// See if the parameter container has defined a property of this name
		int n = pc->FindPropertyIndex(tag.Name());
		if (n >= 0)
		{
			// get the property
			FECoreBase* pp = pc->GetProperty(n);
			if (pp)
			{
				if (!tag.isleaf())
				{
					++tag;
					do
					{
						if (ReadParameter(tag, pp) == false) throw XMLReader::InvalidTag(tag);
						++tag;
					}
					while (!tag.isend());
				}
				else
				{
					// there should be one parameter with the same name as the tag
					if (ReadParameter(tag, pp) == false) throw XMLReader::InvalidTag(tag);
				}
			}
			else throw XMLReader::InvalidTag(tag);
			return true;
		}
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

//-----------------------------------------------------------------------------
void FEFEBioImport::AddPlotVariable(const char* szvar, vector<int>& item)
{
	PlotVariable var(szvar, item);
	m_plot.push_back(var);
}

//-----------------------------------------------------------------------------
void FEFEBioImport::SetPlotCompression(int n)
{
	m_nplot_compression = n;
}
