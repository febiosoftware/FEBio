// FEBioImport.cpp: implementation of the FEBioImport class.
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
#include "FEBioMeshDataSection.h"
#include "FEBioCodeSection.h"
#include "FEBioRigidSection.h"
#include "FECore/DataStore.h"
#include "FECore/Image.h"
#include "FECore/FEModel.h"
#include "FECore/FECoreKernel.h"
#include <FECore/FESurfaceMap.h>
#include <FECore/FEFunction1D.h>
#include <FECore/tens3d.h>
#include "FECore/DOFS.h"
#include <string.h>
#include <stdarg.h>

//-----------------------------------------------------------------------------
void FEBioImport::Exception::SetErrorString(const char* sz, ...)
{
	// get a pointer to the argument list
	va_list	args;

	// make the message
	va_start(args, sz);
	vsprintf(m_szerr, sz, args);
	va_end(args);
}

//-----------------------------------------------------------------------------
FEBioImport::InvalidVersion::InvalidVersion()
{
	SetErrorString("Invalid version");
}

//-----------------------------------------------------------------------------
FEBioImport::InvalidMaterial::InvalidMaterial(int nel)
{
	SetErrorString("Element %d has an invalid material type", nel);
}

//-----------------------------------------------------------------------------
FEBioImport::InvalidDomainType::InvalidDomainType()	
{
	SetErrorString("Invalid domain type");
}

//-----------------------------------------------------------------------------
FEBioImport::InvalidDomainMaterial::InvalidDomainMaterial()
{
	SetErrorString("Invalid domain material");
}

//-----------------------------------------------------------------------------
FEBioImport::FailedCreatingDomain::FailedCreatingDomain()
{
	SetErrorString("Failed creating domain");
}

//-----------------------------------------------------------------------------
FEBioImport::InvalidElementType::InvalidElementType()
{
	SetErrorString("Invalid element type\n");
}

//-----------------------------------------------------------------------------
FEBioImport::FailedLoadingPlugin::FailedLoadingPlugin(const char* szfile)
{
	SetErrorString("failed loading plugin %s\n", szfile);
}

//-----------------------------------------------------------------------------
FEBioImport::DuplicateMaterialSection::DuplicateMaterialSection()
{
	SetErrorString("Material section has already been defined");
}

//-----------------------------------------------------------------------------
FEBioImport::MissingMaterialProperty::MissingMaterialProperty(const char* szmat, const char* szprop)
{
	SetErrorString("Material \"%s\" needs to have property \"%s\" defined", szmat, szprop);
}

//-----------------------------------------------------------------------------
FEBioImport::FailedAllocatingSolver::FailedAllocatingSolver(const char* sztype)
{
	SetErrorString("Failed allocating solver \"%s\"", sztype);
}

//-----------------------------------------------------------------------------
FEBioImport::InvalidNodeID::InvalidNodeID()
{
	SetErrorString("Invalid node ID");
}

//-----------------------------------------------------------------------------
FEBioImport::MissingSlaveSurface::MissingSlaveSurface()
{
	SetErrorString("Missing contact slave surface");
}

//-----------------------------------------------------------------------------
FEBioImport::MissingMasterSurface::MissingMasterSurface()
{
	SetErrorString("Missing contact master surface");
}

//-----------------------------------------------------------------------------
FEBioImport::DataGeneratorError::DataGeneratorError()
{
	SetErrorString("Error in data generation");
}

//-----------------------------------------------------------------------------
FEBioImport::FailedBuildingPart::FailedBuildingPart(const std::string& partName)
{
	SetErrorString("Failed building part %s", partName.c_str());
}

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
FEBioImport::PlotVariable::PlotVariable(const FEBioImport::PlotVariable& pv)
{
	strcpy(m_szvar, pv.m_szvar);
    strcpy(m_szdom, pv.m_szdom);
	m_item = pv.m_item;
}

FEBioImport::PlotVariable::PlotVariable(const char* szvar, vector<int>& item, const char* szdom)
{
    strcpy(m_szvar, szvar);
    m_item = item;
    strcpy(m_szdom, szdom);
}

//=============================================================================
//  The FEBioImport class imports an XML formatted FEBio input file.
//  The actual file is parsed using the XMLReader class.
//
//-----------------------------------------------------------------------------
void FEBioImport::ClearParams()
{
	m_Param.clear();
}

//-----------------------------------------------------------------------------
FEBioImport::XMLParam* FEBioImport::FindParameter(const char* sz)
{
	for (size_t i=0; i<m_Param.size(); ++i)
	{
		XMLParam& p = m_Param[i];
		if (strcmp(p.m_szname, sz) == 0) return &p;
	}
	return 0;
}

//-----------------------------------------------------------------------------
void FEBioImport::AddParameter(const char* szname, const char* szval)
{
	XMLParam p;
	strcpy(p.m_szname, szname);
	strcpy(p.m_szval , szval);
	m_Param.push_back(p);
}

//-----------------------------------------------------------------------------
const char* FEBioImport::get_value_string(XMLTag& tag)
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
void FEBioImport::value(XMLTag& tag, int& n)
{
	const char* sz = get_value_string(tag);
	n = atoi(sz);
}

//-----------------------------------------------------------------------------
void FEBioImport::value(XMLTag& tag, double& g)
{
	const char* sz = get_value_string(tag);
	g = atof(sz);
}

//-----------------------------------------------------------------------------
void FEBioImport::value(XMLTag& tag, bool& b)
{
	const char* sz = get_value_string(tag);
	int n=0; 
	sscanf(sz, "%d", &n); 
	b = (n != 0); 
}

//-----------------------------------------------------------------------------
void FEBioImport::value(XMLTag& tag, vec3d& v)
{
	const char* sz = get_value_string(tag);
	int n = sscanf(sz, "%lg,%lg,%lg", &v.x, &v.y, &v.z);
	if (n != 3) throw XMLReader::XMLSyntaxError();
}

//-----------------------------------------------------------------------------
void FEBioImport::value(XMLTag& tag, mat3d& m)
{
	const char* sz = get_value_string(tag);
	double xx, xy, xz, yx, yy, yz, zx, zy, zz;
	int n = sscanf(sz, "%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg", &xx, &xy, &xz, &yx, &yy, &yz, &zx, &zy, &zz);
	if (n != 9) throw XMLReader::XMLSyntaxError();
	m = mat3d(xx, xy, xz, yx, yy, yz, zx, zy, zz);
}

//-----------------------------------------------------------------------------
void FEBioImport::value(XMLTag& tag, mat3ds& m)
{
	const char* sz = get_value_string(tag);
	double x, y, z, xy, yz, xz;
	int n = sscanf(sz, "%lg,%lg,%lg,%lg,%lg,%lg", &x, &y, &z, &xy, &yz, &xz);
	if (n != 6) throw XMLReader::XMLSyntaxError();
	m = mat3ds(x, y, z, xy, yz, xz);
}

//-----------------------------------------------------------------------------
void FEBioImport::value(XMLTag& tag, tens3drs& m)
{
	double v[18];
	int n = value(tag, v, 18);
	if (n != 18) throw XMLReader::InvalidValue(tag);
	for (int i=0; i<18; ++i) m.d[i] = v[i];
}

//-----------------------------------------------------------------------------
void FEBioImport::value(XMLTag& tag, char* szstr)
{
	const char* sz = get_value_string(tag);
	strcpy(szstr, sz); 
}

//-----------------------------------------------------------------------------
int FEBioImport::value(XMLTag& tag, int* pi, int n)
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
int FEBioImport::value(XMLTag& tag, double* pf, int n)
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
FEBioImport::FEBioImport()
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
	m_map["Contact"    ] = new FEBioContactSection    (this);
	m_map["Discrete"   ] = new FEBioDiscreteSection   (this);
	m_map["Code"       ] = new FEBioCodeSection       (this);	// added in FEBio 2.4 (experimental feature!)

	// version 2.5 only
	m_map["MeshData"] = new FEBioMeshDataSection(this);
	m_map["Rigid"   ] = new FEBioRigidSection   (this); // added in FEBio 2.6 (experimental feature!)
}

//-----------------------------------------------------------------------------
FEBioImport::~FEBioImport()
{
}

//-----------------------------------------------------------------------------
bool FEBioImport::Load(FEModel& fem, const char* szfile)
{
	// keep a pointer to the fem object
	m_pfem = &fem;

	// keep a pointer to the mesh
	m_pMesh = &fem.GetMesh();

	// intialize some variables
	m_pStep = 0;	// zero step pointer
	m_nsteps = 0;	// reset step section counter
	m_szmod[0] = 0;
	m_nversion = -1;
	m_szdmp[0] = 0;
	m_szlog[0] = 0;
	m_szplt[0] = 0;

	// plot output
	m_szplot_type[0] = 0;
	m_plot.clear();
	m_nplot_compression = 0;

	m_data.clear();

	// default element type
	m_ntet4  = FE_TET4G1;
	m_nhex8  = FE_HEX8G8;
	m_ntet10 = FE_TET10G8;
	m_ntet15 = FE_TET15G15;
	m_ntet20 = FE_TET20G15;
	m_ntri6  = FE_TRI6G7;
	m_ntri3  = FE_TRI3G3;
	m_ntri7  = FE_TRI7G7;
	m_ntri10 = FE_TRI10G7;

	// 3-field formulation flags
	m_b3field_hex = true;
	m_b3field_tet = false;

	// UT4 formulation off by default
	m_but4 = false;

	// extract the path
	strcpy(m_szpath, szfile);
	char* ch = strrchr(m_szpath, '\\');
	if (ch==0) ch = strrchr(m_szpath, '/');
	if (ch==0) m_szpath[0] = 0; else *(ch+1)=0;

	// clean up
	ClearParams();
	ClearDataArrays();

	// read the file
	return ReadFile(szfile);
}

//-----------------------------------------------------------------------------
// This function parses the XML input file. The broot parameter is used to indicate
// if this is the master file or an included file. 
bool FEBioImport::ReadFile(const char* szfile, bool broot)
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
		return errf("An error occured while finding the febio_spec tag.\nIs this a valid FEBio input file?\n\n");
	}

	// parse the file
	try
	{
		// get the version number
		ParseVersion(tag);

		// FEBio2 only supports file version 1.2, 2.0, and 2.5
		if ((m_nversion != 0x0102) && 
			(m_nversion != 0x0200) && 
			(m_nversion != 0x0205)) throw InvalidVersion();

		// For versions before 2.5 we need to allocate all the degrees of freedom beforehand. 
		// This is necessary because the Module section doesn't have to defined until a Control section appears.
		// That means that model components that depend on DOFs can be defined before the Module tag (e.g. in multi-step analyses) and this leads to problems.
		// In 2.5 this is solved by requiring that the Module tag is defined at the top of the file. 
		if (broot && (m_nversion < 0x0205))
		{
			// We need to define a default Module type since before 2.5 this tag is optional for structural mechanics model definitions.
			sprintf(m_szmod, "solid");

			// Reset degrees of
			FEModel& fem = *GetFEModel();
			DOFS& dofs = fem.GetDOFS();
			dofs.Reset();

			// Add the default variables and degrees of freedom
			int varD = dofs.AddVariable("displacement", VAR_VEC3);
			dofs.SetDOFName(varD, 0, "x");
			dofs.SetDOFName(varD, 1, "y");
			dofs.SetDOFName(varD, 2, "z");
			int varQ = dofs.AddVariable("rotation", VAR_VEC3);
			dofs.SetDOFName(varQ, 0, "u");
			dofs.SetDOFName(varQ, 1, "v");
			dofs.SetDOFName(varQ, 2, "w");
			int varP = dofs.AddVariable("fluid pressure");
			dofs.SetDOFName(varP, 0, "p");
            int varSP = dofs.AddVariable("shell fluid pressure");
            dofs.SetDOFName(varSP, 0, "q");
			int varQR = dofs.AddVariable("rigid rotation", VAR_VEC3);
			dofs.SetDOFName(varQR, 0, "Ru");
			dofs.SetDOFName(varQR, 1, "Rv");
			dofs.SetDOFName(varQR, 2, "Rw");
			int varT = dofs.AddVariable("temperature");
			dofs.SetDOFName(varT, 0, "T");
			int varV = dofs.AddVariable("velocity", VAR_VEC3);
			dofs.SetDOFName(varV, 0, "vx");
			dofs.SetDOFName(varV, 1, "vy");
			dofs.SetDOFName(varV, 2, "vz");
			int varE = dofs.AddVariable("fluid dilation");
			dofs.SetDOFName(varE, 0, "e");
			int varQP = dofs.AddVariable("previous rotation", VAR_VEC3);
			dofs.SetDOFName(varQP, 0, "up");
			dofs.SetDOFName(varQP, 1, "vp");
			dofs.SetDOFName(varQP, 2, "wp");
			int varQV = dofs.AddVariable("shell velocity", VAR_VEC3);
			dofs.SetDOFName(varQV, 0, "vu");
			dofs.SetDOFName(varQV, 1, "vv");
			dofs.SetDOFName(varQV, 2, "vw");
			int varQA = dofs.AddVariable("shell acceleration", VAR_VEC3);
			dofs.SetDOFName(varQA, 0, "au");
			dofs.SetDOFName(varQA, 1, "av");
			dofs.SetDOFName(varQA, 2, "aw");
			int varQVP = dofs.AddVariable("previous shell velocity", VAR_VEC3);
			dofs.SetDOFName(varQVP, 0, "vup");
			dofs.SetDOFName(varQVP, 1, "vvp");
			dofs.SetDOFName(varQVP, 2, "vwp");
			int varQAP = dofs.AddVariable("previous shell acceleration", VAR_VEC3);
			dofs.SetDOFName(varQAP, 0, "aup");
			dofs.SetDOFName(varQAP, 1, "avp");
			dofs.SetDOFName(varQAP, 2, "awp");
			// must be last variable definition!!
			int varC = dofs.AddVariable("concentration", VAR_ARRAY); // we start with zero concentrations
            // must be last variable definition!!
            int varSC = dofs.AddVariable("shell concentration", VAR_ARRAY); // we start with zero concentrations
		}

		// parse the file
		++tag;

		// From version 2.5 and up the first tag of the master file has to be the Module tag.
		if (broot && (m_nversion >= 0x0205))
		{
			if (tag != "Module")
			{
				return errf("First tag must be the Module section.\n\n");
			}

			// try to find a section parser
			FEBioFileSectionMap::iterator is = m_map.find(tag.Name());

			// make sure we found a section reader
			if (is == m_map.end()) throw XMLReader::InvalidTag(tag);

			// parse the module tag
			is->second->Parse(tag);

			// Now that the Module tag is read in, we'll want to create an analysis step.
			// Creating an analysis step will allocate a solver class (based on the module) 
			// and this in turn will allocate the degrees of freedom.
			// TODO: This is kind of a round-about way and I really want to find a better solution.
			GetStep();

			// let's get the next tag
			++tag;
		}

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
					char sz[512] = {0};
					sprintf(sz, "febio_spec/%s", tag.Name());
					if (xml2.FindTag(sz, tag2) == false) return errf("FATAL ERROR: Couldn't find %s section in file %s.\n\n", tag.Name(), szinc);

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
	catch (XMLReader::Error& e)
	{
		return errf("FATAL ERROR: %s (line %d)\n", e.GetErrorString(), xml.GetCurrentLine());
	}
	// --- FEBioImport Exceptions ---
	catch (FEBioImport::Exception& e)
	{
		return errf("FATAL ERROR: %s (line %d)\n", e.GetErrorString(), xml.GetCurrentLine());
	}
	// --- Exception from DataStore ---
	catch (UnknownDataField& e)
	{
		return errf("Fatal Error: \"%s\" is not a valid field variable name (line %d)\n", e.m_szdata, xml.GetCurrentLine()-1);
	}
	// --- Unknown exceptions ---
	catch (...)
	{
		return errf("FATAL ERROR: unrecoverable error (line %d)\n", xml.GetCurrentLine());
		return false;
	}

	// close the XML file
	xml.Close();

	// we're done!
	return true;
}

//-----------------------------------------------------------------------------
FEAnalysis* FEBioImport::GetStep()
{
	if (m_pStep == 0)
	{
		m_pStep = CreateNewStep();
		m_pfem->AddStep(m_pStep);	
		if (m_pfem->Steps() == 1) 
		{
			m_pfem->SetCurrentStep(m_pStep);
			m_pfem->SetCurrentStepIndex(0);
		}
	}
	return m_pStep;
}

//-----------------------------------------------------------------------------
FESolver* FEBioImport::BuildSolver(const char* sztype, FEModel& fem)
{
	FESolver* ps = fecore_new<FESolver>(FESOLVER_ID, sztype, &fem);
	return ps;
}

//-----------------------------------------------------------------------------
FEAnalysis* FEBioImport::CreateNewStep()
{
	FEAnalysis* pstep = new FEAnalysis(m_pfem);

	// make sure we have a solver defined
	FESolver* psolver = pstep->GetFESolver();
	if (psolver == 0)
	{
		psolver = BuildSolver(m_szmod, *GetFEModel());
		if (psolver == 0) throw FEBioImport::FailedAllocatingSolver(m_szmod);
		pstep->SetFESolver(psolver);
	}
	return pstep;
}

//-----------------------------------------------------------------------------
//! This function parses the febio_spec tag for the version number
void FEBioImport::ParseVersion(XMLTag &tag)
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
bool FEBioImport::ReadParameter(XMLTag& tag, FEParameterList& pl, const char* szparam)
{
	// see if we can find this parameter
	FEParam* pp = pl.Find((szparam == 0 ? tag.Name() : szparam));
	if (pp == 0) return false;
	
	if (pp->dim() == 1)
	{
		switch (pp->type())
		{
		case FE_PARAM_DOUBLE  : value(tag, pp->value<double  >()); break;
		case FE_PARAM_INT     : value(tag, pp->value<int     >()); break;
		case FE_PARAM_BOOL    : value(tag, pp->value<bool    >()); break;
		case FE_PARAM_VEC3D   : value(tag, pp->value<vec3d   >()); break;
		case FE_PARAM_MAT3D   : value(tag, pp->value<mat3d   >()); break;
		case FE_PARAM_MAT3DS  : value(tag, pp->value<mat3ds  >()); break;
		case FE_PARAM_TENS3DRS: value(tag, pp->value<tens3drs>()); break;
		case FE_PARAM_STRING  : value(tag, pp->cvalue()); break;
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
		case FE_PARAM_DATA_ARRAY:
			{
				// get the surface map
				FEDataArray& map = pp->value<FEDataArray>();

				// Make sure that the tag is a leaf
				if (!tag.isleaf()) throw XMLReader::InvalidValue(tag);

				// read the surface map data
				const char* szmap = tag.AttributeValue("surface_data", true);
				if (szmap)
				{
					FEDataArray* pdata = FindDataArray(szmap);
					if (pdata == 0) throw XMLReader::InvalidAttributeValue(tag, "surface_data");

					// make sure the types match
					if (map.DataType() != pdata->DataType()) throw XMLReader::InvalidAttributeValue(tag, "surface_data", szmap);

					// copy data
					map = *pdata;
				}
				else 
				{
					const char* szmap = tag.AttributeValue("node_data", true);
					if (szmap)
					{
						FEDataArray* pdata = FindDataArray(szmap);
						if (pdata == 0) throw XMLReader::InvalidAttributeValue(tag, "node_data");

						// make sure the types match
						if (map.DataType() != pdata->DataType()) throw XMLReader::InvalidAttributeValue(tag, "node_data", szmap);

						// copy data
						map = *pdata;
					}
					else 
					{
						if (map.DataType() == FE_DOUBLE)
						{
							double v;
							tag.value(v);
							map.set(v);
						}
						else if (map.DataType() == FE_VEC2D)
						{
							double v[2] = {0};
							tag.value(v, 2);
							map.set(vec2d(v[0], v[1]));
						}
						else if (map.DataType() == FE_VEC3D)
						{
							double v[3] = {0};
							tag.value(v, 3);
							map.set(vec3d(v[0], v[1], v[2]));
						}
					}
				}
			};
			break;
		case FE_PARAM_FUNC1D:
			{
				int lc = -1;
				tag.AttributeValue("lc", lc, true);
				double v = 1.0;
				tag.value(v);

				FEFunction1D& f = pp->value<FEFunction1D>();
				f.SetLoadCurveIndex(lc - 1, v);
			}
			break;
		default:
			assert(false);
			return false;
		}
	}
	else
	{
		switch (pp->type())
		{
		case FE_PARAM_INT   : value(tag, pp->pvalue<int   >(), pp->dim()); break;
		case FE_PARAM_DOUBLE: value(tag, pp->pvalue<double>(), pp->dim()); break;
		}
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
				switch (pp->type())
				{
				case FE_PARAM_INT   : pp->SetLoadCurve(lc); break;
				case FE_PARAM_BOOL  : pp->SetLoadCurve(lc); break;
				case FE_PARAM_DOUBLE: pp->SetLoadCurve(lc, pp->value<double>()); break;
				case FE_PARAM_VEC3D : pp->SetLoadCurve(lc, pp->value<vec3d >()); break;
				case FE_PARAM_FUNC1D: break; // don't do anything for 1D functions since the lc attribute is already processed.
				default:
					assert(false);
				}
			}
/*			else 
			{
				throw XMLReader::InvalidAttributeValue(tag, szat, tag.m_att[i].m_szatv);
			}
*/		}
		// This is not true. Parameters can have attributes that are used for other purposed. E.g. The local fiber option.
//		else felog.printf("WARNING: attribute \"%s\" of parameter \"%s\" ignored (line %d)\n", szat, tag.Name(), tag.m_ncurrent_line-1);
	}

	// give the parameter container a chance to do additional processing
	pl.GetContainer()->SetParameter(*pp);

	return true;
}

//-----------------------------------------------------------------------------
//! This function parses a parameter list
bool FEBioImport::ReadParameter(XMLTag& tag, FECoreBase* pc, const char* szparam)
{
	// get the parameter list
	FEParameterList& pl = pc->GetParameterList();

	// see if we can find this parameter
	if (ReadParameter(tag, pl, szparam) == false)
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
				// process attributes first
				for (int i = 0; i<tag.m_natt; ++i)
				{
					XMLAtt& att = tag.m_att[i];
					if (pp->SetAttribute(att.m_szatt, att.m_szatv) == false) throw XMLReader::InvalidAttributeValue(tag, att.m_szatt);
				}

				// process the parameter lists
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
		else return false;
	}
	return true;
}

//-----------------------------------------------------------------------------
void FEBioImport::ReadParameterList(XMLTag& tag, FEParameterList& pl)
{
	// Make sure there is something to read
	if (tag.isleaf() || tag.isempty()) return;

	// parse the child tags
	++tag;
	do
	{
		if (ReadParameter(tag, pl) == false) throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());
}

//-----------------------------------------------------------------------------
void FEBioImport::ReadList(XMLTag& tag, vector<int>& l)
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
void FEBioImport::SetDumpfileName(const char* sz) { sprintf(m_szdmp, sz); }
void FEBioImport::SetLogfileName (const char* sz) { sprintf(m_szlog, sz); }
void FEBioImport::SetPlotfileName(const char* sz) { sprintf(m_szplt, sz); }

//-----------------------------------------------------------------------------
void FEBioImport::AddDataRecord(DataRecord* pd)
{
	m_data.push_back(pd);
}

//-----------------------------------------------------------------------------
void FEBioImport::AddPlotVariable(const char* szvar, vector<int>& item, const char* szdom)
{
    PlotVariable var(szvar, item, szdom);
    m_plot.push_back(var);
}

//-----------------------------------------------------------------------------
void FEBioImport::SetPlotCompression(int n)
{
	m_nplot_compression = n;
}

//-----------------------------------------------------------------------------
// This tag parses a node set.
FENodeSet* FEBioImport::ParseNodeSet(XMLTag& tag, const char* szatt)
{
	FEMesh& mesh = GetFEModel()->GetMesh();
	FENodeSet* pns = 0;

	// see if the set attribute is defined
	const char* szset = tag.AttributeValue(szatt, true);
	if (szset)
	{
		// Make sure this is an empty tag
		if (tag.isempty() == false) throw XMLReader::InvalidValue(tag);

		// find the node set
		pns = mesh.FindNodeSet(szset);
		if (pns == 0) throw XMLReader::InvalidAttributeValue(tag, szatt, szset);
	}
	else
	{
		// This defines a node set, so we need a name tag
		// (For now this name is optional)
		const char* szname = tag.AttributeValue("name", true);
		if (szname == 0) szname = "_unnamed";

		// create a new node set
		pns = new FENodeSet(&mesh);
		pns->SetName(szname);

		// add the nodeset to the mesh
		mesh.AddNodeSet(pns);

		// read the nodes
		if (tag.isleaf())
		{
			// This format is deprecated
			vector<int> l;
			ReadList(tag, l);
			for (int i=0; i<l.size(); ++i) pns->add(FindNodeFromID(l[i]));
		}
		else
		{
			// read the nodes
			++tag;
			do
			{
				if (tag == "node")
				{
					int nid = -1;
					tag.AttributeValue("id", nid);

					nid = FindNodeFromID(nid);
					pns->add(nid);
				}
				else if (tag == "NodeSet")
				{
					const char* szset = tag.AttributeValue(szatt);
					
					// Make sure this is an empty tag
					if (tag.isempty() == false) throw XMLReader::InvalidValue(tag);

					// find the node set
					FENodeSet* ps = mesh.FindNodeSet(szset);
					if (ps == 0) throw XMLReader::InvalidAttributeValue(tag, szatt, szset);

					// add the node set
					pns->add(*ps);
				}
				else if (tag == "node_list")
				{
					vector<int> nl;
					ReadList(tag, nl);
					for (int i=0; i<nl.size(); ++i) pns->add(FindNodeFromID(nl[i]));
				}
				else throw XMLReader::InvalidTag(tag);
				++tag;
			}
			while (!tag.isend());
		}
	}

	return pns;
}

//-----------------------------------------------------------------------------
FESurface* FEBioImport::ParseSurface(XMLTag& tag, const char* szatt)
{
	FEMesh& m = GetFEModel()->GetMesh();

	// create new surface
	FESurface* psurf = new FESurface(&m);

	// see if the surface is referenced by a set of defined explicitly
	const char* szset = tag.AttributeValue(szatt, true);
	if (szset)
	{
		// make sure this tag does not have any children
		if (!tag.isleaf()) throw XMLReader::InvalidTag(tag);

		// see if we can find the facet set
		FEMesh& m = GetFEModel()->GetMesh();
		FEFacetSet* ps = 0;
		for (int i=0; i<m.FacetSets(); ++i)
		{
			FEFacetSet& fi = m.FacetSet(i);
			if (strcmp(fi.GetName(), szset) == 0)
			{
				ps = &fi;
				break;
			}
		}

		// create a surface from the facet set
		if (ps)
		{
			if (BuildSurface(*psurf, *ps) == false) throw XMLReader::InvalidTag(tag);
		}
		else throw XMLReader::InvalidAttributeValue(tag, "set", szset);
	}
	else
	{
		// count how many pressure cards there are
		int npr = tag.children();
		psurf->Create(npr);

		++tag;
		int nf[FEElement::MAX_NODES ], N;
		for (int i=0; i<npr; ++i)
		{
			FESurfaceElement& el = psurf->Element(i);

			if      (tag == "quad4") el.SetType(FE_QUAD4G4);
			else if (tag == "tri3" ) el.SetType(m_ntri3);
			else if (tag == "tri6" ) el.SetType(m_ntri6);
			else if (tag == "tri7" ) el.SetType(m_ntri7);
			else if (tag == "tri10") el.SetType(m_ntri10);
			else if (tag == "quad8") el.SetType(FE_QUAD8G9);
			else if (tag == "quad9") el.SetType(FE_QUAD9G9);
			else throw XMLReader::InvalidTag(tag);

			N = el.Nodes();
			tag.value(nf, N);
			for (int j=0; j<N; ++j) el.m_node[j] = nf[j]-1;

			++tag;
		}
	}

	return psurf;
}

//-----------------------------------------------------------------------------
bool FEBioImport::BuildSurface(FESurface& s, FEFacetSet& fs)
{
	FEModel& fem = *GetFEModel();
	FEMesh& m = fem.GetMesh();
	int NN = m.Nodes();

	// count nr of faces
	int faces = fs.Faces();

	// allocate storage for faces
	s.Create(faces);

	// read faces
	for (int i=0; i<faces; ++i)
	{
		FESurfaceElement& el = s.Element(i);
		FEFacetSet::FACET& fi = fs.Face(i);

		if      (fi.ntype == 4) el.SetType(FE_QUAD4G4);
		else if (fi.ntype == 3) el.SetType(m_ntri3);
		else if (fi.ntype == 6) el.SetType(m_ntri6);
		else if (fi.ntype == 7) el.SetType(m_ntri7);
		else if (fi.ntype == 10) el.SetType(m_ntri10);
		else if (fi.ntype == 8) el.SetType(FE_QUAD8G9);
		else if (fi.ntype == 9) el.SetType(FE_QUAD9G9);
		else return false;

		int N = el.Nodes(); assert(N == fi.ntype);
		for (int j=0; j<N; ++j) el.m_node[j] = fi.node[j];
	}

	// copy the name
	s.SetName(fs.GetName());

	return true;
}


//-----------------------------------------------------------------------------
void FEBioImport::ClearDataArrays()
{
	// clear the surface maps
	for (int i=0; i<(int) m_Data.size(); ++i) delete m_Data[i].second;
	m_Data.clear();
}

//-----------------------------------------------------------------------------
void FEBioImport::AddDataArray(const char* szname, FEDataArray* map)
{
	m_Data.push_back(pair<string,FEDataArray*>(string(szname), map));
}

//-----------------------------------------------------------------------------
FEDataArray* FEBioImport::FindDataArray(const char* szmap)
{
	string name(szmap);
	for (int i=0; i<(int) m_Data.size(); ++i)
	{
		if (m_Data[i].first == name) return m_Data[i].second;
	}
	return 0;
}

//-----------------------------------------------------------------------------
void FEBioImport::ParseDataArray(XMLTag& tag, FEDataArray& map, const char* sztag)
{
	int dataType = map.DataType();

	if (dataType == FE_DOUBLE)
	{
		++tag;
		do
		{
			if (tag == sztag)
			{
				int nid;
				tag.AttributeValue("lid", nid);

				double v;
				tag.value(v);

				map.set(nid - 1, v);
			}
			else throw XMLReader::InvalidTag(tag);
			++tag;
		}
		while (!tag.isend());
	}
	else if (dataType == FE_VEC3D)
	{
		++tag;
		do
		{
			if (tag == sztag)
			{
				int nid;
				tag.AttributeValue("lid", nid);

				double v[3];
				tag.value(v, 3);

				map.set(nid - 1, vec3d(v[0], v[1], v[2]));
			}
			else throw XMLReader::InvalidTag(tag);
			++tag;
		}
		while (!tag.isend());
	}
}

//-----------------------------------------------------------------------------
void FEBioImport::BuildNodeList()
{
	// find the min, max ID
	// (We assume that they are given by the first and last node)
	FEMesh& mesh = GetFEModel()->GetMesh();
	int NN = mesh.Nodes();
	int nmin = mesh.Node(0   ).GetID();
	int nmax = mesh.Node(NN-1).GetID();
	assert(nmax >= nmin);

	// get the range
	int nn = nmax - nmin+1;

	// allocate list
	m_node_off = nmin;
	m_node_list.assign(nn, -1);

	// build the list
	for (int i=0; i<NN; ++i)
	{
		int nid = mesh.Node(i).GetID();
		m_node_list[nid - m_node_off] = i;
	}
}

//-----------------------------------------------------------------------------
int FEBioImport::FindNodeFromID(int nid)
{
	int N = (int) m_node_list.size();
	if (N > 0)
	{
		int n = nid - m_node_off;
		if ((n>=0)&&(n<N)) return m_node_list[n];
	}
	throw FEBioImport::InvalidNodeID();
	return -1;
}

//-----------------------------------------------------------------------------
void FEBioImport::GlobalToLocalID(int* l, int n, vector<int>& m)
{
	assert((int) m.size()==n);
	for (int i=0; i<n; ++i)
	{
		m[i] = FindNodeFromID(l[i]);
	}
}

//-----------------------------------------------------------------------------
int FEBioImport::ReadNodeID(XMLTag& tag)
{
	int n = atoi(tag.AttributeValue("id"));
	return FindNodeFromID(n);
}

//-----------------------------------------------------------------------------
FEBioImport::SurfacePair* FEBioImport::FindSurfacePair(const char* szname)
{
	for (int i=0; i<m_surfacePair.size(); ++i)
		if (strcmp(m_surfacePair[i].szname, szname) == 0) return &m_surfacePair[i];
	return 0;
}

//-----------------------------------------------------------------------------
FEBioImport::NodeSetPair* FEBioImport::FindNodeSetPair(const char* szname)
{
	for (int i = 0; i<m_nsetPair.size(); ++i)
	if (strcmp(m_nsetPair[i].szname, szname) == 0) return &m_nsetPair[i];
	return 0;
}

//-----------------------------------------------------------------------------
FEBioImport::NodeSetSet* FEBioImport::FindNodeSetSet(const char* szname)
{
	for (int i = 0; i<m_nsetSet.size(); ++i)
	if (strcmp(m_nsetSet[i].szname, szname) == 0) return &m_nsetSet[i];
	return 0;
}
