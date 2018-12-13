// FileImport.cpp: implementation of the FEFileImport class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FileImport.h"
#include <FECore/Image.h>
#include <FECore/FENodeDataMap.h>
#include <FECore/FESurfaceMap.h>
#include <FECore/FEFunction1D.h>
#include <FECore/FEModel.h>
#include <FECore/FEMaterial.h>
#include <FECore/FEModelParam.h>
#include <FECore/FESurface.h>
#include <FECore/FESurfaceLoad.h>
#include <FECore/FEBodyLoad.h>
#include <FECore/FEDomainMap.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include "FEBioImport.h"

//-----------------------------------------------------------------------------
FEFileException::FEFileException()
{
	m_szerr[0] = 0;
}

//-----------------------------------------------------------------------------
FEFileException::FEFileException(const char* sz, ...)
{
	// get a pointer to the argument list
	va_list	args;

	// make the message
	va_start(args, sz);
	vsprintf(m_szerr, sz, args);
	va_end(args);
}

//-----------------------------------------------------------------------------
void FEFileException::SetErrorString(const char* sz, ...)
{
	// get a pointer to the argument list
	va_list	args;

	// make the message
	va_start(args, sz);
	vsprintf(m_szerr, sz, args);
	va_end(args);
}

//-----------------------------------------------------------------------------
FEModel* FEFileSection::GetFEModel() { return m_pim->GetFEModel(); }

//-----------------------------------------------------------------------------
FEModelBuilder* FEFileSection::GetBuilder() { return m_pim->GetBuilder(); }

//-----------------------------------------------------------------------------
const char* FEFileSection::get_value_string(XMLTag& tag)
{
	const char* sz = tag.szvalue();
	if (sz[0] == '@')
	{
		FEFileParam* p = m_pim->FindFileParameter(sz + 1);
		if (p == 0) throw XMLReader::InvalidValue(tag);
		sz = p->m_szval;
	}
	return sz;
}


//-----------------------------------------------------------------------------
void FEFileSection::value(XMLTag& tag, int& n)
{
	const char* sz = get_value_string(tag);
	n = atoi(sz);
}

//-----------------------------------------------------------------------------
void FEFileSection::value(XMLTag& tag, double& g)
{
	const char* sz = get_value_string(tag);
	g = atof(sz);
}

//-----------------------------------------------------------------------------
void FEFileSection::value(XMLTag& tag, bool& b)
{
	const char* sz = get_value_string(tag);
	int n = 0;
	sscanf(sz, "%d", &n);
	b = (n != 0);
}

//-----------------------------------------------------------------------------
void FEFileSection::value(XMLTag& tag, vec3d& v)
{
	const char* sz = get_value_string(tag);
	int n = sscanf(sz, "%lg,%lg,%lg", &v.x, &v.y, &v.z);
	if (n != 3) throw XMLReader::XMLSyntaxError();
}

//-----------------------------------------------------------------------------
void FEFileSection::value(XMLTag& tag, mat3d& m)
{
	const char* sz = get_value_string(tag);
	double xx, xy, xz, yx, yy, yz, zx, zy, zz;
	int n = sscanf(sz, "%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg", &xx, &xy, &xz, &yx, &yy, &yz, &zx, &zy, &zz);
	if (n != 9) throw XMLReader::XMLSyntaxError();
	m = mat3d(xx, xy, xz, yx, yy, yz, zx, zy, zz);
}

//-----------------------------------------------------------------------------
void FEFileSection::value(XMLTag& tag, mat3ds& m)
{
	const char* sz = get_value_string(tag);
	double x, y, z, xy, yz, xz;
	int n = sscanf(sz, "%lg,%lg,%lg,%lg,%lg,%lg", &x, &y, &z, &xy, &yz, &xz);
	if (n != 6) throw XMLReader::XMLSyntaxError();
	m = mat3ds(x, y, z, xy, yz, xz);
}

//-----------------------------------------------------------------------------
void FEFileSection::value(XMLTag& tag, tens3drs& m)
{
	double v[18];
	int n = value(tag, v, 18);
	if (n != 18) throw XMLReader::InvalidValue(tag);
	for (int i = 0; i<18; ++i) m.d[i] = v[i];
}

//-----------------------------------------------------------------------------
void FEFileSection::value(XMLTag& tag, char* szstr)
{
	const char* sz = get_value_string(tag);
	strcpy(szstr, sz);
}

//-----------------------------------------------------------------------------
void FEFileSection::value(XMLTag& tag, std::string& v)
{
	v = get_value_string(tag);
}

//-----------------------------------------------------------------------------
int FEFileSection::value(XMLTag& tag, int* pi, int n)
{
	const char* sz = get_value_string(tag);
	int nr = 0;
	for (int i = 0; i<n; ++i)
	{
		const char* sze = strchr(sz, ',');

		pi[i] = atoi(sz);
		nr++;

		if (sze) sz = sze + 1;
		else break;
	}
	return nr;
}

//-----------------------------------------------------------------------------
int FEFileSection::value(XMLTag& tag, double* pf, int n)
{
	const char* sz = get_value_string(tag);
	int nr = 0;
	for (int i = 0; i<n; ++i)
	{
		const char* sze = strchr(sz, ',');

		pf[i] = atof(sz);
		nr++;

		if (sze) sz = sze + 1;
		else break;
	}
	return nr;
}

//-----------------------------------------------------------------------------
int FEFileSection::ReadNodeID(XMLTag& tag)
{
	int n = atoi(tag.AttributeValue("id"));
	int node = GetBuilder()->FindNodeFromID(n);
	if (node == -1) throw FEFileException("Invalid node ID");
	return node;
}

//-----------------------------------------------------------------------------
bool parseEnumParam(FEParam* pp, const char* val)
{
	assert(pp->type() == FE_PARAM_INT);

	const char* ch = pp->enums();
	int n = 0;
	bool bfound = false;
	while (ch && *ch)
	{
		if (strcmp(ch, val) == 0)
		{
			pp->value<int>() = n;
			bfound = true;
			break;
		}
		ch = strchr(ch, '\0');
		if (ch) ch++;
		n++;
	}
	
	return bfound;
}

//-----------------------------------------------------------------------------
//! This function parese a parameter list
bool FEFileSection::ReadParameter(XMLTag& tag, FEParameterList& pl, const char* szparam, FECoreBase* pc, bool parseAttributes)
{
	// see if we can find this parameter
	FEParam* pp = pl.FindFromName((szparam == 0 ? tag.Name() : szparam));
	if (pp == 0) return false;

	if (pp->dim() == 1)
	{
		switch (pp->type())
		{
		case FE_PARAM_DOUBLE: value(tag, pp->value<double  >()); break;
		case FE_PARAM_INT: 
			{
				const char* szenum = pp->enums();
				if (szenum == 0)
				{
					value(tag, pp->value<int>());
				}
				else
				{
					bool bfound = parseEnumParam(pp, get_value_string(tag));
					if (bfound == false) throw XMLReader::InvalidValue(tag);
				}
			}
			break;
		case FE_PARAM_BOOL: value(tag, pp->value<bool    >()); break;
		case FE_PARAM_VEC3D: value(tag, pp->value<vec3d   >()); break;
		case FE_PARAM_MAT3D: value(tag, pp->value<mat3d   >()); break;
		case FE_PARAM_MAT3DS: value(tag, pp->value<mat3ds  >()); break;
		case FE_PARAM_TENS3DRS: value(tag, pp->value<tens3drs>()); break;
		case FE_PARAM_STRING: value(tag, pp->cvalue()); break;
		case FE_PARAM_STD_STRING: value(tag, pp->value<string>()); break;
		case FE_PARAM_IMAGE_3D:
		{
			// get the file name
			const char* szfile = tag.AttributeValue("file");

			++tag;
			int n[3] = { 0 };
			bool bend = false;
			Image::ImageFormat fmt = Image::RAW8;
			do
			{
				if (tag == "size") tag.value(n, 3);
				else if (tag == "format")
				{
					const char* szfmt = tag.szvalue();
					// figure out the image format
					if      (strcmp(szfmt, "RAW8"  ) == 0) fmt = Image::RAW8;
					else if (strcmp(szfmt, "RAW16U") == 0) fmt = Image::RAW16U;
					else throw XMLReader::InvalidValue(tag);
				}
				else if (tag == "endianness") tag.value(bend);
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
			if (ch == 0) ch = strrchr(szin, '/');
			if (ch == 0)
			{
				// pre-pend the name with the input path
				sprintf(szin, "%s%s", m_pim->GetFilePath(), szfile);
			}

			// Try to load the image file
			if (im.Load(szin, fmt, bend) == false) throw XMLReader::InvalidValue(tag);
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
				FESurfaceMap* pmap = dynamic_cast<FESurfaceMap*>(&map);
				if (pmap == 0) throw XMLReader::InvalidTag(tag);

				FESurfaceMap* pdata = dynamic_cast<FESurfaceMap*>(GetFEModel()->FindDataArray(szmap));
				if (pdata == 0) throw XMLReader::InvalidAttributeValue(tag, "surface_data");

				// make sure the types match
				if (map.DataSize() != pdata->DataSize()) throw XMLReader::InvalidAttributeValue(tag, "surface_data", szmap);

				// copy data
				*pmap = *pdata;
			}
			else
			{
				const char* szmap = tag.AttributeValue("node_data", true);
				if (szmap)
				{
					FENodeDataMap* pmap = dynamic_cast<FENodeDataMap*>(&map);
					if (pmap == 0) throw XMLReader::InvalidTag(tag);

					FENodeDataMap* pdata = dynamic_cast<FENodeDataMap*>(GetFEModel()->FindDataArray(szmap));
					if (pdata == 0) throw XMLReader::InvalidAttributeValue(tag, "node_data");

					// make sure the types match
					if (map.DataSize() != pdata->DataSize()) throw XMLReader::InvalidAttributeValue(tag, "node_data", szmap);

					// copy data
					*pmap = *pdata;
				}
				else
				{
					FEDataType dataType = map.DataType();
					if (dataType == FE_DOUBLE)
					{
						double v;
						tag.value(v);
						map.fillValue(v);
					}
					else if (dataType == FE_VEC2D)
					{
						double v[2] = { 0 };
						tag.value(v, 2);
						map.fillValue(vec2d(v[0], v[1]));
					}
					else if (dataType == FE_VEC3D)
					{
						double v[3] = { 0 };
						tag.value(v, 3);
						map.fillValue(vec3d(v[0], v[1], v[2]));
					}
				}
			}
		};
		break;
		case FE_PARAM_STD_VECTOR_VEC2D:
		{
			// make sure this is leaf
			if (tag.isempty()) throw XMLReader::InvalidValue(tag);

			std::vector<vec2d>& data = pp->value< std::vector<vec2d> >();

			// Note that this parameter is read in point per point, not all at once!
			double d[2];
			tag.value(d, 2);
			data.push_back(vec2d(d[0], d[1]));
		}
		break;
		case FE_PARAM_STD_VECTOR_STRING:
		{
			// make sure this is leaf
			if (tag.isempty()) throw XMLReader::InvalidValue(tag);

			std::vector<string>& data = pp->value< std::vector<string> >();

			// Note that this parameter is read in item per item, not all at once!
			string s;
			tag.value(s);
			data.push_back(s);
		}
		break;
		case FE_PARAM_DOUBLE_MAPPED:
		{
			// get the model parameter
			FEParamDouble& p = pp->value<FEParamDouble>();

			// get the type
			const char* sztype = tag.AttributeValue("type", true);
			if (sztype == 0) sztype = "const";

			// allocate valuator
			FEScalarValuator* val = fecore_new<FEScalarValuator>(sztype, GetFEModel());
			if (val == nullptr) throw XMLReader::InvalidAttributeValue(tag, "type", sztype);

			// read the parameter list
			ReadParameterList(tag, val);

			// assign the valuator to the parameter
			p.setValuator(val);

			// do the initialization.
			// TODO: Is this a good place to do this?
			if (val->Init() == false) throw XMLReader::InvalidTag(tag);
		}
		break;
		case FE_PARAM_VEC3D_MAPPED:
		{
			// get the model parameter
			FEParamVec3& p = pp->value<FEParamVec3>();

			// get the type
			const char* sztype = tag.AttributeValue("type", true);
			if (sztype == 0) sztype = "vector";

			// allocate valuator
			FEVec3dValuator* val = fecore_new<FEVec3dValuator>(sztype, GetFEModel());
			if (val == nullptr) throw XMLReader::InvalidAttributeValue(tag, "type", sztype);

			// read the parameter list
			ReadParameterList(tag, val);

			// assign the valuator to the parameter
			p.setValuator(val);

			// do the initialization.
			// TODO: Is this a good place to do this?
			if (val->Init() == false) throw XMLReader::InvalidTag(tag);
		}
		break;
		case FE_PARAM_MAT3D_MAPPED:
		{
			// get the model parameter
			FEParamMat3d& p = pp->value<FEParamMat3d>();

			// get the type
			const char* sztype = tag.AttributeValue("type", true);
			if (sztype == 0) sztype = "const";

			// allocate valuator
			FEMat3dValuator* val = fecore_new<FEMat3dValuator>(sztype, GetFEModel());
			if (val == nullptr) throw XMLReader::InvalidAttributeValue(tag, "type", sztype);

			// read the parameter list
			ReadParameterList(tag, val);

			// assign the valuator to the parameter
			p.setValuator(val);

			// do the initialization.
			// TODO: Is this a good place to do this?
			if (val->Init() == false) throw XMLReader::InvalidTag(tag);
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
		case FE_PARAM_INT: value(tag, pp->pvalue<int   >(), pp->dim()); break;
		case FE_PARAM_DOUBLE: value(tag, pp->pvalue<double>(), pp->dim()); break;
		case FE_PARAM_DOUBLE_MAPPED:
		{
			std::vector<double> v(pp->dim());
			int m = tag.value(&v[0], pp->dim());
			if (m != pp->dim()) throw XMLReader::InvalidValue(tag);

			for (int i = 0; i < pp->dim(); ++i)
			{
				FEParamDouble& pi = pp->value<FEParamDouble>(i);
				pi = v[i];
			}
		}
		break;
        default: 
			assert(false);
			throw XMLReader::InvalidValue(tag);
			break;
		}
	}

	if (parseAttributes)
	{
		int nattr = tag.m_natt;
		for (int i = 0; i < nattr; ++i)
		{
			const char* szat = tag.m_att[i].m_szatt;
			if (pl.GetContainer()->SetParameterAttribute(*pp, szat, tag.m_att[i].m_szatv) == false)
			{
				// If we get here, the container did not understand the attribute.
				// If the attribute is a "lc", we interpret it as a load curve
				if (strcmp(szat, "lc") == 0)
				{
					int lc = atoi(tag.m_att[i].m_szatv) - 1;
					if (lc < 0) throw XMLReader::InvalidAttributeValue(tag, szat, tag.m_att[i].m_szatv);
					GetFEModel()->AttachLoadController(pp, lc);
				}
				/*			else
				{
				throw XMLReader::InvalidAttributeValue(tag, szat, tag.m_att[i].m_szatv);
				}
				*/
			}
			// This is not true. Parameters can have attributes that are used for other purposes. E.g. The local fiber option.
			//		else felog.printf("WARNING: attribute \"%s\" of parameter \"%s\" ignored (line %d)\n", szat, tag.Name(), tag.m_ncurrent_line-1);
		}
	}

	// give the parameter container a chance to do additional processing
	pl.GetContainer()->SetParameter(*pp);

	return true;
}

//-----------------------------------------------------------------------------
void FEFileSection::ReadAttributes(XMLTag& tag, FECoreBase* pc)
{
	FEParameterList& pl = pc->GetParameterList();

	// process all the other attributes
	for (int i = 0; i<tag.m_natt; ++i)
	{
		XMLAtt& att = tag.m_att[i];
		const char* szatt = att.name();
		const char* szval = att.cvalue();
		if ((att.m_bvisited == false) && szatt && szval)
		{
			if      (strcmp("name", szatt) == 0) pc->SetName(szval);
			else if (strcmp("id"  , szatt) == 0) pc->SetID(atoi(szval));
			else if (strcmp("lc"  , szatt) == 0) { /* don't do anything. Loadcurves are processed elsewhere. */ }
			else
			{
				FEParam* param = pl.FindFromName(szatt);
				if (param && (param->GetFlags() & FE_PARAM_ATTRIBUTE))
				{
					switch (param->type())
					{
					case FE_PARAM_INT:
					{
						if (param->enums() == nullptr)
							param->value<int>() = atoi(szval);
						else
						{
							if (parseEnumParam(param, szval) == false) throw XMLReader::InvalidAttributeValue(tag, szatt, szval);
						}
						break;
					}
					case FE_PARAM_DOUBLE: param->value<double>() = atof(szval); break;
					default:
						throw XMLReader::InvalidAttributeValue(tag, szatt, szval);
					}
				}
				else
				{
					throw XMLReader::InvalidAttribute(tag, szatt);
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! This function parses a parameter list
bool FEFileSection::ReadParameter(XMLTag& tag, FECoreBase* pc, const char* szparam, bool parseAttributes)
{
	// get the parameter list
	FEParameterList& pl = pc->GetParameterList();

	// see if we can find this parameter
	if (ReadParameter(tag, pl, szparam, pc, parseAttributes) == false)
	{
		// if we get here, the parameter is not found.
		// See if the parameter container has defined a property of this name
		int n = pc->FindPropertyIndex(tag.Name());
		if (n >= 0)
		{
			FEProperty* prop = pc->PropertyClass(n);
			const char* sztype = tag.AttributeValue("type", true);

			// If the type attribute is omitted we assume the tag's name is the type
			if (sztype == 0) sztype = tag.Name();

			// try to allocate the class
			FECoreBase* pp = fecore_new<FECoreBase>(prop->GetClassID(), sztype, GetFEModel());
			if (pp == nullptr) throw XMLReader::InvalidAttributeValue(tag, "type", sztype);

			prop->SetProperty(pp);

			// read the property data
			if (tag.isleaf() == false)
			{
				ReadParameterList(tag, pp);
			}
			else if (tag.isempty() == false)
			{
				// There should be a parameter with the same name as the type
				if (ReadParameter(tag, pp->GetParameterList(), sztype, pp) == false)
					throw XMLReader::InvalidValue(tag);
			}
			return true;
		}
		else return false;
	}
	return true;
}

//-----------------------------------------------------------------------------
void FEFileSection::ReadParameterList(XMLTag& tag, FEParameterList& pl)
{
	// Make sure there is something to read
	if (tag.isleaf() || tag.isempty()) return;

	// parse the child tags
	++tag;
	do
	{
		if (ReadParameter(tag, pl, 0, 0) == false) throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());
}

//-----------------------------------------------------------------------------
void FEFileSection::ReadParameterList(XMLTag& tag, FECoreBase* pc)
{
	// parse attributes first
	ReadAttributes(tag, pc);

	// process the parameter lists
	if (!tag.isleaf())
	{
		++tag;
		do
		{
			if (ReadParameter(tag, pc) == false) throw XMLReader::InvalidTag(tag);
			++tag;
		}
		while (!tag.isend());
	}
	else if (tag.isempty() == false)
	{
		// there should be one parameter with the same name as the tag
		// Notice that we don't process attributes, since this situation should be synonymous 
		// to an additional tag that defines the parameter.
		if (ReadParameter(tag, pc, 0, false) == false)
		{
			// try a parameter with the type string as name
			if (ReadParameter(tag, pc, pc->GetTypeStr(), false) == false)
				throw XMLReader::InvalidTag(tag);
		}
	}

	// validate the class
	int NP = pc->PropertyClasses();
	for (int i = 0; i<NP; ++i)
	{
		FEProperty* pi = pc->PropertyClass(i);
		bool a = pi->IsRequired();
		bool b = (pi->size() == 0);
		if (a && b)
		{
			throw FEBioImport::MissingProperty(pc->GetName(), pi->GetName());
		}
	}

}

//-----------------------------------------------------------------------------
FEFileSectionMap::~FEFileSectionMap()
{
	Clear();
}

//-----------------------------------------------------------------------------
void FEFileSectionMap::Clear()
{
	// clear the map
	FEFileSectionMap::iterator is;
	for (is = begin(); is != end(); ++is)
	{
		FEFileSection* ps = is->second; delete ps;
	}
	clear();
}

//-----------------------------------------------------------------------------
void FEFileSectionMap::Parse(XMLTag& tag)
{
	if (tag.isleaf()) return;

	++tag;
	do
	{
		std::map<string, FEFileSection*>::iterator is = find(tag.Name());
		if (is != end()) is->second->Parse(tag);
		else throw XMLReader::InvalidTag(tag);

		++tag;
	}
	while (!tag.isend());
}

//-----------------------------------------------------------------------------
//! class constructor
FEFileImport::FEFileImport()
{
	m_fp = 0;
	m_szerr[0] = 0;
	m_fem = 0;
	m_builder = 0;
	m_nversion = 0;
}

//-----------------------------------------------------------------------------
// class destructor. Closes file on call.
FEFileImport::~FEFileImport()
{
	// make sure to close the file
	Close();
}

//-----------------------------------------------------------------------------
//! Open a file and store the file name and file pointer.
bool FEFileImport::Open(const char* szfile, const char* szmode)
{
	m_fp = fopen(szfile, szmode);
	if (m_fp == 0) return false;

	strcpy(m_szfile, szfile);

	return true;
}

//-----------------------------------------------------------------------------
//! Close the file
void FEFileImport::Close()
{
	if (m_fp)
	{
		fclose(m_fp);
		m_fp = 0;
	}

	m_fem = 0;
	delete m_builder;
	m_builder = 0;
}

//-----------------------------------------------------------------------------
//! parse the file
bool FEFileImport::ParseFile(XMLTag& tag)
{
	// parse the file
	try
	{
		m_map.Parse(tag);
	}
	catch (XMLReader::Error& e)
	{
		fprintf(stderr, "FATAL ERROR: %s (line %d)\n", e.GetErrorString(), tag.m_ncurrent_line);
		return false;
	}
	catch (...)
	{
		fprintf(stderr, "FATAL ERROR: unknown exception occured while parsing file.\n\n");
		return false;
	}

	return true;
}

//-----------------------------------------------------------------------------
//! Get the FEModel that is being processed. This is the model that was passed in the Load function
FEModel* FEFileImport::GetFEModel()
{
	return m_fem;
}

//-----------------------------------------------------------------------------
//! Get the model builder
FEModelBuilder* FEFileImport::GetBuilder()
{
	return m_builder;
}

//-----------------------------------------------------------------------------
// return the file path
const char* FEFileImport::GetFilePath() { return m_szpath; }

//-----------------------------------------------------------------------------
// set file version
void FEFileImport::SetFileVerion(int nversion)
{
	FECoreKernel& fecore = FECoreKernel::GetInstance();
	fecore.SetSpecID(nversion);
	m_nversion = nversion;
}

//-----------------------------------------------------------------------------
// get file version
int FEFileImport::GetFileVersion() const
{
	return m_nversion;
}

//-----------------------------------------------------------------------------
//! Get the error message. Errors message are stored when calling the errf function.
void FEFileImport::GetErrorMessage(char* szerr)
{
	strcpy(szerr, m_szerr);
}

//-----------------------------------------------------------------------------
//! Call this function to report an error message. The user can retrieve the 
//! error message with the GetErrorMessage member function.
bool FEFileImport::errf(const char* szerr, ...)
{
	// get a pointer to the argument list
	va_list	args;

	// copy to string
	va_start(args, szerr);
	vsprintf(m_szerr, szerr, args);
	va_end(args);

	// close the file
	Close();

	return false;
}

//-----------------------------------------------------------------------------
void FEFileImport::ClearFileParams()
{
	m_Param.clear();
}

//-----------------------------------------------------------------------------
FEFileParam* FEFileImport::FindFileParameter(const char* sz)
{
	for (size_t i = 0; i<m_Param.size(); ++i)
	{
		FEFileParam& p = m_Param[i];
		if (strcmp(p.m_szname, sz) == 0) return &p;
	}
	return 0;
}

//-----------------------------------------------------------------------------
void FEFileImport::AddFileParameter(const char* szname, const char* szval)
{
	FEFileParam p;
	strcpy(p.m_szname, szname);
	strcpy(p.m_szval, szval);
	m_Param.push_back(p);
}
