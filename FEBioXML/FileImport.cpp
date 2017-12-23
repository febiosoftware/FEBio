// FileImport.cpp: implementation of the FEFileImport class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FileImport.h"
#include <FECore/Image.h>
#include <FECore/FEDataArray.h>
#include <FECore/FEFunction1D.h>
#include <FECore/FEModel.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>

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
//! This function parese a parameter list
bool FEFileSection::ReadParameter(XMLTag& tag, FEParameterList& pl, const char* szparam)
{
	// see if we can find this parameter
	FEParam* pp = pl.Find((szparam == 0 ? tag.Name() : szparam));
	if (pp == 0) return false;

	if (pp->dim() == 1)
	{
		switch (pp->type())
		{
		case FE_PARAM_DOUBLE: value(tag, pp->value<double  >()); break;
		case FE_PARAM_INT: value(tag, pp->value<int     >()); break;
		case FE_PARAM_BOOL: value(tag, pp->value<bool    >()); break;
		case FE_PARAM_VEC3D: value(tag, pp->value<vec3d   >()); break;
		case FE_PARAM_MAT3D: value(tag, pp->value<mat3d   >()); break;
		case FE_PARAM_MAT3DS: value(tag, pp->value<mat3ds  >()); break;
		case FE_PARAM_TENS3DRS: value(tag, pp->value<tens3drs>()); break;
		case FE_PARAM_STRING: value(tag, pp->cvalue()); break;
		case FE_PARAM_IMAGE_3D:
		{
			const char* szfile = tag.AttributeValue("file");
			++tag;
			int n[3] = { 0 };
			do
			{
				if (tag == "size") tag.value(n, 3);
				else throw XMLReader::InvalidTag(tag);
				++tag;
			} while (!tag.isend());
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
				FEDataArray* pdata = GetFEModel()->FindDataArray(szmap);
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
					FEDataArray* pdata = GetFEModel()->FindDataArray(szmap);
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
						double v[2] = { 0 };
						tag.value(v, 2);
						map.set(vec2d(v[0], v[1]));
					}
					else if (map.DataType() == FE_VEC3D)
					{
						double v[3] = { 0 };
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
		case FE_PARAM_INT: value(tag, pp->pvalue<int   >(), pp->dim()); break;
		case FE_PARAM_DOUBLE: value(tag, pp->pvalue<double>(), pp->dim()); break;
        default: break;
		}
	}

	int nattr = tag.m_natt;
	for (int i = 0; i<nattr; ++i)
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
				switch (pp->type())
				{
				case FE_PARAM_INT: pp->SetLoadCurve(lc); break;
				case FE_PARAM_BOOL: pp->SetLoadCurve(lc); break;
				case FE_PARAM_DOUBLE: pp->SetLoadCurve(lc, pp->value<double>()); break;
				case FE_PARAM_VEC3D: pp->SetLoadCurve(lc, pp->value<vec3d >()); break;
				case FE_PARAM_FUNC1D: break; // don't do anything for 1D functions since the lc attribute is already processed.
				default:
					assert(false);
				}
			}
			/*			else
			{
			throw XMLReader::InvalidAttributeValue(tag, szat, tag.m_att[i].m_szatv);
			}
			*/
		}
		// This is not true. Parameters can have attributes that are used for other purposed. E.g. The local fiber option.
		//		else felog.printf("WARNING: attribute \"%s\" of parameter \"%s\" ignored (line %d)\n", szat, tag.Name(), tag.m_ncurrent_line-1);
	}

	// give the parameter container a chance to do additional processing
	pl.GetContainer()->SetParameter(*pp);

	return true;
}

//-----------------------------------------------------------------------------
//! This function parses a parameter list
bool FEFileSection::ReadParameter(XMLTag& tag, FECoreBase* pc, const char* szparam)
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
					} while (!tag.isend());
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
void FEFileSection::ReadParameterList(XMLTag& tag, FEParameterList& pl)
{
	// Make sure there is something to read
	if (tag.isleaf() || tag.isempty()) return;

	// parse the child tags
	++tag;
	do
	{
		if (ReadParameter(tag, pl) == false) throw XMLReader::InvalidTag(tag);
		++tag;
	} while (!tag.isend());
}

//-----------------------------------------------------------------------------
void FEFileSection::ReadParameterList(XMLTag& tag, FECoreBase* pc)
{
	FEParameterList& pl = pc->GetParameterList();
	ReadParameterList(tag, pl);
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
bool FEFileImport::Load(FEModel& fem, const char* szfile)
{
	m_fem = &fem;
	m_builder = new FEModelBuilder(fem);

	return Parse(szfile);
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
