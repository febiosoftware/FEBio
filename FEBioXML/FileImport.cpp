/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/

#include "stdafx.h"
#include "FileImport.h"
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
#include <FECore/FEPointFunction.h>
#include <FECore/FEGlobalData.h>
#include <FECore/log.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <sstream>
#include "FEBioImport.h"

#ifndef WIN32
#define strnicmp strncasecmp
#endif

//-----------------------------------------------------------------------------
FEObsoleteParamHandler::FEObsoleteParamHandler(XMLTag& tag, FECoreBase* pc) : m_pc(pc) 
{
	m_root = tag.Name();
}

void FEObsoleteParamHandler::AddParam(const char* oldName, const char* newName, FEParamType paramType)
{
	FEObsoleteParam p = { oldName, newName, paramType };
	m_param.push_back(p);
}

bool FEObsoleteParamHandler::ProcessTag(XMLTag& tag)
{
	std::string path = tag.relpath(m_root.c_str());
	for (int i = 0; i < m_param.size(); ++i)
	{
		if (path == m_param[i].oldName)
		{
			m_param[i].readIn = true;

			switch (m_param[i].paramType)
			{
			case FE_PARAM_INVALID: break; // parameter will be ignored
			case FE_PARAM_BOOL   : tag.value(m_param[i].bVal); break;
			case FE_PARAM_INT    : tag.value(m_param[i].iVal); break;
			case FE_PARAM_DOUBLE : tag.value(m_param[i].gVal); break;
			default:
				assert(false);
			}
			return true;
		}
	}

	return false;
}

void FEObsoleteParamHandler::MapParameters()
{
	FEModel* fem = m_pc->GetFEModel();
	feLogEx(fem, "\n");
	for (FEObsoleteParam& p : m_param)
	{
		if (p.readIn)
		{
			if (p.newName == nullptr)
			{
				feLogWarningEx(fem, "Obsolete parameter %s is ignored!", p.oldName);
			}
			else
			{
				ParamString ps(p.newName);
				FEParam* pp = m_pc->FindParameter(ps);
				if (pp == nullptr)
				{
					feLogErrorEx(fem, "Failed to map obsolete parameter %s. Could not find new parameter %s", p.oldName, p.newName);
				}
				else
				{
					if (pp->type() == p.paramType)
					{
						switch (pp->type())
						{
						case FE_PARAM_BOOL  : pp->value<bool>  () = p.bVal; break;
						case FE_PARAM_INT   : pp->value<int>   () = p.iVal; break;
						case FE_PARAM_DOUBLE: pp->value<double>() = p.gVal; break;
						}
						feLogEx(fem, "Successfully mapped obsolete parameter %s to %s\n", p.oldName, p.newName);
					}
					else if ((pp->type() == FE_PARAM_DOUBLE_MAPPED) && (p.paramType == FE_PARAM_DOUBLE))
					{
						FEParamDouble& v = pp->value<FEParamDouble>();
						v = p.gVal;
						feLogEx(fem, "Successfully mapped obsolete parameter %s to %s\n", p.oldName, p.newName);
					}
					else
					{
						feLogErrorEx(fem, "Failed to map obsolete parameter %s. New parameter %s has different type.", p.oldName, p.newName);
					}
				}
			}
		}
	}
}

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
FEModel* FEFileSection::GetFEModel() { return &GetBuilder()->GetFEModel(); }

//-----------------------------------------------------------------------------
FEModelBuilder* FEFileSection::GetBuilder() { return m_pim->GetBuilder(); }

//-----------------------------------------------------------------------------
void FEFileSection::SetInvalidTagHandler(FEInvalidTagHandler* ith)
{
	m_ith = ith;
}

//-----------------------------------------------------------------------------
void FEFileSection::value(XMLTag& tag, int& n)
{
	n = atoi(tag.szvalue());
}

//-----------------------------------------------------------------------------
void FEFileSection::value(XMLTag& tag, double& g)
{
	g = atof(tag.szvalue());
}

//-----------------------------------------------------------------------------
void FEFileSection::value(XMLTag& tag, bool& b)
{
	const char* sz = tag.szvalue();
	int n = 0;
	sscanf(sz, "%d", &n);
	b = (n != 0);
}

//-----------------------------------------------------------------------------
void FEFileSection::value(XMLTag& tag, vec3d& v)
{
	const char* sz = tag.szvalue();
	int n = sscanf(sz, "%lg,%lg,%lg", &v.x, &v.y, &v.z);
	if (n != 3) throw XMLReader::XMLSyntaxError(tag.m_nstart_line);
}

//-----------------------------------------------------------------------------
void FEFileSection::value(XMLTag& tag, mat3d& m)
{
	const char* sz = tag.szvalue();
	double xx, xy, xz, yx, yy, yz, zx, zy, zz;
	int n = sscanf(sz, "%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg", &xx, &xy, &xz, &yx, &yy, &yz, &zx, &zy, &zz);
	if (n != 9) throw XMLReader::XMLSyntaxError(tag.m_nstart_line);
	m = mat3d(xx, xy, xz, yx, yy, yz, zx, zy, zz);
}

//-----------------------------------------------------------------------------
void FEFileSection::value(XMLTag& tag, mat3ds& m)
{
	const char* sz = tag.szvalue();
	double x, y, z, xy, yz, xz;
	int n = sscanf(sz, "%lg,%lg,%lg,%lg,%lg,%lg", &x, &y, &z, &xy, &yz, &xz);
	if (n != 6) throw XMLReader::XMLSyntaxError(tag.m_nstart_line);
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
	const char* sz = tag.szvalue();
	strcpy(szstr, sz);
}

//-----------------------------------------------------------------------------
void FEFileSection::value(XMLTag& tag, std::string& v)
{
	v = tag.szvalue();
}

//-----------------------------------------------------------------------------
void FEFileSection::value(XMLTag& tag, std::vector<int>& v)
{
	v.clear();
	const char* sz = tag.szvalue();
	do
	{
		const char* sze = strchr(sz, ',');

		int d = atoi(sz);
		v.push_back(d);

		if (sze) sz = sze + 1;
		else sz = nullptr;

	}
	while (sz);
}

//-----------------------------------------------------------------------------
void FEFileSection::value(XMLTag& tag, std::vector<double>& v)
{
	v.clear();
	const char* sz = tag.szvalue();
	do
	{
		const char* sze = strchr(sz, ',');

		const char* szc = strchr(sz, ':');
		if ((szc == nullptr) || (sze && (szc > sze)))
		{
			double d = atof(sz);
			v.push_back(d);
		}
		else
		{
			double d0 = atof(sz);
			double d1 = atof(szc + 1);
			double di = 1.0;
			// see if there is a second colon
			const char* szc2 = strchr(szc + 1, ':');
			if (szc2 && ((sze == nullptr) || (szc2 < sze))) di = atof(szc2 + 1);
			if (di > 0)
			{
				const double eps = 1e-12;
				for (double d = d0; d <= d1 + eps; d += di)
				{
					v.push_back(d);
				}
			}
			else throw XMLReader::InvalidValue(tag);
		}

		if (sze) sz = sze + 1;
		else sz = nullptr;
	} 
	while (sz);
}

//-----------------------------------------------------------------------------
int FEFileSection::value(XMLTag& tag, int* pi, int n)
{
	const char* sz = tag.szvalue();
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
	const char* sz = tag.szvalue();
	int nr = 0;
	for (int i = 0; i < n; ++i)
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
int enumValue(const char* val, const char* szenum)
{
	if ((val == nullptr) || (szenum == nullptr)) return -1;

	// get the string's length. 
	// there could be a comma, so correct for that.
	size_t L = strlen(val);
	const char* c = strchr(val, ',');
	if (c) L = c - val;

	const char* ch = szenum;

	int n = 0;
	while (ch && *ch)
	{
		size_t l = strlen(ch);
		int nval = n;
		// see if the value of the enum is overridden
		const char* ce = strrchr(ch, '=');
		if (ce)
		{
			l = ce - ch;
			nval = atoi(ce + 1);
		}

		if ((L==l) && (strnicmp(ch, val, l) == 0))
		{
			return nval;
		}
		ch = strchr(ch, '\0');
		if (ch) ch++;
		n++;
	}
	return -1;
}

//-----------------------------------------------------------------------------
// helper function to see if a string is a number
bool is_number(const char* sz)
{
	char* cend;
	double tmp = strtod(sz, &cend);
	return ((cend == nullptr) || (cend[0] == 0));
}

//-----------------------------------------------------------------------------
bool FEFileSection::parseEnumParam(FEParam* pp, const char* val)
{
	// get the enums
	const char* szenums = pp->enums();
	if (szenums == nullptr) return false;

	// special cases
	if (szenums[0] == '$')
	{
		FEModel* fem = GetFEModel();

		char var[256] = { 0 };
		const char* chl = strchr(szenums, '('); assert(chl);
		const char* chr = strchr(szenums, ')'); assert(chr);
		strncpy(var, chl + 1, chr - chl - 1);

		if (strncmp(var, "dof_list", 8) == 0)
		{
			DOFS& dofs = GetFEModel()->GetDOFS();
			const char* szvar = nullptr;
			if (var[8] == ':') szvar = var + 9;

			if (pp->type() == FE_PARAM_INT)
			{
				int ndof = dofs.GetDOF(val, szvar);
				if (ndof < 0) return false;
				pp->value<int>() = ndof;
				return true;
			}
			else if (pp->type() == FE_PARAM_STD_VECTOR_INT)
			{
				std::vector<int>& v = pp->value<std::vector<int> >();
				return dofs.ParseDOFString(val, v, szvar);
			}
			else return false;
		}
		else if (strcmp(var, "solutes") == 0)
		{
			int n = -1;
			if (is_number(val)) n = atoi(val);
			else
			{
				FEGlobalData* pd = fem->FindGlobalData(val);
				if (pd == nullptr) return false;
				n = pd->GetID(); assert(n > 0);
			}
				
			pp->value<int>() = n;
			return true;
		}
		else if (strcmp(var, "sbms") == 0)
		{
			int n = -1;
			if (is_number(val)) n = atoi(val);
			else
			{
				FEGlobalData* pd = fem->FindGlobalData(val);
				if (pd == nullptr) return false;
				n = pd->GetID(); assert(n > 0);
			}

			pp->value<int>() = n;
			return true;
		}
		else if (strcmp(var, "species") == 0)
		{
			int n = -1;
			if (is_number(val)) n = atoi(val);
			else
			{
				// NOTE: This assumes that the solutes are defined before the SBMS!
				int m = fem->FindGlobalDataIndex(val);
				if (m == -1) return false;
				n = m + 1;
			}
			pp->value<int>() = n;
			return true;
		}
		else if (strcmp(var, "rigid_materials") == 0)
		{
			if (is_number(val))
			{
				int n = atoi(val);
				pp->value<int>() = n;
			}
			else
			{
				FEMaterial* mat = fem->FindMaterial(val);
				if (mat == nullptr) return false;
				int n = mat->GetID();
				pp->value<int>() = n;
			}
			return true;
		}
		else return false;
	}
	else if (strncmp(szenums, "@factory_list", 13) == 0)
	{
		int classID = atoi(szenums + 14);

		FECoreKernel& fecore = FECoreKernel::GetInstance();
		for (int i = 0; i < fecore.FactoryClasses(); ++i)
		{
			const FECoreFactory* fac = fecore.GetFactoryClass(i);
			if (fac->GetSuperClassID() == classID)
			{
				if (strcmp(fac->GetTypeStr(), val) == 0)
				{
					pp->value<int>() = i;
					return true;
				}
			}
		}
		return false;
	}

	switch (pp->type())
	{
	case FE_PARAM_INT:
	{
		int n = enumValue(val, szenums);
		if (n != -1) pp->value<int>() = n;
		else
		{
			// see if the value is an actual number
			if (is_number(val))
			{
				n = atoi(val);
				pp->value<int>() = n;
			}
		}
		return (n != -1);
	}
	break;
	case FE_PARAM_STD_VECTOR_INT:
	{
		std::vector<int>& v = pp->value<std::vector<int> >();
		v.clear();
		const char* tmp = val;
		while (tmp)
		{
			int n = enumValue(tmp, szenums);
			v.push_back(n);
			tmp = strchr(tmp, ',');
			if (tmp) tmp++;
		}

		return true;
	}
	break;
	};

	return false;
}

//-----------------------------------------------------------------------------
std::vector<std::string> split_string(const std::string& s, char delim)
{
	std::stringstream ss(s);
	std::string tmp;
	std::vector<std::string> vs;
	while (std::getline(ss, tmp, delim))
	{
		vs.push_back(tmp);
	}
	return vs;
}

//-----------------------------------------------------------------------------
//! This function parses a parameter list
bool FEFileSection::ReadParameter(XMLTag& tag, FEParameterList& pl, const char* szparam, FECoreBase* pc, bool parseAttributes)
{
	FEParam* pp = nullptr;
	if (tag == "add_param")
	{
		// get the name and value
		double v = 0.0; tag.value(v);
		const char* szname = tag.AttributeValue("name");

		// make sure this parameter does not exist yet
		pp = pl.FindFromName(szname);
		if (pp) throw XMLReader::InvalidTag(tag);

		// add a new user parameter
		pp = pl.AddParameter(new double(v), FE_PARAM_DOUBLE, 1, strdup(szname));
		pp->SetFlags(FEParamFlag::FE_PARAM_USER);
	}
	else
	{
		// see if we can find this parameter
		pp = pl.FindFromName((szparam == 0 ? tag.Name() : szparam));
	}
	if (pp == 0) return false;

	if (pp->dim() == 1)
	{
		switch (pp->type())
		{
			case FE_PARAM_DOUBLE: 
			{
				// make sure the type attribute is not defined 
				// This most likely means a user thinks this parameter can be mapped
				// but the corresponding parameter is not a FEModelParam
				const char* sztype = tag.AttributeValue("type", true);
				if (sztype) throw XMLReader::InvalidAttribute(tag, "type");

				value(tag, pp->value<double  >());
			}
			break;
		case FE_PARAM_INT: 
			{
				const char* szenum = pp->enums();
				if (szenum == 0)
				{
					value(tag, pp->value<int>());
				}
				else
				{
					bool bfound = parseEnumParam(pp, tag.szvalue());
					if (bfound == false)
					{
						if ((m_ith == nullptr) || (m_ith->ProcessTag(tag) == false)) {
							throw XMLReader::InvalidValue(tag);
						}
					}
				}
			}
			break;
		case FE_PARAM_STD_VECTOR_INT:
			{
				const char* szenum = pp->enums();
				if (szenum == 0)
				{
					std::vector<int>& v = pp->value< std::vector<int> >();
					value(tag, v);
				}
				else
				{
					bool bfound = parseEnumParam(pp, tag.szvalue());
					if (bfound == false) throw XMLReader::InvalidValue(tag);
				}
			}
			break;
		case FE_PARAM_STD_VECTOR_DOUBLE:
			{
				std::vector<double>& v = pp->value< std::vector<double> >();
				value(tag, v);
			}
			break;
		case FE_PARAM_BOOL: value(tag, pp->value<bool    >()); break;
		case FE_PARAM_VEC3D: 
			{
				// make sure the type attribute is not defined 
				// This most likely means a user thinks this parameter can be mapped
				// but the corresponding parameter is not a FEModelParam
//				const char* sztype = tag.AttributeValue("type", true);
//				if (sztype) throw XMLReader::InvalidAttribute(tag, "type");

				value(tag, pp->value<vec3d   >());
			}
			break;
		case FE_PARAM_MAT3D: value(tag, pp->value<mat3d   >()); break;
		case FE_PARAM_MAT3DS: value(tag, pp->value<mat3ds  >()); break;
		case FE_PARAM_TENS3DRS: value(tag, pp->value<tens3drs>()); break;
		case FE_PARAM_STRING: value(tag, pp->cvalue()); break;
		case FE_PARAM_STD_STRING: value(tag, pp->value<string>()); break;
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

				FESurfaceMap* pdata = dynamic_cast<FESurfaceMap*>(GetFEModel()->GetMesh().FindDataMap(szmap));
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

					FENodeDataMap* pdata = dynamic_cast<FENodeDataMap*>(GetFEModel()->GetMesh().FindDataMap(szmap));
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
			std::vector<vec2d>& data = pp->value< std::vector<vec2d> >();
			data.clear();

			double d[2];
			++tag;
			do
			{
				int nread = tag.value(d, 2);
				if (nread != 2) throw XMLReader::InvalidValue(tag);
				data.push_back(vec2d(d[0], d[1]));
				++tag;
			}
			while (!tag.isend());
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

			// if the type is not specified, we'll try to determine if 
			// it's a math expression or a const
			if (sztype == 0)
			{
				const char* szval = tag.szvalue();
				sztype = (is_number(szval) ? "const" : "math");
			}

			// allocate valuator
			FEScalarValuator* val = fecore_new<FEScalarValuator>(sztype, GetFEModel());
			if (val == nullptr) throw XMLReader::InvalidAttributeValue(tag, "type", sztype);

			// Figure out the item list
			FEItemList* itemList = nullptr;
			if (dynamic_cast<FESurfaceLoad*>(pc))
			{
				FESurfaceLoad* psl = dynamic_cast<FESurfaceLoad*>(pc);
				itemList = psl->GetSurface().GetFacetSet();
			}

			p.SetItemList(itemList);

			// mapped values require special treatment
			// The value is just the name of the map, but the problem is that 
			// these maps may not be defined yet.
			// So, we add them to the FEBioModel, which will process mapped 
			// parameters after the rest of the file is processed
			if (strcmp(sztype, "map") == 0)
			{
				GetBuilder()->AddMappedParameter(pp, pc, tag.szvalue());
			}
			else {
				// read the parameter list
				ReadParameterList(tag, val);
			}

			// assign the valuator to the parameter
			p.setValuator(val);
		}
		break;
		case FE_PARAM_VEC3D_MAPPED:
		{
			// get the model parameter
			FEParamVec3& p = pp->value<FEParamVec3>();

			// get the type
			const char* sztype = tag.AttributeValue("type", true);

			// ignore "user" types
			if (sztype && strcmp(sztype, "user") == 0)
			{
				return true;
			}

			if (sztype == 0) sztype = "vector";

			// allocate valuator
			FEVec3dValuator* val = fecore_new<FEVec3dValuator>(sztype, GetFEModel());
			if (val == nullptr) throw XMLReader::InvalidAttributeValue(tag, "type", sztype);

			// mapped values require special treatment
			// The value is just the name of the map, but the problem is that 
			// these maps may not be defined yet.
			// So, we add them to the FEBioModel, which will process mapped 
			// parameters after the rest of the file is processed
			if (strcmp(sztype, "map") == 0)
			{
				GetBuilder()->AddMappedParameter(pp, pc, tag.szvalue());
			}
			else {
				// read the parameter list
				ReadParameterList(tag, val);
			}

			// assign the valuator to the parameter
			p.setValuator(val);
		}
		break;
		case FE_PARAM_MAT3D_MAPPED:
		{
			// get the model parameter
			FEParamMat3d& p = pp->value<FEParamMat3d>();

			// get the type
			const char* sztype = tag.AttributeValue("type", true);
			if (sztype == 0) sztype = "const";

			// ignore user type
			if (sztype && (strcmp(sztype, "user") == 0)) return true;

			// allocate valuator
			FEMat3dValuator* val = fecore_new<FEMat3dValuator>(sztype, GetFEModel());
			if (val == nullptr) throw XMLReader::InvalidAttributeValue(tag, "type", sztype);

			// mapped values require special treatment
			// The value is just the name of the map, but the problem is that 
			// these maps may not be defined yet.
			// So, we add them to the FEBioModel, which will process mapped 
			// parameters after the rest of the file is processed
			if (strcmp(sztype, "map") == 0)
			{
				GetBuilder()->AddMappedParameter(pp, pc, tag.szvalue());
			}
			else {
				// read the parameter list
				ReadParameterList(tag, val);
			}

			// assign the valuator to the parameter
			p.setValuator(val);
		}
		break;
		case FE_PARAM_MAT3DS_MAPPED:
		{
			// get the model parameter
			FEParamMat3ds& p = pp->value<FEParamMat3ds>();

			// get the type
			const char* sztype = tag.AttributeValue("type", true);
			if (sztype == 0) sztype = "const";

			// allocate valuator
			FEMat3dsValuator* val = fecore_new<FEMat3dsValuator>(sztype, GetFEModel());
			if (val == nullptr) throw XMLReader::InvalidAttributeValue(tag, "type", sztype);

			// mapped values require special treatment
			// The value is just the name of the map, but the problem is that 
			// these maps may not be defined yet.
			// So, we add them to the FEBioModel, which will process mapped 
			// parameters after the rest of the file is processed
			if (strcmp(sztype, "map") == 0)
			{
				GetBuilder()->AddMappedParameter(pp, pc, tag.szvalue());
			}
			else {
				// read the parameter list
				ReadParameterList(tag, val);
			}

			// assign the valuator to the parameter
			p.setValuator(val);

			// do the initialization.
			// TODO: Is this a good place to do this?
			if (val->Init() == false) throw XMLReader::InvalidTag(tag);
		}
		break;		default:
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
			if (tag.isleaf())
			{
				// find the type attribute
				const char* sztype = tag.AttributeValue("type", true);
				if (sztype == nullptr) sztype = "const";

				if (strcmp(sztype, "const") == 0)
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
				else if (strcmp(sztype, "math") == 0)
				{
					string sval = tag.szvalue();

					vector<string> s = split_string(sval, ',');
					if (s.size() != pp->dim()) throw XMLReader::InvalidValue(tag);

					for (int i = 0; i < pp->dim(); ++i)
					{
						FEParamDouble& pi = pp->value<FEParamDouble>(i);
						FEMathValue* v = fecore_alloc(FEMathValue, GetFEModel());
						v->setMathString(s[i]);
						pi.setValuator(v);
					}
				}
				else if (strcmp(sztype, "map") == 0)
				{ 
					string sval = tag.szvalue();

					vector<string> s = split_string(sval, ',');
					if (s.size() != pp->dim()) throw XMLReader::InvalidValue(tag);

					for (int i = 0; i < pp->dim(); ++i)
					{
						FEParamDouble& pi = pp->value<FEParamDouble>(i);

						// allocate valuator
						FEScalarValuator* val = fecore_alloc(FEMappedValue, GetFEModel());

						// Figure out the item list
						FEItemList* itemList = nullptr;
						if (dynamic_cast<FESurfaceLoad*>(pc))
						{
							FESurfaceLoad* psl = dynamic_cast<FESurfaceLoad*>(pc);
							itemList = psl->GetSurface().GetFacetSet();
						}
						pi.SetItemList(itemList);

						// mapped values require special treatment
						// The value is just the name of the map, but the problem is that 
						// these maps may not be defined yet.
						// So, we add them to the FEBioModel, which will process mapped 
						// parameters after the rest of the file is processed
						GetBuilder()->AddMappedParameter(pp, pc, tag.szvalue(), i);

						// assign the valuator to the parameter
						pi.setValuator(val);
					}
				}
				else throw XMLReader::InvalidAttributeValue(tag, "type", sztype);
			}
			else
			{
				int n = 0;
				++tag;
				do {
					// find the type attribute
					const char* sztype = tag.AttributeValue("type", true);
					if (sztype == nullptr)
					{
						// if the type is not specified, we'll try to determine if 
						// it's a math expression or a const
						const char* szval = tag.szvalue();
						sztype = (is_number(szval) ? "const" : "math");
					}

					if (strcmp(sztype, "const") == 0)
					{
						double v = 0.0;
						tag.value(v);
						FEParamDouble& pi = pp->value<FEParamDouble>(n);
						pi = v;
					}
					else if (strcmp(sztype, "math") == 0)
					{
						string sval = tag.szvalue();

						FEParamDouble& pi = pp->value<FEParamDouble>(n);
						FEMathValue* v = fecore_alloc(FEMathValue, GetFEModel());
						v->setMathString(sval);
						pi.setValuator(v);

						// do the initialization.
						// TODO: Is this a good place to do this?
						if (v->Init() == false) throw XMLReader::InvalidTag(tag);
					}
					else throw XMLReader::InvalidAttributeValue(tag, "type", sztype);
					++tag;
					n++;
					if (n > pp->dim()) throw XMLReader::InvalidTag(tag);
				}
				while (!tag.isend());
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
		for (XMLAtt& att : tag.m_att)
		{
			const char* szat = att.m_name.c_str();
			// If we get here, the container did not understand the attribute.
			// If the attribute is a "lc", we interpret it as a load curve
			if (strcmp(szat, "lc") == 0)
			{
				int lc = atoi(att.m_val.c_str()) - 1;
				if (lc < 0) throw XMLReader::InvalidAttributeValue(tag, szat, att.m_val.c_str());

				// make sure the parameter is volatile
				if (pp->IsVolatile() == false)
				{
					throw XMLReader::InvalidAttribute(tag, szat);
				}

				GetFEModel()->AttachLoadController(pp, lc);
			}
/*			else
			{
				throw XMLReader::InvalidAttributeValue(tag, szat, tag.m_att[i].m_szatv);
			}
*/
			// This is not true. Parameters can have attributes that are used for other purposes. E.g. The local fiber option.
			//		else felog.printf("WARNING: attribute \"%s\" of parameter \"%s\" ignored (line %d)\n", szat, tag.Name(), tag.m_ncurrent_line-1);
		}
	}

	// Set the watch flag since the parameter was read in successfully
	// (This requires that the parameter was declared with a watch variable)
	pp->SetWatchFlag(true);

	return true;
}

//-----------------------------------------------------------------------------
void FEFileSection::ReadAttributes(XMLTag& tag, FECoreBase* pc)
{
	FEParameterList& pl = pc->GetParameterList();

	// process all the other attributes
	for (XMLAtt& att : tag.m_att)
	{
		const char* szatt = att.name();
		const char* szval = att.cvalue();
		if ((att.m_bvisited == false) && szatt && szval)
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
				case FE_PARAM_STD_STRING: param->value<std::string>() = szval; break;
				default:
					throw XMLReader::InvalidAttributeValue(tag, szatt, szval);
				}
			}
			else if (strcmp("name", szatt) == 0) pc->SetName(szval);
			else if (strcmp("id"  , szatt) == 0) pc->SetID(atoi(szval));
			else if (strcmp("lc"  , szatt) == 0) { /* don't do anything. Loadcurves are processed elsewhere. */ }
			else
			{
				if (m_pim->StopOnUnknownAttribute())
				{
					throw XMLReader::InvalidAttribute(tag, szatt);
				}
				else
				{
					fprintf(stderr, "WARNING: Unknown attribute %s in tag %s\n", szatt, tag.Name());
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

			if (prop->IsReference())
			{
				// get the reference. It is either defined by the ref attribute
				// or the value of the tag.
				if (tag.isleaf() == false) throw XMLReader::InvalidValue(tag);

				const char* szref = tag.AttributeValue("ref", true);
				if (szref == nullptr) szref = tag.szvalue();

				const char* sztag = tag.Name();

				FEMesh& mesh = GetFEModel()->GetMesh();

				// This property should reference an existing class
				SUPER_CLASS_ID classID = prop->GetSuperClassID();
/*				if (classID == FEITEMLIST_ID)
				{
					FENodeSet* nodeSet = mesh.FindNodeSet(szref);
					if (nodeSet == nullptr) throw XMLReader::InvalidValue(tag);
					prop->SetProperty(nodeSet);
				}
				else */if (classID == FESURFACE_ID)
				{
					FEModelBuilder* builder = GetBuilder();
					FEFacetSet* facetSet = mesh.FindFacetSet(szref);
					if (facetSet == nullptr) throw XMLReader::InvalidValue(tag);

					FESurface* surface = fecore_alloc(FESurface, GetFEModel());
					GetBuilder()->BuildSurface(*surface, *facetSet);
					mesh.AddSurface(surface);

					prop->SetProperty(surface);
				}
				else throw XMLReader::InvalidTag(tag);
			}
			else
			{
				// see if the property is already allocated
				if ((prop->IsArray() == false) && (prop->get(0)))
				{
					// If so, let's just read the parameters
					FECoreBase* pc = prop->get(0);
					if (tag.isleaf() == false) ReadParameterList(tag, pc);
				}
				else
				{
					const char* sztype = tag.AttributeValue("type", true);

					// If the type attribute is omitted we try the property's default type,
					// otherwise assume the tag's name is the default type
					if (sztype == nullptr)
					{
						if (prop->GetDefaultType()) sztype = prop->GetDefaultType();
						else sztype = tag.Name();
					}

					// HACK for getting passed the old "user" fiber type.
					if (strcmp(sztype, "user") == 0) sztype = "map";

					// HACK for mapping load curves to FEFunction1D
					const char* szlc = tag.AttributeValue("lc", true);
					if (szlc && (tag.m_att.size() == 1) && (prop->GetSuperClassID() == FEFUNCTION1D_ID))
					{
						double v = 1;
						tag.value(v);
						FEPointFunction* f = fecore_alloc(FEPointFunction, GetFEModel()); assert(f);
						prop->SetProperty(f);

						int lc = atoi(szlc) - 1;
						GetBuilder()->MapLoadCurveToFunction(f, lc, v);
					}
					else
					{
						// try to allocate the class
						FECoreBase* pp = fecore_new<FECoreBase>(prop->GetSuperClassID(), sztype, GetFEModel());
						if (pp == nullptr) throw XMLReader::InvalidAttributeValue(tag, "type", sztype);

						prop->SetProperty(pp);

						// read the property data
						if (tag.isleaf() == false)
						{
							ReadParameterList(tag, pp);
						}
						else if (tag.isempty() == false)
						{
							if ((tag.szvalue() != nullptr) && (tag.szvalue()[0] != 0))
							{
								// parse attributes first
								ReadAttributes(tag, pp);

								// There should be a parameter with the same name as the type
								if (ReadParameter(tag, pp->GetParameterList(), sztype, pp) == false)
									throw XMLReader::InvalidValue(tag);
							}
						}
						else
						{
							// we get here if the property was defined with an empty tag.
							// We should still validate it.
							int NP = pp->PropertyClasses();
							for (int i = 0; i < NP; ++i)
							{
								FEProperty* pi = pp->PropertyClass(i);
								bool a = pi->IsRequired();
								bool b = (pi->size() == 0);
								if (a && b)
								{
									std::string name = pp->GetName();
									if (name.empty()) name = prop->GetName();
									throw FEBioImport::MissingProperty(name, pi->GetName());
								}
							}

						}
					}
				}
				return true;
			}
		}
		else
		{
			// backward compatibility hack for older formats (< v 3.0)
/*			if (strcmp(tag.Name(), "fiber") == 0)
			{
				return ReadParameter(tag, pc, "mat_axis", parseAttributes);
			}
			else return false;
*/			return false;
		}
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
			if (ReadParameter(tag, pc) == false)
			{
				if ((m_ith == nullptr) || (m_ith->ProcessTag(tag) == false))
					throw XMLReader::InvalidTag(tag);
			}
			++tag;
		}
		while (!tag.isend());
	}
	else if ((tag.isempty() == false) && (pc->Parameters() > 0))
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
	++tag;
	while (!tag.isend())
	{
		std::map<string, FEFileSection*>::iterator is = find(tag.Name());
		if (is != end()) is->second->Parse(tag);
		else throw XMLReader::InvalidTag(tag);

		++tag;
	};
}

//-----------------------------------------------------------------------------
//! class constructor
FEFileImport::FEFileImport()
{
	m_fp = 0;
	m_szerr[0] = 0;
	m_builder = 0;
	m_nversion = 0;

	m_stopOnUnknownAttribute = false;
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

	delete m_builder;
	m_builder = nullptr;
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
		fprintf(stderr, "FATAL ERROR: %s (line %d)\n", e.what(), tag.m_ncurrent_line);
		return false;
	}
	catch (FEBioImport::MissingProperty e)
	{
		fprintf(stderr, "FATAL ERROR: %s\n\n", e.GetErrorString());
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
FEModel* FEFileImport::GetFEModel()
{
	return &m_builder->GetFEModel();
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
void FEFileImport::SetStopOnUnknownAttribute(bool b)
{
	m_stopOnUnknownAttribute = b;
}

//-----------------------------------------------------------------------------
// throw exception if an unknown attribute is found
bool FEFileImport::StopOnUnknownAttribute() const
{
	return m_stopOnUnknownAttribute;
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
