/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2020 University of Utah, The Trustees of Columbia University in 
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
#include "FEBioMeshDataSection.h"
#include "FECore/FEModel.h"
#include "FECore/DOFS.h"
#include <FECore/FEDataGenerator.h>
#include <FECore/FECoreKernel.h>
#include <FECore/FEDataMathGenerator.h>
#include <FECore/FEMaterial.h>
#include <FECore/FEModelParam.h>
#include <FECore/FEDomainMap.h>
#include <sstream>

//-----------------------------------------------------------------------------
void FEBioMeshDataSection::Parse(XMLTag& tag)
{
	// Make sure there is something in this tag
	if (tag.isleaf()) return;

	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// get the total nr of elements
	int nelems = mesh.Elements();

	//make sure we've read the element section
	if (nelems == 0) throw XMLReader::InvalidTag(tag);

	// get the largest ID
	int max_id = 0;
	for (int i = 0; i<mesh.Domains(); ++i)
	{
		FEDomain& d = mesh.Domain(i);
		for (int j = 0; j<d.Elements(); ++j)
		{
			FEElement& el = d.ElementRef(j);
			if (el.GetID() > max_id) max_id = el.GetID();
		}
	}

	// create the pelem array
	m_pelem.assign(max_id, static_cast<FEElement*>(0));
	for (int nd = 0; nd<mesh.Domains(); ++nd)
	{
		FEDomain& d = mesh.Domain(nd);
		for (int i = 0; i<d.Elements(); ++i)
		{
			FEElement& el = d.ElementRef(i);
			assert(m_pelem[el.GetID() - 1] == 0);
			m_pelem[el.GetID() - 1] = &el;
		}
	}

	FEBioImport* feb = GetFEBioImport();

	// loop over all mesh data
	++tag;
	do
	{
		if (tag == "ElementData")
		{
			const char* szset = tag.AttributeValue("elem_set");
			FEElementSet* part = mesh.FindElementSet(szset);
			if (part == 0) throw XMLReader::InvalidAttributeValue(tag, "elem_set", szset);

			// see if the data will be generated or tabulated
			const char* szgen = tag.AttributeValue("generator", true);

			if (szgen == 0)
			{
				// data is tabulated and mapped directly to a variable.
				const char* szvar = tag.AttributeValue("var", true);
				if (szvar)
				{
					if (strcmp(szvar, "shell thickness") == 0) ParseShellThickness(tag, *part);
					else if (strcmp(szvar, "fiber") == 0) ParseMaterialFibers(tag, *part);
					else if (strcmp(szvar, "mat_axis") == 0) ParseMaterialAxes(tag, *part);
					else if (strstr(szvar, ".fiber") != 0) ParseMaterialFiberProperty(tag, *part);
					else if (strstr(szvar, ".mat_axis") != 0) ParseMaterialAxesProperty(tag, *part);
					else ParseMaterialData(tag, *part, szvar);
				}
				else
				{
					const char* szname = tag.AttributeValue("name");
					FEDomainMap* map = new FEDomainMap(FE_DOUBLE);
					map->Create(part);
					map->SetName(szname);

					// add it to the mesh
					mesh.AddDataMap(map);
				}
			}
			else
			{
				// data will be generated
				FEDataGenerator* gen = fecore_new<FEDataGenerator>(szgen, &fem);
				if (gen == 0) throw XMLReader::InvalidAttributeValue(tag, "generator", szgen);

				// get the variable or name
				const char* szvar = tag.AttributeValue("var", true);
				const char* szname = (szvar == nullptr ? tag.AttributeValue("name") : nullptr);

				// make sure the parameter is valid
				FEParamDouble* pp = nullptr;
				if (szvar)
				{
					// find the variable
					ParamString paramName(szvar);
					FEParam* pv = fem.FindParameter(paramName);
					if (pv == nullptr) throw XMLReader::InvalidAttributeValue(tag, "var", szvar);

					// make sure it's a mapped parameter
					if (pv->type() != FE_PARAM_DOUBLE_MAPPED)throw XMLReader::InvalidAttributeValue(tag, "var", szvar);

					// if it's an array parameter, get the right index
					if (pv->dim() > 1)
					{
						ParamString l = paramName.last();
						int m = l.Index();
						assert(m >= 0);
						pp = &(pv->value<FEParamDouble>(m));
					}
					else pp = &(pv->value<FEParamDouble>());
				}

				// read the parameters of the generator
				ReadParameterList(tag, gen);

				// give the generator a chance to validate itself
				if (gen->Init() == false) throw FEBioImport::DataGeneratorError();

				// create a new domain map (only scalars for now!)
				FEDomainMap* map = new FEDomainMap(FE_DOUBLE);
				map->Create(part);

				// generate the data
				if (gen->Generate(*map) == false)
				{
					delete map;
					throw FEBioImport::DataGeneratorError();
				}

				// either map it directly to a variable or just store it
				if (szvar)
				{
					FEParamDouble& p = *pp;

					// see if this map is already defined
					FEDomainMap* oldMap = dynamic_cast<FEDomainMap*>(mesh.FindDataMap(szvar));
					if (oldMap)
					{
						// it is, so merge it
						oldMap->Merge(*map);

						// we can now delete this map
						delete map;
					}
					else
					{
						// nope, so add it
						map->SetName(szvar);
						mesh.AddDataMap(map);

						// apply the map
						FEMappedValue* val = fecore_alloc(FEMappedValue, &fem);
						val->setDataMap(map);
						p.setValuator(val);
					}
				}
				else
				{
					// see if this map is already defined
					FEDomainMap* oldMap = dynamic_cast<FEDomainMap*>(mesh.FindDataMap(szname));
					if (oldMap)
					{
						// it is, so merge it
						oldMap->Merge(*map);

						// we can now delete this map
						delete map;
					}
					else
					{
						// set the name
						map->SetName(szname);

						// add it to the mesh
						mesh.AddDataMap(map);
					}
				}
			}
		}
		else if (tag == "SurfaceData")
		{
			const char* szsurf = tag.AttributeValue("surface");
			FEFacetSet* psurf = mesh.FindFacetSet(szsurf);
			if (psurf == 0) throw XMLReader::InvalidAttributeValue(tag, "surface", szsurf);

			const char* sztype = tag.AttributeValue("data_type", true);
			if (sztype == 0) sztype = "scalar";

			FEDataType dataType;
			if      (strcmp(sztype, "scalar") == 0) dataType = FE_DOUBLE;
			else if (strcmp(sztype, "vec2"  ) == 0) dataType = FE_VEC2D;
			else if (strcmp(sztype, "vec3"  ) == 0) dataType = FE_VEC3D;
			else throw XMLReader::InvalidAttributeValue(tag, "data_type", sztype);

			const char* szname = tag.AttributeValue("name");
			FESurfaceMap* pdata = new FESurfaceMap(dataType);
			mesh.AddDataMap(pdata);

			pdata->SetName(szname);
			pdata->Create(psurf);

			const char* szgen = tag.AttributeValue("generator", true);
			if (szgen)
			{
				if (dataType != FE_DOUBLE) throw XMLReader::InvalidAttributeValue(tag, "generator", szgen);

				FEDataGenerator* gen = fecore_new<FEDataGenerator>(szgen, &fem);
				if (gen == nullptr) throw XMLReader::InvalidAttributeValue(tag, "generator", szgen);

				ReadParameterList(tag, gen);

				gen->Generate(*pdata);
			}
			else
				ParseDataArray(tag, *pdata, "face");
		}
		else if (tag == "EdgeData")
		{
			const char* szedge = tag.AttributeValue("edge");
			FESegmentSet* pset = mesh.FindSegmentSet(szedge);
			if (pset == 0) throw XMLReader::InvalidAttributeValue(tag, "edge", szedge);

			const char* sztype = tag.AttributeValue("data_type", true);
			if (sztype == 0) sztype = "scalar";

			int dataType = -1;
			if (strcmp(sztype, "scalar") == 0) dataType = FE_DOUBLE;
			else if (strcmp(sztype, "vec2") == 0) dataType = FE_VEC2D;
			else if (strcmp(sztype, "vec3") == 0) dataType = FE_VEC3D;
			if (dataType == -1) throw XMLReader::InvalidAttributeValue(tag, "data_type", sztype);
			/*
			const char* szname = tag.AttributeValue("name");
			FEDataArray* pdata = new FEDataArray(dataType);
			fem.AddDataArray(szname, pdata);

			pdata->resize(pset->Segments());
			ParseDataArray(tag, *pdata, "edge");
			*/
		}
		else if (tag == "NodeData")
		{
			const char* szset = tag.AttributeValue("node_set");
			FENodeSet* nodeSet = mesh.FindNodeSet(szset);
			if (nodeSet == 0) throw XMLReader::InvalidAttributeValue(tag, "node_set", szset);

			const char* sztype = tag.AttributeValue("data_type", true);
			if (sztype == 0) sztype = "scalar";

			FEDataType dataType;
			if      (strcmp(sztype, "scalar") == 0) dataType = FE_DOUBLE;
			else if (strcmp(sztype, "vec2"  ) == 0) dataType = FE_VEC2D;
			else if (strcmp(sztype, "vec3"  ) == 0) dataType = FE_VEC3D;
			else throw XMLReader::InvalidAttributeValue(tag, "data_type", sztype);

			const char* szname = tag.AttributeValue("name");

			FENodeDataMap* pdata = new FENodeDataMap(dataType);
			pdata->SetName(szname);
			mesh.AddDataMap(pdata);

			pdata->Create(nodeSet);

			const char* szgen = tag.AttributeValue("generator", true);
			if (szgen)
			{
				if (dataType != FE_DOUBLE) throw XMLReader::InvalidAttributeValue(tag, "generator", szgen);

				++tag;
				do
				{
					if (tag == "math")
					{
						FEDataMathGenerator gen(&fem);
						gen.setExpression(tag.szvalue());
						if (gen.Init() == false)  throw XMLReader::InvalidValue(tag);
						if (gen.Generate(*pdata) == false) throw XMLReader::InvalidValue(tag);
					}
					else throw XMLReader::InvalidTag(tag);
					++tag;
				} while (!tag.isend());
			}
			else
			{
				ParseDataArray(tag, *pdata, "node");
			}
		}
		else throw XMLReader::InvalidTag(tag);
		++tag;
	} while (!tag.isend());
}

//-----------------------------------------------------------------------------
void FEBioMeshDataSection::ParseShellThickness(XMLTag& tag, FEElementSet& set)
{
	if (tag.isleaf())
	{
		FEMesh& mesh = GetFEModel()->GetMesh();
		double h[FEElement::MAX_NODES];
		int nval = tag.value(h, FEElement::MAX_NODES);

		for (int i = 0; i<set.Elements(); ++i)
		{
			FEShellElement* pel = dynamic_cast<FEShellElement*>(m_pelem[set[i] - 1]);
			if (pel == 0) throw XMLReader::InvalidValue(tag);

			if (pel->Nodes() != nval) throw XMLReader::InvalidValue(tag);
			for (int j = 0; j<nval; ++j) pel->m_h0[j] = h[j];
		}
	}
	else
	{
		vector<ELEMENT_DATA> data;
		ParseElementData(tag, set, data, FEElement::MAX_NODES);
		for (int i = 0; i<(int)data.size(); ++i)
		{
			ELEMENT_DATA& di = data[i];
			if (di.nval > 0)
			{
				FEElement& el = *m_pelem[set[i] - 1];

				if (el.Class() != FE_ELEM_SHELL) throw XMLReader::InvalidTag(tag);
				FEShellElement& shell = static_cast<FEShellElement&> (el);

				int ne = shell.Nodes();
				if (ne != di.nval) throw XMLReader::InvalidTag(tag);
				for (int j = 0; j<ne; ++j) shell.m_h0[j] = di.val[j];
			}
		}
	}
}

//-----------------------------------------------------------------------------
// Defined in FEBioGeometrySection.cpp
void set_element_fiber(FEElement& el, const vec3d& v, int ncomp);
void set_element_mat_axis(FEElement& el, const vec3d& v1, const vec3d& v2, int ncomp);

//-----------------------------------------------------------------------------
void FEBioMeshDataSection::ParseMaterialFibers(XMLTag& tag, FEElementSet& set)
{
	// find the domain with the same name
	string name = set.GetName();

	FEMesh* mesh = const_cast<FEMesh*>(set.GetMesh());
	FEDomain* dom = mesh->FindDomain(name);
	if (dom == nullptr) throw XMLReader::InvalidAttributeValue(tag, "elem_set", name.c_str());

	// get the material
	FEMaterial* mat = dom->GetMaterial();
	if (mat == nullptr) throw XMLReader::InvalidAttributeValue(tag, "elem_set", name.c_str());

	// get the fiber parameter
	ParamString s("fiber");
	FEParam* param = mat->FindParameter(s);
	if (param == nullptr) throw XMLReader::InvalidAttributeValue(tag, "elem_set", name.c_str());
	if (param->type() != FE_PARAM_VEC3D_MAPPED) throw XMLReader::InvalidAttributeValue(tag, "elem_set", name.c_str());

	// get the parameter
	FEParamVec3& p = param->value<FEParamVec3>();

	// create a domain map
	FEDomainMap* map = new FEDomainMap(FE_VEC3D, FMT_ITEM);
	map->Create(&set);
	FEMappedValueVec3* val = fecore_new<FEMappedValueVec3>("map", GetFEModel());
	val->setDataMap(map);
	p.setValuator(val);

	vector<ELEMENT_DATA> data;
	ParseElementData(tag, set, data, 3);
	for (int i = 0; i<(int)data.size(); ++i)
	{
		ELEMENT_DATA& di = data[i];
		if (di.nval > 0)
		{
			FEElement& el = *m_pelem[set[i] - 1];

			if (di.nval != 3) throw XMLReader::InvalidTag(tag);
			vec3d v(di.val[0], di.val[1], di.val[2]);
			v.unit();
			map->set<vec3d>(i, v);
		}
	}
}

//-----------------------------------------------------------------------------
void FEBioMeshDataSection::ParseMaterialFiberProperty(XMLTag& tag, FEElementSet& set)
{
	const char* szvar = tag.AttributeValue("var");
	char szbuf[256] = { 0 };
	strcpy(szbuf, szvar);
	char* ch = strstr(szbuf, ".fiber");
	if (ch == 0) return;
	*ch = 0;
	ch = strrchr(szbuf, ']');
	if (ch == 0) return;
	*ch = 0;
	ch = strchr(szbuf, '[');
	if (ch == 0) return;
	*ch++ = 0;
	int nindex = atoi(ch);

	vector<ELEMENT_DATA> data;
	ParseElementData(tag, set, data, 3);
	for (int i = 0; i<(int)data.size(); ++i)
	{
		ELEMENT_DATA& di = data[i];
		if (di.nval > 0)
		{
			FEElement& el = *m_pelem[set[i] - 1];

			if (di.nval != 3) throw XMLReader::InvalidTag(tag);
			vec3d v(di.val[0], di.val[1], di.val[2]);

			set_element_fiber(el, v, nindex);
		}
	}
}

//-----------------------------------------------------------------------------
void FEBioMeshDataSection::ParseMaterialAxesProperty(XMLTag& tag, FEElementSet& set)
{
	const char* szvar = tag.AttributeValue("var");
	char szbuf[256] = { 0 };
	strcpy(szbuf, szvar);
	char* ch = strstr(szbuf, ".mat_axis");
	if (ch == 0) return;
	*ch = 0;
	ch = strrchr(szbuf, ']');
	if (ch == 0) return;
	*ch = 0;
	ch = strchr(szbuf, '[');
	if (ch == 0) return;
	*ch++ = 0;
	int nindex = atoi(ch);

	// get the total nr of elements
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	++tag;
	do
	{
		if (tag == "elem")
		{
			// get the local element number
			const char* szlid = tag.AttributeValue("lid");
			int lid = atoi(szlid) - 1;

			// make sure the number is valid
			if ((lid<0) || (lid >= set.Elements())) throw XMLReader::InvalidAttributeValue(tag, "lid", szlid);

			// get the element
			FEElement* el = mesh.FindElementFromID(set[lid]);
			if (el == 0) throw XMLReader::InvalidAttributeValue(tag, "lid", szlid);

			// read parameters
			double a[3] = { 0 };
			double d[3] = { 0 };
			++tag;
			do
			{
				if (tag == "a") tag.value(a, 3);
				else if (tag == "d") tag.value(d, 3);
				else throw XMLReader::InvalidTag(tag);
				++tag;
			} while (!tag.isend());

			vec3d v1(a[0], a[1], a[2]);
			vec3d v2(d[0], d[1], d[2]);
			set_element_mat_axis(*el, v1, v2, nindex);
		}
		else throw XMLReader::InvalidTag(tag);
		++tag;
	} while (!tag.isend());
}

//-----------------------------------------------------------------------------
void FEBioMeshDataSection::ParseMaterialAxes(XMLTag& tag, FEElementSet& set)
{
	// get the set's name
	const char* szname = set.GetName().c_str();

	// get the domain list
	const FEDomainList& domList = set.GetDomainList();

	// NOTE For now, this only works for single domains
	if (domList.Domains() != 1) throw XMLReader::InvalidAttributeValue(tag, "elem_set", szname);

	FEDomain* dom = const_cast<FEDomain*>(domList.GetDomain(0));
	if (dom == nullptr) throw XMLReader::InvalidAttributeValue(tag, "elem_set", szname);

	FEMesh* mesh = dom->GetMesh();

	// get the material
	FEMaterial* mat = dom->GetMaterial();
	if (mat == nullptr) throw XMLReader::InvalidAttributeValue(tag, "elem_set", szname);

	// get the mat_axis parameter
	ParamString s("mat_axis");
	FEParam* param = mat->FindParameter(s);
	if (param == nullptr) throw XMLReader::InvalidAttributeValue(tag, "elem_set", szname);
	if (param->type() != FE_PARAM_MAT3D_MAPPED) throw XMLReader::InvalidAttributeValue(tag, "elem_set", szname);

	// get the parameter
	FEParamMat3d& p = param->value<FEParamMat3d>();

	// create the map's name: material_name.mat_axis
	stringstream ss;
	ss << "material" << mat->GetID() << ".mat_axis";
	string mapName = ss.str();

	// create a domain map
	FEDomainMap* map = new FEDomainMap(FE_MAT3D, FMT_ITEM);
	map->SetName(mapName);
	map->Create(&set);

	++tag;
	do
	{
		if (tag == "elem")
		{
			// get the local element number
			const char* szlid = tag.AttributeValue("lid");
			int lid = atoi(szlid) - 1;

			// make sure the number is valid
			if ((lid<0) || (lid >= set.Elements())) throw XMLReader::InvalidAttributeValue(tag, "lid", szlid);

			// get the element
			FEElement* el = mesh->FindElementFromID(set[lid]);
			if (el == 0) throw XMLReader::InvalidAttributeValue(tag, "lid", szlid);

			// read parameters
			double a[3] = { 0 };
			double d[3] = { 0 };
			++tag;
			do
			{
				if (tag == "a") tag.value(a, 3);
				else if (tag == "d") tag.value(d, 3);
				else throw XMLReader::InvalidTag(tag);
				++tag;
			} while (!tag.isend());

			vec3d v1(a[0], a[1], a[2]);
			vec3d v2(d[0], d[1], d[2]);

			vec3d e1(v1);
			vec3d e3 = v1^v2;
			vec3d e2 = e3^e1;

			// normalize
			e1.unit();
			e2.unit();
			e3.unit();

			// set the value
			mat3d Q(e1, e2, e3);
			map->set<mat3d>(lid, Q);
		}
		else throw XMLReader::InvalidTag(tag);
		++tag;
	} while (!tag.isend());

	// see if this map already exists
	FEDomainMap* oldMap = dynamic_cast<FEDomainMap*>(mesh->FindDataMap(mapName));
	if (oldMap)
	{
		// It does, so merge it
		oldMap->Merge(*map);
		delete map;
	}
	else
	{
		// It does not, so add it
		FEMappedValueMat3d* val = fecore_alloc(FEMappedValueMat3d, GetFEModel());
		val->setDataMap(map);
		p.setValuator(val);
		mesh->AddDataMap(map);
	}
}

//-----------------------------------------------------------------------------
void FEBioMeshDataSection::ParseMaterialData(XMLTag& tag, FEElementSet& set, const string& pname)
{
	// get the (optional) scale factor
	double scale = 1.0;
	const char* szscale = tag.AttributeValue("scale", true);
	if (szscale) scale = atof(szscale);

	// find the parameter
	FEModel& fem = *GetFEModel();
	ParamString PS(pname.c_str());
	FEParam* p = fem.FindParameter(PS);
	if (p == nullptr)
	{
		printf("Can't find parameter %s\n", pname.c_str());
		throw XMLReader::InvalidAttributeValue(tag, "var", pname.c_str());
	}

	if (p->type() == FE_PARAM_DOUBLE_MAPPED)
	{
		vector<ELEMENT_DATA> data;
		ParseElementData(tag, set, data, 1);

		FEParamDouble& param = p->value<FEParamDouble>();
		param.SetItemList(&set);

		FEMappedValue* val = fecore_alloc(FEMappedValue, &fem);
		if (val == nullptr)
		{
			printf("Something went horribly wrong.");
			return;
		}

		FEDomainMap* map = new FEDomainMap(FEDataType::FE_DOUBLE, FMT_ITEM);
		map->Create(&set);
		val->setDataMap(map);

		for (int i = 0; i < data.size(); ++i) map->set<double>(i, data[i].val[0] * scale);

		param.setValuator(val);
	}
	else if (p->type() == FE_PARAM_MAT3DS_MAPPED)
	{
		vector<ELEMENT_DATA> data;
		ParseElementData(tag, set, data, 6);

		FEParamMat3ds& param = p->value<FEParamMat3ds>();
		param.SetItemList(&set);

		FEMappedValueMat3ds* val = fecore_alloc(FEMappedValueMat3ds, &fem);
		if (val == nullptr)
		{
			printf("Something went horribly wrong.");
			return;
		}

		FEDomainMap* map = new FEDomainMap(FEDataType::FE_MAT3DS, FMT_ITEM);
		map->Create(&set);
		val->setDataMap(map);

		for (int i = 0; i < data.size(); ++i)
		{
			double* d = data[i].val;
			mat3ds m(d[0], d[1], d[2], d[3], d[4], d[5]);
			map->set<mat3ds>(i, m * scale);
		}

		param.setValuator(val);
	}
	else
	{
		printf("A mesh data map cannot be assigned to this parameter.");
		throw XMLReader::InvalidAttribute(tag, "var");
	}
}

//-----------------------------------------------------------------------------
void FEBioMeshDataSection::ParseElementData(XMLTag& tag, FEElementSet& set, vector<ELEMENT_DATA>& values, int nvalues)
{
	// get the total nr of elements
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	int nelems = set.Elements();

	// resize the array
	values.resize(nelems);
	for (int i = 0; i<nelems; ++i) values[i].nval = 0;

	++tag;
	do
	{
		if (tag == "elem")
		{
			// get the local element number
			const char* szlid = tag.AttributeValue("lid");
			int n = atoi(szlid) - 1;

			// make sure the number is valid
			if ((n<0) || (n >= nelems)) throw XMLReader::InvalidAttributeValue(tag, "lid", szlid);

			ELEMENT_DATA& data = values[n];
			data.nval = tag.value(data.val, nvalues);
		}
		else throw XMLReader::InvalidTag(tag);
		++tag;
	} while (!tag.isend());
}

//-----------------------------------------------------------------------------
void FEBioMeshDataSection::ParseDataArray(XMLTag& tag, FEDataArray& map, const char* sztag)
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

				map.setValue(nid - 1, v);
			}
			else throw XMLReader::InvalidTag(tag);
			++tag;
		} while (!tag.isend());
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

				map.setValue(nid - 1, vec3d(v[0], v[1], v[2]));
			}
			else throw XMLReader::InvalidTag(tag);
			++tag;
		} while (!tag.isend());
	}
}
