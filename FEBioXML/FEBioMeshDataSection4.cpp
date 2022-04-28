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
#include "FEBioMeshDataSection.h"
#include "FECore/FEModel.h"
#include "FECore/DOFS.h"
#include <FECore/FEDataGenerator.h>
#include <FECore/FECoreKernel.h>
#include <FECore/FEDataMathGenerator.h>
#include <FECore/FEMaterial.h>
#include <FECore/FEModelParam.h>
#include <FECore/FEDomainMap.h>
#include <FECore/FESurfaceLoad.h>
#include <FECore/FEBodyLoad.h>
#include <FECore/FEPrescribedDOF.h>
#include <FECore/FEMaterialPointProperty.h>
#include <FECore/FEConstDataGenerator.h>
#include <FECore/FEConstValueVec3.h>
#include <sstream>

// defined in FEBioMeshDataSection3.cpp
extern FEDataType str2datatype(const char* szdataType);

//-----------------------------------------------------------------------------
#ifdef WIN32
#define szcmp    _stricmp
#else
#define szcmp    strcmp
#endif

//-----------------------------------------------------------------------------
void FEBioMeshDataSection4::Parse(XMLTag& tag)
{
	// Make sure there is something in this tag
	if (tag.isleaf()) return;

	// loop over all mesh data section
	++tag;
	do
	{
		if      (tag == "NodeData"   ) ParseNodalData  (tag);
		else if (tag == "SurfaceData") ParseSurfaceData(tag);
		else if (tag == "ElementData") ParseElementData(tag);
		else throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());
}

//-----------------------------------------------------------------------------
void FEBioMeshDataSection4::ParseNodalData(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// find the element set
	const char* szset = tag.AttributeValue("node_set");

	// find the element set in the mesh
	FENodeSet* nset = GetBuilder()->FindNodeSet(szset);
	if (nset == nullptr) throw XMLReader::InvalidAttributeValue(tag, "node_set", szset);

	// get the type
	const char* sztype = tag.AttributeValue("type");

	// get the name (required!)
	const char* szname = tag.AttributeValue("name");

	// allocate the data generator
	FENodeDataGenerator* gen = fecore_new<FENodeDataGenerator>(sztype, &fem);
	if (gen == 0) throw XMLReader::InvalidAttributeValue(tag, "type", sztype);

	gen->SetName(szname);

	GetBuilder()->GetFEModel().AddMeshDataGenerator(gen);

	// read the parameters
	ReadParameterList(tag, gen);
}

void FEBioMeshDataSection4::ParseSurfaceData(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// find the element set
	const char* szset = tag.AttributeValue("surface");

	// find the element set in the mesh
	FEFacetSet* surf = mesh.FindFacetSet(szset);
	if (surf == nullptr) throw XMLReader::InvalidAttributeValue(tag, "surface", szset);

	// get the type
	const char* sztype = tag.AttributeValue("type", true);

	// get the name (required!)
	const char* szname = tag.AttributeValue("name");

	// see if there is a generator
	if (sztype)
	{
		FEFaceDataGenerator* gen = fecore_new<FEFaceDataGenerator>(sztype, &fem);
		if (gen == 0) throw XMLReader::InvalidAttributeValue(tag, "type", sztype);

		gen->SetName(szname);

		GetBuilder()->GetFEModel().AddMeshDataGenerator(gen);

		// read the parameters
		ReadParameterList(tag, gen);
	}
	else
	{
		// get the data type
		const char* szdataType = tag.AttributeValue("datatype", true);
		if (szdataType == nullptr) szdataType = "scalar";
		FEDataType dataType = str2datatype(szdataType);
		if (dataType == FEDataType::FE_INVALID_TYPE) throw XMLReader::InvalidAttributeValue(tag, "datatype", szdataType);

		// create the data map
		FESurfaceMap* map = new FESurfaceMap(dataType);
		map->Create(surf);
		map->SetName(szname);

		// add it to the mesh
		mesh.AddDataMap(map);

		// read the data
		ParseSurfaceData(tag, *map);
	}
}

void FEBioMeshDataSection4::ParseElementData(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// find the element set
	const char* szset = tag.AttributeValue("elem_set");

	// find the element set in the mesh
	FEElementSet* elset = mesh.FindElementSet(szset);
	if (elset == nullptr) throw XMLReader::InvalidAttributeValue(tag, "elem_set", szset);

	// get the type
	const char* sztype = tag.AttributeValue("type", true);

	const char* szname = tag.AttributeValue("name", true);

	if (sztype)
	{
		if (strcmp(sztype, "shell thickness") == 0) ParseShellThickness(tag, *elset);
		else
		{
			// get the name or var (required!)
			const char* szname = tag.AttributeValue("name");

			FEElemDataGenerator* gen = fecore_new<FEElemDataGenerator>(sztype, &fem);
			gen->SetElementSet(elset);
			gen->SetName(szname);

			GetBuilder()->GetFEModel().AddMeshDataGenerator(gen);

			// read the parameters
			ReadParameterList(tag, gen);

			// Add it to the list (will be applied after the rest of the model was read in)
			GetBuilder()->AddMeshDataGenerator(gen, nullptr, nullptr);
		}
	}
	else
	{
		if (szname == nullptr) throw XMLReader::InvalidAttributeValue(tag, "name");
		string name = szname;

		// get the data type
		const char* szdataType = tag.AttributeValue("datatype", true);
		if (szdataType == nullptr) szdataType = "scalar";
		FEDataType dataType = str2datatype(szdataType);
		if (dataType == FEDataType::FE_INVALID_TYPE) throw XMLReader::InvalidAttributeValue(tag, "datatype", szdataType);

		// default format
		Storage_Fmt fmt = (((dataType == FE_MAT3D) || (dataType == FE_MAT3DS)) ? Storage_Fmt::FMT_ITEM : Storage_Fmt::FMT_MULT);

		// format overrider?
		const char* szfmt = tag.AttributeValue("format", true);
		if (szfmt)
		{
			if (szcmp(szfmt, "MAT_POINTS") == 0) fmt = Storage_Fmt::FMT_MATPOINTS;
			else if (szcmp(szfmt, "ITEM") == 0) fmt = Storage_Fmt::FMT_ITEM;
			else throw XMLReader::InvalidAttributeValue(tag, "format", szfmt);
		}

		// create the data map
		FEDomainMap* map = new FEDomainMap(dataType, fmt);
		map->Create(elset);
		map->SetName(name);

		if (tag.isleaf())
		{
			if (dataType == FE_DOUBLE)
			{
				double v = 0.0;
				tag.value(v);
				map->set(v);
			}
			else throw XMLReader::InvalidValue(tag);
		}
		else
		{
			// read the data
			ParseElementData(tag, *map);
		}

		// see if this map already exsits 
		FEDomainMap* oldMap = dynamic_cast<FEDomainMap*>(mesh.FindDataMap(name));
		if (oldMap)
		{
			oldMap->Merge(*map);
			delete map;
		}
		else
		{
			map->SetName(name);
			mesh.AddDataMap(map);
		}
	}
}

//-----------------------------------------------------------------------------
void FEBioMeshDataSection4::ParseSurfaceData(XMLTag& tag, FESurfaceMap& map)
{
	const FEFacetSet* set = map.GetFacetSet();
	if (set == nullptr) throw XMLReader::InvalidTag(tag);

	// get the total nr of elements
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	int nelems = set->Faces();

	FEDataType dataType = map.DataType();
	int dataSize = map.DataSize();
	int m = map.MaxNodes();
	double data[3 * FEElement::MAX_NODES]; // make sure this array is large enough to store any data map type (current 3 for FE_VEC3D)

	++tag;
	do
	{
		// get the local element number
		const char* szlid = tag.AttributeValue("lid");
		int n = atoi(szlid) - 1;

		// make sure the number is valid
		if ((n < 0) || (n >= nelems)) throw XMLReader::InvalidAttributeValue(tag, "lid", szlid);

		int nread = tag.value(data, m * dataSize);
		if (nread == dataSize)
		{
			switch (dataType)
			{
			case FE_DOUBLE:	map.setValue(n, data[0]); break;
			case FE_VEC2D:	map.setValue(n, vec2d(data[0], data[1])); break;
			case FE_VEC3D:	map.setValue(n, vec3d(data[0], data[1], data[2])); break;
			default:
				assert(false);
			}
		}
		else if (nread == m * dataSize)
		{
			double* pd = data;
			for (int i = 0; i < m; ++i, pd += dataSize)
			{
				switch (dataType)
				{
				case FE_DOUBLE:	map.setValue(n, i, pd[0]); break;
				case FE_VEC2D:	map.setValue(n, i, vec2d(pd[0], pd[1])); break;
				case FE_VEC3D:	map.setValue(n, i, vec3d(pd[0], pd[1], pd[2])); break;
				default:
					assert(false);
				}
			}
		}
		else throw XMLReader::InvalidValue(tag);
		++tag;
	} while (!tag.isend());
}

//-----------------------------------------------------------------------------
void FEBioMeshDataSection4::ParseShellThickness(XMLTag& tag, FEElementSet& set)
{
	if (tag.isleaf())
	{
		FEMesh& mesh = GetFEModel()->GetMesh();
		double h[FEElement::MAX_NODES];
		int nval = tag.value(h, FEElement::MAX_NODES);

		for (int i = 0; i < set.Elements(); ++i)
		{
			FEShellElement* pel = dynamic_cast<FEShellElement*>(&set.Element(i));
			if (pel == 0) throw XMLReader::InvalidValue(tag);

			if (pel->Nodes() != nval) throw XMLReader::InvalidValue(tag);
			for (int j = 0; j < nval; ++j) pel->m_h0[j] = h[j];
		}
	}
	else
	{
		vector<ELEMENT_DATA> data;
		ParseElementData(tag, set, data, FEElement::MAX_NODES);
		for (int i = 0; i < (int)data.size(); ++i)
		{
			ELEMENT_DATA& di = data[i];
			if (di.nval > 0)
			{
				FEElement& el = set.Element(i);

				if (el.Class() != FE_ELEM_SHELL) throw XMLReader::InvalidTag(tag);
				FEShellElement& shell = static_cast<FEShellElement&> (el);

				int ne = shell.Nodes();
				if (ne != di.nval) throw XMLReader::InvalidTag(tag);
				for (int j = 0; j < ne; ++j) shell.m_h0[j] = di.val[j];
			}
		}
	}
}

//-----------------------------------------------------------------------------
void FEBioMeshDataSection4::ParseElementData(XMLTag& tag, FEElementSet& set, vector<ELEMENT_DATA>& values, int nvalues)
{
	// get the total nr of elements
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	int nelems = set.Elements();

	// resize the array
	values.resize(nelems);
	for (int i = 0; i < nelems; ++i) values[i].nval = 0;

	++tag;
	do
	{
		// get the local element number
		const char* szlid = tag.AttributeValue("lid");
		int n = atoi(szlid) - 1;

		// make sure the number is valid
		if ((n < 0) || (n >= nelems)) throw XMLReader::InvalidAttributeValue(tag, "lid", szlid);

		ELEMENT_DATA& data = values[n];
		data.nval = tag.value(data.val, nvalues);
		++tag;
	} while (!tag.isend());
}

//-----------------------------------------------------------------------------
void FEBioMeshDataSection4::ParseElementData(XMLTag& tag, FEDomainMap& map)
{
	const FEElementSet* set = map.GetElementSet();
	if (set == nullptr) throw XMLReader::InvalidTag(tag);

	// get the total nr of elements
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	int nelems = set->Elements();

	FEDataType dataType = map.DataType();
	int dataSize = map.DataSize();
	int m = map.MaxNodes();
	double data[3 * FEElement::MAX_NODES]; // make sure this array is large enough to store any data map type (current 3 for FE_VEC3D)

	// TODO: For vec3d values, I sometimes need to normalize the vectors (e.g. for fibers). How can I do this?

	int ncount = 0;
	++tag;
	do
	{
		// get the local element number
		const char* szlid = tag.AttributeValue("lid");
		int n = atoi(szlid) - 1;

		// make sure the number is valid
		if ((n < 0) || (n >= nelems)) throw XMLReader::InvalidAttributeValue(tag, "lid", szlid);

		int nread = tag.value(data, m * dataSize);
		if (nread == dataSize)
		{
			double* v = data;
			switch (dataType)
			{
			case FE_DOUBLE:	map.setValue(n, v[0]); break;
			case FE_VEC2D:	map.setValue(n, vec2d(v[0], v[1])); break;
			case FE_VEC3D:	map.setValue(n, vec3d(v[0], v[1], v[2])); break;
			case FE_MAT3D: map.setValue(n, mat3d(v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8])); break;
			case FE_MAT3DS: map.setValue(n, mat3ds(v[0], v[1], v[2], v[3], v[4], v[5])); break;
			default:
				assert(false);
			}
		}
		else if (nread == m * dataSize)
		{
			double* v = data;
			for (int i = 0; i < m; ++i, v += dataSize)
			{
				switch (dataType)
				{
				case FE_DOUBLE:	map.setValue(n, i, v[0]); break;
				case FE_VEC2D:	map.setValue(n, i, vec2d(v[0], v[1])); break;
				case FE_VEC3D:	map.setValue(n, i, vec3d(v[0], v[1], v[2])); break;
				default:
					assert(false);
				}
			}
		}
		else throw XMLReader::InvalidValue(tag);
		++tag;

		ncount++;
	} while (!tag.isend());

	if (ncount != nelems) throw FEBioImport::MeshDataError();
}