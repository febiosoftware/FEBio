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

	// make sure the MeshDomain section was processed. 
	FEMesh& mesh = GetFEModel()->GetMesh();
	if (mesh.Domains() == 0)
	{
		throw FEFileException("MeshData must appear after MeshDomain section.");
	}

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
	const char* sztype = tag.AttributeValue("type", true);

	// get the name (required!)
	const char* szname = tag.AttributeValue("name");

	// see if there is a generator
	if (sztype)
	{
		// allocate the data generator
		FENodeDataGenerator* gen = dynamic_cast<FENodeDataGenerator*>(fecore_new<FEMeshDataGenerator>(sztype, &fem));
		if (gen == 0) throw XMLReader::InvalidAttributeValue(tag, "type", sztype);

		gen->SetName(szname);
		gen->SetNodeSet(nset);

		GetBuilder()->GetFEModel().AddMeshDataGenerator(gen);

		// read the parameters
		ReadParameterList(tag, gen);

		// Add it to the list (will be applied after the rest of the model was read in)
		GetBuilder()->AddMeshDataGenerator(gen, nullptr, nullptr);
	}
	else
	{
		// get the data type
		const char* szdataType = tag.AttributeValue("data_type", true);
		if (szdataType == nullptr) szdataType = "scalar";
		FEDataType dataType = str2datatype(szdataType);
		if (dataType == FEDataType::FE_INVALID_TYPE) throw XMLReader::InvalidAttributeValue(tag, "data_type", szdataType);

		// create the data map
		FENodeDataMap* map = new FENodeDataMap(dataType);
		map->Create(nset);
		map->SetName(szname);

		// add it to the mesh
		mesh.AddDataMap(map);

		// read the data
		ParseNodeData(tag, *map);
	}
}

//-----------------------------------------------------------------------------
void FEBioMeshDataSection4::ParseNodeData(XMLTag& tag, FENodeDataMap& map)
{
	// get the total nr of nodes
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	int nodes = map.DataCount();

	FEDataType dataType = map.DataType();
	int dataSize = map.DataSize();
	double data[3]; // make sure this array is large enough to store any data map type (current 3 for FE_VEC3D)

	++tag;
	do
	{
		// get the local element number
		const char* szlid = tag.AttributeValue("lid");
		int n = atoi(szlid) - 1;

		// make sure the number is valid
		if ((n < 0) || (n >= nodes)) throw XMLReader::InvalidAttributeValue(tag, "lid", szlid);

		int nread = tag.value(data, dataSize);
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
		else throw XMLReader::InvalidValue(tag);
		++tag;
	}
	while (!tag.isend());
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
		FEFaceDataGenerator* gen = dynamic_cast<FEFaceDataGenerator*>(fecore_new<FEMeshDataGenerator>(sztype, &fem));
		if (gen == 0) throw XMLReader::InvalidAttributeValue(tag, "type", sztype);

		gen->SetFacetSet(surf);
		gen->SetName(szname);

		GetBuilder()->GetFEModel().AddMeshDataGenerator(gen);

		// read the parameters
		ReadParameterList(tag, gen);

		// Add it to the list (will be applied after the rest of the model was read in)
		GetBuilder()->AddMeshDataGenerator(gen, nullptr, nullptr);
	}
	else
	{
		// get the data type
		const char* szdataType = tag.AttributeValue("data_type", true);
		if (szdataType == nullptr) szdataType = "scalar";
		FEDataType dataType = str2datatype(szdataType);
		if (dataType == FEDataType::FE_INVALID_TYPE) throw XMLReader::InvalidAttributeValue(tag, "data_type", szdataType);

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

	if (sztype)
	{
		if      (strcmp(sztype, "shell thickness") == 0) ParseShellThickness(tag, *elset);
		else if (strcmp(sztype, "mat_axis"       ) == 0) ParseMaterialAxes  (tag, *elset);
		else if (strcmp(sztype, "fiber"          ) == 0) ParseMaterialFibers(tag, *elset);
		else
		{
			// get the name (required!)
			const char* szname = tag.AttributeValue("name");

			// allocate generator
			FEElemDataGenerator* gen = dynamic_cast<FEElemDataGenerator*>(fecore_new<FEMeshDataGenerator>(sztype, &fem));
			if (gen == nullptr) throw XMLReader::InvalidAttributeValue(tag, "type", sztype);

			// set the name and element set
			gen->SetElementSet(elset);
			gen->SetName(szname);

			// add it to the model
			GetBuilder()->GetFEModel().AddMeshDataGenerator(gen);

			// read the parameters
			ReadParameterList(tag, gen);

			// Add it to the list (will be applied after the rest of the model was read in)
			GetBuilder()->AddMeshDataGenerator(gen, nullptr, nullptr);
		}
	}
	else
	{
		// get the name (required!)
		const char* szname = tag.AttributeValue("name");

		if (szname == nullptr) throw XMLReader::InvalidAttributeValue(tag, "name");
		string name = szname;

		// get the data type
		const char* szdataType = tag.AttributeValue("data_type", true);
		if (szdataType == nullptr) szdataType = "scalar";
		FEDataType dataType = str2datatype(szdataType);
		if (dataType == FEDataType::FE_INVALID_TYPE) throw XMLReader::InvalidAttributeValue(tag, "data_type", szdataType);

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

//-----------------------------------------------------------------------------
void FEBioMeshDataSection4::ParseMaterialFibers(XMLTag& tag, FEElementSet& set)
{
	// find the domain with the same name
	string name = set.GetName();

	FEMesh* mesh = const_cast<FEMesh*>(set.GetMesh());
	FEDomain* dom = mesh->FindDomain(name);
	if (dom == nullptr) throw XMLReader::InvalidAttributeValue(tag, "elem_set", name.c_str());

	// get the material
	FEMaterial* mat = dom->GetMaterial();
	if (mat == nullptr) throw XMLReader::InvalidAttributeValue(tag, "elem_set", name.c_str());

	// get the fiber property
	FEProperty* fiber = mat->FindProperty("fiber");
	if (fiber == nullptr) throw XMLReader::InvalidAttributeValue(tag, "elem_set", name.c_str());
	if (fiber->GetSuperClassID() != FEVEC3DVALUATOR_ID) throw XMLReader::InvalidAttributeValue(tag, "elem_set", name.c_str());

	// create a domain map
	FEDomainMap* map = new FEDomainMap(FE_VEC3D, FMT_ITEM);
	map->Create(&set);
	FEMappedValueVec3* val = fecore_new<FEMappedValueVec3>("map", GetFEModel());
	val->setDataMap(map);
	fiber->SetProperty(val);

	vector<ELEMENT_DATA> data;
	ParseElementData(tag, set, data, 3);
	for (int i = 0; i < (int)data.size(); ++i)
	{
		ELEMENT_DATA& di = data[i];
		if (di.nval > 0)
		{
			FEElement& el = set.Element(i);

			if (di.nval != 3) throw XMLReader::InvalidTag(tag);
			vec3d v(di.val[0], di.val[1], di.val[2]);
			v.unit();
			map->set<vec3d>(i, v);
		}
	}
}


//-----------------------------------------------------------------------------
void FEBioMeshDataSection4::ParseMaterialAxes(XMLTag& tag, FEElementSet& set)
{
	// find the domain with the same name
	string name = set.GetName();
	const char* szname = name.c_str();

	FEMesh* mesh = const_cast<FEMesh*>(set.GetMesh());

	// find the domain
	string domName = set.GetName();
	FEDomainList& DL = set.GetDomainList();
	if (DL.Domains() != 1)
	{
		throw XMLReader::InvalidAttributeValue(tag, "elem_set", domName.c_str());
	}
	FEDomain* dom = DL.GetDomain(0);

	// get the material
	FEMaterial* mat = dom->GetMaterial();
	if (mat == nullptr) throw XMLReader::InvalidAttributeValue(tag, "elem_set", szname);

	Storage_Fmt fmt = FMT_ITEM;
	const char* szfmt = tag.AttributeValue("format", true);
	if (szfmt)
	{
		if (szcmp(szfmt, "mat_points") == 0) fmt = FMT_MATPOINTS;
	}

	// get the mat_axis property
	FEProperty* pQ = mat->FindProperty("mat_axis", true);
	if (pQ == nullptr)
	{
		// if the material does not have the mat_axis property, we'll assign it directly to the material points
		// This only works for ITEM storage
		if (fmt != FMT_ITEM) throw XMLReader::InvalidAttributeValue(tag, "format", szfmt);

		++tag;
		do
		{
			if ((tag == "e") || (tag == "elem"))
			{
				// get the local element number
				const char* szlid = tag.AttributeValue("lid");
				int lid = atoi(szlid) - 1;

				// make sure the number is valid
				if ((lid < 0) || (lid >= set.Elements())) throw XMLReader::InvalidAttributeValue(tag, "lid", szlid);

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
				vec3d e3 = v1 ^ v2;
				vec3d e2 = e3 ^ e1;

				// normalize
				e1.unit();
				e2.unit();
				e3.unit();

				// set the value
				mat3d A(e1, e2, e3);

				// convert to quaternion
				quatd Q(A);

				// assign to all material points
				int ni = el->GaussPoints();
				for (int n=0; n<ni; ++n)
				{
					FEMaterialPoint* mp = el->GetMaterialPoint(n);
					mp->m_Q = Q;
				}
			}
			else throw XMLReader::InvalidTag(tag);
			++tag;
		} while (!tag.isend());
		return;
	}

	if (pQ->GetSuperClassID() != FEMAT3DVALUATOR_ID) throw XMLReader::InvalidAttributeValue(tag, "elem_set", szname);

	// create the map's name: material_name.mat_axis
	stringstream ss;
	ss << "material" << mat->GetID() << ".mat_axis";
	string mapName = ss.str();

	// the domain map we're about to create
	FEDomainMap* map = nullptr;

	// see if the generator is defined
	const char* szgen = tag.AttributeValue("generator", true);
	if (szgen)
	{
		// create a domain map
		map = new FEDomainMap(FE_MAT3D, fmt);
		map->SetName(mapName);
		map->Create(&set);

		// data will be generated
		FEModel* fem = GetFEModel();
		FEElemDataGenerator* gen = 0;
		if (strcmp(szgen, "const") == 0) gen = new FEConstDataGenerator<mat3d, FEElemDataGenerator>(fem);
		else
		{
			gen = fecore_new<FEElemDataGenerator>(szgen, fem);
		}
		if (gen == 0) throw XMLReader::InvalidAttributeValue(tag, "generator", szgen);

		// read the parameters
		ReadParameterList(tag, gen);

		// initialize the generator
		if (gen->Init() == false) throw FEBioImport::DataGeneratorError();

		// generate the data
		if (gen->Generate(*map) == false) throw FEBioImport::DataGeneratorError();
	}
	else
	{
		// This only works for ITEM storage
		if (fmt != FMT_ITEM) throw XMLReader::InvalidAttributeValue(tag, "format", szfmt);

		// create a domain map
		map = new FEDomainMap(FE_MAT3D, FMT_ITEM);
		map->SetName(mapName);
		map->Create(&set);

		++tag;
		do
		{
			if ((tag == "e") || (tag == "elem"))
			{
				// get the local element number
				const char* szlid = tag.AttributeValue("lid");
				int lid = atoi(szlid) - 1;

				// make sure the number is valid
				if ((lid < 0) || (lid >= set.Elements())) throw XMLReader::InvalidAttributeValue(tag, "lid", szlid);

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
				vec3d e3 = v1 ^ v2;
				vec3d e2 = e3 ^ e1;

				// normalize
				e1.unit();
				e2.unit();
				e3.unit();

				// set the value
				mat3d Q(e1, e2, e3);
				map->setValue(lid, Q);
			}
			else throw XMLReader::InvalidTag(tag);
			++tag;
		} while (!tag.isend());
	}
	assert(map);

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
		pQ->SetProperty(val);
		mesh->AddDataMap(map);
	}
}
