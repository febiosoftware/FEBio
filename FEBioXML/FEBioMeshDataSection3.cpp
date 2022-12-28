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

//-----------------------------------------------------------------------------
#ifdef WIN32
#define szcmp    _stricmp
#else
#define szcmp    strcmp
#endif

//-----------------------------------------------------------------------------
// helper function for converting a datatype attribute to FEDataType
FEDataType str2datatype(const char* szdataType)
{
	FEDataType dataType = FEDataType::FE_INVALID_TYPE;
	if      (strcmp(szdataType, "scalar") == 0) dataType = FEDataType::FE_DOUBLE;
	else if (strcmp(szdataType, "vec2"  ) == 0) dataType = FEDataType::FE_VEC2D;
	else if (strcmp(szdataType, "vec3"  ) == 0) dataType = FEDataType::FE_VEC3D;
	else if (strcmp(szdataType, "mat3"  ) == 0) dataType = FEDataType::FE_MAT3D;
    else if (strcmp(szdataType, "mat3s" ) == 0) dataType = FEDataType::FE_MAT3DS;
	return dataType;
}

//-----------------------------------------------------------------------------
void FEBioMeshDataSection3::Parse(XMLTag& tag)
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
void FEBioMeshDataSection3::ParseNodalData(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// find the element set
	const char* szset = tag.AttributeValue("node_set");

	// find the element set in the mesh
	FENodeSet* nset = GetBuilder()->FindNodeSet(szset);
	if (nset == nullptr) throw XMLReader::InvalidAttributeValue(tag, "node_set", szset);

	// get the data type
	const char* szdataType = tag.AttributeValue("datatype", true);
	if (szdataType == nullptr) szdataType = "scalar";
	FEDataType dataType = str2datatype(szdataType);
	if (dataType == FEDataType::FE_INVALID_TYPE) throw XMLReader::InvalidAttributeValue(tag, "datatype", szdataType);

	// get the name (required!)
	const char* szname = tag.AttributeValue("name");

	// create the data map
	FENodeDataMap* map = new FENodeDataMap(dataType);
	map->Create(nset);
	map->SetName(szname);

	// add it to the mesh
	mesh.AddDataMap(map);

	// see if there is a generator
	const char* szgen = tag.AttributeValue("generator", true);
	if (szgen)
	{
		FENodeDataGenerator* gen = 0;
		// data will be generated
		if (strcmp(szgen, "const") == 0)
		{
			if      (dataType == FE_DOUBLE) gen = new FEConstDataGenerator<double, FENodeDataGenerator>(&fem);
			else if (dataType == FE_VEC3D ) gen = new FEConstDataGenerator<vec3d, FENodeDataGenerator>(&fem);
			else if (dataType == FE_MAT3D ) gen = new FEConstDataGenerator<mat3d, FENodeDataGenerator>(&fem);
            else if (dataType == FE_MAT3DS) gen = new FEConstDataGenerator<mat3ds, FENodeDataGenerator>(&fem);
		}
		else if (strcmp(szgen, "math") == 0)
		{
			gen = new FEDataMathGenerator(&fem);
		}
		else
		{
			gen = fecore_new<FENodeDataGenerator>(szgen, &fem);
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
		// read the data
		ParseNodeData(tag, *map);
	}
}

void FEBioMeshDataSection3::ParseSurfaceData(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// find the element set
	const char* szset = tag.AttributeValue("surface");

	// find the element set in the mesh
	FEFacetSet* surf = mesh.FindFacetSet(szset);
	if (surf == nullptr) throw XMLReader::InvalidAttributeValue(tag, "surface", szset);

	// get the data type
	const char* szdataType = tag.AttributeValue("datatype", true);
	if (szdataType == nullptr) szdataType = "scalar";
	FEDataType dataType = str2datatype(szdataType);
	if (dataType == FEDataType::FE_INVALID_TYPE) throw XMLReader::InvalidAttributeValue(tag, "datatype", szdataType);

	// get the name (required!)
	const char* szname = tag.AttributeValue("name");

	// create the data map
	FESurfaceMap* map = new FESurfaceMap(dataType);
	map->Create(surf);
	map->SetName(szname);

	// add it to the mesh
	mesh.AddDataMap(map);

	// see if there is a generator
	const char* szgen = tag.AttributeValue("generator", true);
	if (szgen)
	{
		FEFaceDataGenerator* gen = 0;
		// data will be generated
		if (strcmp(szgen, "const") == 0)
		{
			if      (dataType == FE_DOUBLE) gen = new FEConstDataGenerator<double, FEFaceDataGenerator>(&fem);
			else if (dataType == FE_VEC3D ) gen = new FEConstDataGenerator<vec3d , FEFaceDataGenerator>(&fem);
			else if (dataType == FE_MAT3D ) gen = new FEConstDataGenerator<mat3d , FEFaceDataGenerator>(&fem);
            else if (dataType == FE_MAT3DS) gen = new FEConstDataGenerator<mat3ds, FEFaceDataGenerator>(&fem);
		}
		else
		{
			gen = fecore_new<FEFaceDataGenerator>(szgen, &fem);
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
		// read the data
		ParseSurfaceData(tag, *map);
	}
}

void FEBioMeshDataSection3::ParseElementData(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// find the element set
	const char* szset = tag.AttributeValue("elem_set");

	// find the element set in the mesh
	FEElementSet* elset = mesh.FindElementSet(szset);
	if (elset == nullptr) throw XMLReader::InvalidAttributeValue(tag, "elem_set", szset);

	// get the data type
	const char* szdataType = tag.AttributeValue("datatype", true);
	if (szdataType == nullptr) szdataType = "scalar";
	FEDataType dataType = str2datatype(szdataType);
	if (dataType == FEDataType::FE_INVALID_TYPE) throw XMLReader::InvalidAttributeValue(tag, "datatype", szdataType);

	// see if a generator is defined
	const char* szgen = tag.AttributeValue("generator", true);

	// get the name or var (required!)
	const char* szvar = tag.AttributeValue("var", true);
	const char* szname = (szvar == nullptr ? tag.AttributeValue("name") : nullptr);

	bool isVar = (szvar != nullptr);
	string mapName = (isVar ? szvar : szname);

	// process some special variables
	if ((szname == nullptr) && (szgen == nullptr))
	{
		if      (strcmp(szvar, "shell thickness") == 0) ParseShellThickness(tag, *elset);
		else if (strcmp(szvar, "fiber"          ) == 0) ParseMaterialFibers(tag, *elset);
		else if (strcmp(szvar, "mat_axis"       ) == 0) ParseMaterialAxes(tag, *elset);
		else throw XMLReader::InvalidAttributeValue(tag, "var");
		return;
	}

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
	map->SetName(mapName);

	// see if there is a generator
	if (szgen)
	{
		// make sure the parameter is valid
		FEParamDouble* pp = nullptr;
		if (isVar)
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

		FEElemDataGenerator* gen = nullptr;		// data will be generated
		if (strcmp(szgen, "const") == 0)
		{
			if      (dataType == FE_DOUBLE) gen = new FEConstDataGenerator<double, FEElemDataGenerator>(&fem);
			else if (dataType == FE_VEC3D ) gen = new FEConstDataGenerator<vec3d , FEElemDataGenerator>(&fem);
			else if (dataType == FE_MAT3D ) gen = new FEConstDataGenerator<mat3d , FEElemDataGenerator>(&fem);
            else if (dataType == FE_MAT3DS) gen = new FEConstDataGenerator<mat3ds, FEElemDataGenerator>(&fem);
		}
		else
		{
			gen = fecore_new<FEElemDataGenerator>(szgen, &fem);
		}
		if (gen == 0) throw XMLReader::InvalidAttributeValue(tag, "generator", szgen);

		// read the parameters
		ReadParameterList(tag, gen);

		// Add it to the list (will be applied after the rest of the model was read in)
		GetBuilder()->AddMeshDataGenerator(gen, map, pp);
	}
	else
	{
		string name = szname;

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
/*
void FEBioMeshDataSection3::ParseModelParameter(XMLTag& tag, FEParamValue param)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// get the parameter name
	const char* szparam = tag.AttributeValue("param");

	// see if the data will be generated or tabulated
	const char* szgen = tag.AttributeValue("generator", true);

	FEParam* pp = param.param();
	if (pp->type() == FE_PARAM_MATERIALPOINT)
	{
		ParseMaterialPointData(tag, param);
		return;
	}

	// make sure it is a mapped param
	if ((pp->type() != FE_PARAM_DOUBLE_MAPPED) &&
		(pp->type() != FE_PARAM_VEC3D_MAPPED) && 
		(pp->type() != FE_PARAM_MAT3D_MAPPED)) throw XMLReader::InvalidAttributeValue(tag, "param", szparam);

	FEDataType dataType = FE_DOUBLE;
	if (pp->type() == FE_PARAM_VEC3D_MAPPED) dataType = FE_VEC3D;
	if (pp->type() == FE_PARAM_MAT3D_MAPPED) dataType = FE_MAT3D;

	// get the parent
	FECoreBase* pc = dynamic_cast<FECoreBase*>(pp->parent());
	if (pc == 0) throw XMLReader::InvalidAttributeValue(tag, "param", szparam);

	// data generator
	// TODO: Make this a shared pointer so it will be deleted properly
	FEDataGenerator* gen = 0;
	if (szgen)
	{
		// data will be generated
		if (strcmp(szgen, "const") == 0)
		{
			if      (dataType == FE_DOUBLE) gen = new FEConstDataGenerator<double>(&fem);
			else if (dataType == FE_VEC3D ) gen = new FEConstDataGenerator<vec3d>(&fem);
			else if (dataType == FE_MAT3D ) gen = new FEConstDataGenerator<mat3d>(&fem);
		}
		else
		{
			gen = fecore_new<FEDataGenerator>(szgen, &fem);
		}
		if (gen == 0) throw XMLReader::InvalidAttributeValue(tag, "generator", szgen);

		// read the parameters
		ReadParameterList(tag, gen);

		// initialize the generator
		if (gen->Init() == false) throw FEBioImport::DataGeneratorError();
	}

	if (dynamic_cast<FEMaterial*>(pc))
	{
		FEMaterial* mat = dynamic_cast<FEMaterial*>(pc->GetAncestor());
		if (mat == 0) throw XMLReader::InvalidAttributeValue(tag, "param", szparam);

		FEDomainList& DL = mat->GetDomainList();
		FEElementSet* set = new FEElementSet(&fem);
		set->Create(DL);
		mesh.AddElementSet(set);

		// create a new domain map
		FEDomainMap* map = nullptr;
		switch (dataType)
		{
		case FE_DOUBLE: map = new FEDomainMap(dataType); break;
		case FE_VEC3D : map = new FEDomainMap(dataType); break;
		case FE_MAT3D : map = new FEDomainMap(dataType, FMT_ITEM); break;
		}

		map->Create(set);
		map->SetName(szparam);
		mesh.AddDataArray(szparam, map);

		if (gen)
		{
			if (gen->Generate(*map, *set) == false) throw FEBioImport::DataGeneratorError();
		}
		else
		{
			ParseElementData(tag, *map);
		}

		if (dataType == FE_DOUBLE)
		{
			FEParamDouble& p = pp->value<FEParamDouble>();
			if (p.isConst() == false) throw FEBioImport::DataGeneratorError();
			FEMappedValue* val = new FEMappedValue(&fem);
			val->setDataMap(map, p.constValue());
			p.setValuator(val);
		}
		else if (dataType == FE_VEC3D)
		{
			FEParamVec3& p = pp->value<FEParamVec3>();
			FEMappedValueVec3* val = new FEMappedValueVec3(&fem);
			val->setDataMap(map);
			p.setValuator(val);
		}
		else if (dataType == FE_MAT3D)
		{
			FEParamMat3d& p = pp->value<FEParamMat3d>();
			FEMappedValueMat3d* val = fecore_alloc(FEMappedValueMat3d, &fem);
			val->setDataMap(map);
			p.setValuator(val);
		}

		int a = 0;
	}
	else if (dynamic_cast<FEBodyLoad*>(pc))
	{
		FEBodyLoad* pbl = dynamic_cast<FEBodyLoad*>(pc);

		FEDomainList& DL = pbl->GetDomainList();
		FEElementSet* set = new FEElementSet(&fem);
		set->Create(DL);
		mesh.AddElementSet(set);

		// create a new domain map
		FEDomainMap* map = new FEDomainMap(dataType);
		map->Create(set);
		map->SetName(szparam);
		mesh.AddDataArray(szparam, map);

		if (gen)
		{
			if (gen->Generate(*map, *set) == false) throw FEBioImport::DataGeneratorError();
		}
		else
		{
			ParseElementData(tag, *map);
		}

		if (dataType == FE_DOUBLE)
		{
			FEParamDouble& p = pp->value<FEParamDouble>();
			if (p.isConst() == false) throw FEBioImport::DataGeneratorError();
			FEMappedValue* val = new FEMappedValue(&fem);
			val->setDataMap(map, p.constValue());
			p.setValuator(val);
		}
		else if (dataType == FE_VEC3D)
		{
			FEParamVec3& p = pp->value<FEParamVec3>();
			if (p.isConst() == false) throw FEBioImport::DataGeneratorError();
			FEMappedValueVec3* val = new FEMappedValueVec3(&fem);
			val->setDataMap(map, p.constValue());
			p.setValuator(val);
		}
	}
	else if (dynamic_cast<FESurfaceLoad*>(pc))
	{
		FESurfaceLoad* psl = dynamic_cast<FESurfaceLoad*>(pc);

		FESurface* set = &psl->GetSurface();

		// create a new domain map
		FESurfaceMap* map = new FESurfaceMap(dataType);
		map->Create(set);
		map->SetName(szparam);
		mesh.AddDataArray(szparam, map);

		if (gen)
		{
			if (gen->Generate(*map, *set->GetFacetSet()) == false) throw FEBioImport::DataGeneratorError();
		}
		else
		{
			ParseSurfaceData(tag, *map);
		}

		if (dataType == FE_DOUBLE)
		{
			FEParamDouble& p = pp->value<FEParamDouble>();
			if (p.isConst() == false) throw FEBioImport::DataGeneratorError();
			FEMappedValue* val = new FEMappedValue(&fem);
			val->setDataMap(map, p.constValue());
			p.setValuator(val);
		}
		else if (dataType == FE_VEC3D)
		{
			FEParamVec3& p = pp->value<FEParamVec3>();
			if (p.isConst() == false) throw FEBioImport::DataGeneratorError();
			FEMappedValueVec3* val = new FEMappedValueVec3(&fem);
			val->setDataMap(map, p.constValue());
			p.setValuator(val);
		}
	}
	else if (dynamic_cast<FEPrescribedDOF*>(pc))
	{
		FEPrescribedDOF* bc = dynamic_cast<FEPrescribedDOF*>(pc);
		// create node set
		const FENodeSet& bc_set = *bc->GetNodeSet();
		int nsize = bc_set.Size();
		FENodeSet* set = new FENodeSet(&fem);
		for (int i = 0; i < nsize; ++i) set->Add(bc_set[i]);

		FENodeDataMap* map = new FENodeDataMap(FE_DOUBLE);
		mesh.AddDataArray(szparam, map);

		if (gen)
		{
			if (gen->Generate(*map, *set) == false) throw FEBioImport::DataGeneratorError();
		}
		else
		{
			map->Create(set->Size());
			ParseNodeData(tag, *map);
		}

		if (dataType == FE_DOUBLE)
		{
			FEParamDouble& p = pp->value<FEParamDouble>();
			if (p.isConst() == false) throw FEBioImport::DataGeneratorError();
			FENodeMappedValue* val = new FENodeMappedValue(&fem);
			val->setDataMap(map, p.constValue());
			p.setValuator(val);
		}
	}
	else throw XMLReader::InvalidAttributeValue(tag, "param", szparam);
}

//-----------------------------------------------------------------------------
// Helper function for setting material point member data
template <class T> void setMaterialPointData(FEElement& el, FEMaterialPointProperty& d, const T& v)
{
	int nint = el.GaussPoints();
	for (int j = 0; j < nint; ++j)
	{
		FEMaterialPoint& mp = *el.GetMaterialPoint(j);
		d.set(mp, v);
	}
}

void FEBioMeshDataSection3::ParseMaterialPointData(XMLTag& tag, FEParamValue param)
{
	if (param.type() != FE_PARAM_MATERIALPOINT) throw XMLReader::InvalidAttributeValue(tag, "param");

	FEParam* pp = param.param();
	FEMaterialPointProperty& matProp = pp->value<FEMaterialPointProperty>();
	FECoreBase* pc = dynamic_cast<FECoreBase*>(pp->parent());
	if (pc == 0) throw XMLReader::InvalidAttributeValue(tag, "param");

	FEDataType dataType = matProp.dataType();

	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// get the parameter name
	const char* szparam = tag.AttributeValue("param");

	// see if the data will be generated or tabulated
	const char* szgen = tag.AttributeValue("generator", true);

	FEMaterial* mat = dynamic_cast<FEMaterial*>(pc->GetAncestor());
	if (mat)
	{
		FEDomainList& DL = mat->GetDomainList();
		FEElementSet* set = new FEElementSet(&fem);
		set->Create(DL);

		if (szgen)
		{
			FEDataGenerator* gen = 0;
			if (strcmp(szgen, "const") == 0)
			{
				if      (dataType == FE_DOUBLE) gen = new FEConstDataGenerator<double>(&fem);
				else if (dataType == FE_VEC2D ) gen = new FEConstDataGenerator<vec2d >(&fem);
				else if (dataType == FE_VEC3D ) gen = new FEConstDataGenerator<vec3d >(&fem);
				else if (dataType == FE_MAT3D ) gen = new FEConstDataGenerator<mat3d >(&fem);
			}
			else
			{
				gen = fecore_new<FEDataGenerator>(szgen, &fem);
			}
			if (gen == 0) throw XMLReader::InvalidAttributeValue(tag, "generator", szgen);

			// read the parameters
			ReadParameterList(tag, gen);

			// initialize the generator
			if (gen->Init() == false) throw FEBioImport::DataGeneratorError();

			for (int i = 0; i < set->Elements(); ++i)
			{
				FEElement& el = set->Element(i);

				int nint = el.GaussPoints();
				for (int n = 0; n < nint; ++n)
				{
					FEMaterialPoint& mp = *el.GetMaterialPoint(n);
					vec3d rn = mp.m_r0;
					switch (dataType)
					{
					case FE_DOUBLE: { double d; gen->value(rn, d); matProp.set(mp, d); } break;
					case FE_VEC2D : { vec2d  d; gen->value(rn, d); matProp.set(mp, d); } break;
					case FE_VEC3D : { vec3d  d; gen->value(rn, d); matProp.set(mp, d); } break;
					case FE_MAT3D : { mat3d  d; gen->value(rn, d); matProp.set(mp, d); } break;
					}
				}
			}
		}
		else
		{
			vector<ELEMENT_DATA> data;
			ParseElementData(tag, *set, data, fecore_data_size(matProp.dataType()));

			for (int i = 0; i<set->Elements(); ++i)
			{
				FEElement& el = set->Element(i);

				ELEMENT_DATA& di = data[i];

				// make sure the correct number of values were read in
				if (di.nval != fecore_data_size(matProp.dataType()))
				{
					throw FEBioImport::MeshDataError();
				}

				double* v = di.val;
				switch (matProp.dataType())
				{
				case FE_DOUBLE: setMaterialPointData(el, matProp, v[0]); break;
				case FE_VEC2D : setMaterialPointData(el, matProp, vec2d(v[0], v[1])); break;
				case FE_VEC3D : setMaterialPointData(el, matProp, vec3d(v[0], v[1], v[2])); break;
				case FE_MAT3D : setMaterialPointData(el, matProp, mat3d(v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8])); break;
				default:
					assert(false);
				}
			}
		}

		delete set;

		return;
	}
	
	throw XMLReader::InvalidAttributeValue(tag, "param", szparam);
}
*/
//-----------------------------------------------------------------------------
void FEBioMeshDataSection3::ParseShellThickness(XMLTag& tag, FEElementSet& set)
{
	if (tag.isleaf())
	{
		FEMesh& mesh = GetFEModel()->GetMesh();
		double h[FEElement::MAX_NODES];
		int nval = tag.value(h, FEElement::MAX_NODES);

		for (int i=0; i<set.Elements(); ++i)
		{
			FEShellElement* pel = dynamic_cast<FEShellElement*>(&set.Element(i));
			if (pel == 0) throw XMLReader::InvalidValue(tag);

			if (pel->Nodes() != nval) throw XMLReader::InvalidValue(tag);
			for (int j=0; j<nval; ++j) pel->m_h0[j] = h[j];
		}
	}
	else
	{
		vector<ELEMENT_DATA> data;
		ParseElementData(tag, set, data, FEElement::MAX_NODES);
		for (int i=0; i<(int)data.size(); ++i)
		{
			ELEMENT_DATA& di = data[i];
			if (di.nval > 0)
			{
				FEElement& el = set.Element(i);

				if (el.Class() != FE_ELEM_SHELL) throw XMLReader::InvalidTag(tag);
				FEShellElement& shell = static_cast<FEShellElement&> (el);
		
				int ne = shell.Nodes();
				if (ne != di.nval) throw XMLReader::InvalidTag(tag);
				for (int j=0; j<ne; ++j) shell.m_h0[j] = di.val[j];
			}
		}
	}
}

//-----------------------------------------------------------------------------
void FEBioMeshDataSection3::ParseMaterialFibers(XMLTag& tag, FEElementSet& set)
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
void FEBioMeshDataSection3::ParseMaterialAxes(XMLTag& tag, FEElementSet& set)
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
				for (int n = 0; n < ni; ++n)
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

//-----------------------------------------------------------------------------
void FEBioMeshDataSection3::ParseNodeData(XMLTag& tag, FENodeDataMap& map)
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

//-----------------------------------------------------------------------------
void FEBioMeshDataSection3::ParseSurfaceData(XMLTag& tag, FESurfaceMap& map)
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

		int nread = tag.value(data, m*dataSize);
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
		else if (nread == m*dataSize)
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
	}
	while (!tag.isend());
}

//-----------------------------------------------------------------------------
void FEBioMeshDataSection3::ParseElementData(XMLTag& tag, FEDomainMap& map)
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

		int nread = tag.value(data, m*dataSize);
		if (nread == dataSize)
		{
			double* v = data;
			switch (dataType)
			{
			case FE_DOUBLE:	map.setValue(n, v[0]); break;
			case FE_VEC2D :	map.setValue(n, vec2d(v[0], v[1])); break;
			case FE_VEC3D :	map.setValue(n, vec3d(v[0], v[1], v[2])); break;
			case FE_MAT3D : map.setValue(n, mat3d(v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8])); break;
            case FE_MAT3DS: map.setValue(n, mat3ds(v[0], v[1], v[2], v[3], v[4], v[5])); break;
			default:
				assert(false);
			}
		}
		else if (nread == m*dataSize)
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
	}
	while (!tag.isend());

	if (ncount != nelems) throw FEBioImport::MeshDataError();
}

//-----------------------------------------------------------------------------
void FEBioMeshDataSection3::ParseElementData(XMLTag& tag, FEElementSet& set, vector<ELEMENT_DATA>& values, int nvalues)
{
	// get the total nr of elements
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	int nelems = set.Elements();

	// resize the array
	values.resize(nelems);
	for (int i=0; i<nelems; ++i) values[i].nval = 0;

	++tag;
	do
	{
		// get the local element number
		const char* szlid = tag.AttributeValue("lid");
		int n = atoi(szlid)-1;

		// make sure the number is valid
		if ((n<0) || (n>=nelems)) throw XMLReader::InvalidAttributeValue(tag, "lid", szlid);

		ELEMENT_DATA& data = values[n];
		data.nval = tag.value(data.val, nvalues);
		++tag;
	}
	while (!tag.isend());
}
