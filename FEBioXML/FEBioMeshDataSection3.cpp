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

//-----------------------------------------------------------------------------
void FEBioMeshDataSection3::Parse(XMLTag& tag)
{
	// Make sure there is something in this tag
	if (tag.isleaf()) return;

	// loop over all mesh data section
	++tag;
	do
	{
		if (tag == "mesh_data")
		{
			ParseMeshDataSection(tag);
		}
		else throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());
}

//-----------------------------------------------------------------------------
void FEBioMeshDataSection3::ParseMeshDataSection(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();

	// get the parameter name
	const char* szparam = tag.AttributeValue("param");

	/*	if (szgen == 0)
		{
			// data is tabulated and mapped directly to a variable.
			if      (strcmp(szparam, "shell thickness") == 0) ParseShellThickness(tag, *part);
		}
	*/
	// get the model parameter
	ParamString ps(szparam);
	FEParamValue param = fem.GetParameterValue(ps);
	if (param.isValid())
	{
		ParseModelParameter(tag, param);
	}
	else throw XMLReader::InvalidAttributeValue(tag, "param", szparam);
}

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
			gen = fecore_new<FEDataGenerator>(FEDATAGENERATOR_ID, szgen, &fem);
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
		FEElementSet* set = new FEElementSet(&mesh);
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
		fem.AddDataArray(szparam, map);

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
			FEMappedValueMat3d* val = new FEMappedValueMat3d(&fem);
			val->setDataMap(map);
			p.setValuator(val);
		}

		int a = 0;
	}
	else if (dynamic_cast<FEBodyLoad*>(pc))
	{
		FEBodyLoad* pbl = dynamic_cast<FEBodyLoad*>(pc);

		FEDomainList& DL = pbl->GetDomainList();
		FEElementSet* set = new FEElementSet(&mesh);
		set->Create(DL);
		mesh.AddElementSet(set);

		// create a new domain map
		FEDomainMap* map = new FEDomainMap(dataType);
		map->Create(set);
		map->SetName(szparam);
		fem.AddDataArray(szparam, map);

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
		fem.AddDataArray(szparam, map);

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
		int nsize = (int)bc->Items();
		FENodeSet* set = new FENodeSet(&mesh);
		for (int i = 0; i < nsize; ++i) set->add(bc->NodeID(i));

		FENodeDataMap* map = new FENodeDataMap(FE_DOUBLE);
		fem.AddDataArray(szparam, map);

		if (gen)
		{
			if (gen->Generate(*map, *set) == false) throw FEBioImport::DataGeneratorError();
		}
		else
		{
			map->Create(set->size());
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
		FEElementSet* set = new FEElementSet(&mesh);
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
				gen = fecore_new<FEDataGenerator>(FEDATAGENERATOR_ID, szgen, &fem);
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
			FEShellElement* pel = dynamic_cast<FEShellElement*>(m_pelem[set[i]-1]);
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
				FEElement& el = *m_pelem[set[i]-1];

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
	const FESurface* set = map.GetSurface();
	if (set == nullptr) throw XMLReader::InvalidTag(tag);

	// get the total nr of elements
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	int nelems = set->Elements();

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
