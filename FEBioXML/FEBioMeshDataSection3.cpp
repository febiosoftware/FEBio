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
#include <FECore/FEVectorGenerator.h>

//-----------------------------------------------------------------------------
// Defined in FEBioGeometrySection.cpp
void set_element_fiber(FEElement& el, const vec3d& v, int ncomp);
void set_element_mat_axis(FEElement& el, const vec3d& v1, const vec3d& v2, int ncomp);

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
			else if (strcmp(szparam, "fiber"    ) == 0) ParseMaterialFibers(tag, *part);
			else if (strcmp(szparam, "mat_axis" ) == 0) ParseMaterialAxes(tag, *part);
			else if (strstr(szparam, ".fiber"   ) != 0) ParseMaterialFiberProperty(tag, *part);
			else if (strstr(szparam, ".mat_axis") != 0) ParseMaterialAxesProperty(tag, *part);
			else ParseMaterialData(tag, *part, szvar);

		}
	*/
	// get the model parameter
	ParamString ps(szparam);
	FEParamValue param = fem.GetParameterValue(ps);
	if (param.isValid())
	{
		ParseModelParameter(tag, param);
	}
	else
	{
		// the param does not reference a model parameter, but it can be something else
		ParseMeshDataField(tag);
	}
}

void FEBioMeshDataSection3::ParseModelParameter(XMLTag& tag, FEParamValue param)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// get the parameter name
	const char* szparam = tag.AttributeValue("param");

	// see if the data will be generated or tabulated
	const char* szgen = tag.AttributeValue("generator", true);

	// make sure it is a mapped param
	if ((param.type() != FE_PARAM_DOUBLE_MAPPED) &&
		(param.type() != FE_PARAM_VEC3D_MAPPED)) throw XMLReader::InvalidAttributeValue(tag, "param", szparam);

	int dataType = (param.type() == FE_PARAM_DOUBLE_MAPPED ? FE_DOUBLE : FE_VEC3D);

	// get the parent
	FECoreBase* pc = dynamic_cast<FECoreBase*>(param.param()->parent());
	if (pc == 0) throw XMLReader::InvalidAttributeValue(tag, "param", szparam);

	// data generator
	FEDataGenerator* gen = 0;
	if (szgen)
	{
		// data will be generated
		gen = fecore_new<FEDataGenerator>(FEDATAGENERATOR_ID, szgen, &fem);
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
			FEParamDouble& p = param.value<FEParamDouble>();
			if (p.isConst() == false) throw FEBioImport::DataGeneratorError();
			p.setValuator(new FEMappedValue(map, p.constValue()));
		}
		else if (dataType == FE_VEC3D)
		{
			FEParamVec3& p = param.value<FEParamVec3>();
			p.setValuator(new FEMappedValueVec3(map));
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
			FEParamDouble& p = param.value<FEParamDouble>();
			if (p.isConst() == false) throw FEBioImport::DataGeneratorError();
			p.setValuator(new FEMappedValue(map, p.constValue()));
		}
		else if (dataType == FE_VEC3D)
		{
			FEParamVec3& p = param.value<FEParamVec3>();
			if (p.isConst() == false) throw FEBioImport::DataGeneratorError();
			p.setValuator(new FEMappedValueVec3(map, p.constValue()));
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
			FEParamDouble& p = param.value<FEParamDouble>();
			if (p.isConst() == false) throw FEBioImport::DataGeneratorError();
			p.setValuator(new FEMappedValue(map, p.constValue()));
		}
		else if (dataType == FE_VEC3D)
		{
			FEParamVec3& p = param.value<FEParamVec3>();
			if (p.isConst() == false) throw FEBioImport::DataGeneratorError();
			p.setValuator(new FEMappedValueVec3(map, p.constValue()));
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
			FEParamDouble& p = param.value<FEParamDouble>();
			if (p.isConst() == false) throw FEBioImport::DataGeneratorError();
			p.setValuator(new FENodeMappedValue(map, p.constValue()));
		}
	}
	else throw XMLReader::InvalidAttributeValue(tag, "param", szparam);
}

//-----------------------------------------------------------------------------
void FEBioMeshDataSection3::ParseMeshDataField(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// get the parameter name
	const char* szparam = tag.AttributeValue("param");

	// see if the data will be generated or tabulated
	const char* szgen = tag.AttributeValue("generator", true);

	// See if this references a proprety
	ParamString ps(szparam);
	FECoreBase* pp = fem.FindComponent(ps);
	if (pp == nullptr) throw XMLReader::InvalidAttributeValue(tag, "param", szparam);

	if (dynamic_cast<FEUserVectorGenerator*>(pp))
	{
		FEMaterial* mat = dynamic_cast<FEMaterial*>(pp->GetAncestor());
		if (mat)
		{
			FEDomainList& dl = mat->GetDomainList();

			FEElementSet* set = new FEElementSet(&mesh);
			set->Create(dl);

			// create new domain map
			FEDomainMap* map = new FEDomainMap(FE_VEC3D);
			map->Create(set);
			fem.AddDataArray(szparam, map);

			ParseElementData(tag, *map);

			FEUserVectorGenerator* vec = dynamic_cast<FEUserVectorGenerator*>(pp);
			vec->SetData(map);

			return;
		}
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
// Defined in FEBioGeometrySection.cpp
void set_element_fiber(FEElement& el, const vec3d& v, int ncomp);
void set_element_mat_axis(FEElement& el, const vec3d& v1, const vec3d& v2, int ncomp);

//-----------------------------------------------------------------------------
void FEBioMeshDataSection3::ParseMaterialFibers(XMLTag& tag, FEElementSet& set)
{
	vector<ELEMENT_DATA> data;
	ParseElementData(tag, set, data, 3);
	for (int i=0; i<(int)data.size(); ++i)
	{
		ELEMENT_DATA& di = data[i];
		if (di.nval > 0)
		{
			FEElement& el = *m_pelem[set[i]-1];

			if (di.nval != 3) throw XMLReader::InvalidTag(tag);
			vec3d v(di.val[0], di.val[1], di.val[2]);

			set_element_fiber(el, v, -1);
		}
	}
}

//-----------------------------------------------------------------------------
void FEBioMeshDataSection3::ParseMaterialFiberProperty(XMLTag& tag, FEElementSet& set)
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
void FEBioMeshDataSection3::ParseMaterialAxesProperty(XMLTag& tag, FEElementSet& set)
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
void FEBioMeshDataSection3::ParseMaterialAxes(XMLTag& tag, FEElementSet& set)
{
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
			double a[3] = {0};
			double d[3] = {0};
			++tag;
			do
			{
				if (tag == "a") tag.value(a, 3);
				else if (tag == "d") tag.value(d, 3);
				else throw XMLReader::InvalidTag(tag);
				++tag;
			}
			while (!tag.isend());

			vec3d v1(a[0], a[1], a[2]);
			vec3d v2(d[0], d[1], d[2]);
			set_element_mat_axis(*el, v1, v2, -1);
		}
		else throw XMLReader::InvalidTag(tag);
		++tag;	
	}
	while (!tag.isend());
}

//-----------------------------------------------------------------------------
void FEBioMeshDataSection3::ParseNodeData(XMLTag& tag, FENodeDataMap& map)
{
	// get the total nr of nodes
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	int nodes = map.DataCount();

	int dataSize = map.DataSize();
	double data[3]; // make sure this array is large enough to store any data map type (current 3 for FE_VEC3D)

	++tag;
	do
	{
		if (tag == "node")
		{
			// get the local element number
			const char* szlid = tag.AttributeValue("lid");
			int n = atoi(szlid) - 1;

			// make sure the number is valid
			if ((n < 0) || (n >= nodes)) throw XMLReader::InvalidAttributeValue(tag, "lid", szlid);

			int nread = tag.value(data, dataSize);
			if (nread == dataSize)
			{
				switch (dataSize)
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
	} while (!tag.isend());
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

	int dataSize = map.DataSize();
	int m = map.MaxNodes();
	double data[3 * FEElement::MAX_NODES]; // make sure this array is large enough to store any data map type (current 3 for FE_VEC3D)

	++tag;
	do
	{
		if (tag == "face")
		{
			// get the local element number
			const char* szlid = tag.AttributeValue("lid");
			int n = atoi(szlid) - 1;

			// make sure the number is valid
			if ((n < 0) || (n >= nelems)) throw XMLReader::InvalidAttributeValue(tag, "lid", szlid);

			int nread = tag.value(data, m*dataSize);
			if (nread == dataSize)
			{
				switch (dataSize)
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
					switch (dataSize)
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
	} while (!tag.isend());
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

	int dataSize = map.DataSize();
	int m = map.MaxNodes();
	double data[3 * FEElement::MAX_NODES]; // make sure this array is large enough to store any data map type (current 3 for FE_VEC3D)

	++tag;
	do
	{
		if (tag == "elem")
		{
			// get the local element number
			const char* szlid = tag.AttributeValue("lid");
			int n = atoi(szlid) - 1;

			// make sure the number is valid
			if ((n < 0) || (n >= nelems)) throw XMLReader::InvalidAttributeValue(tag, "lid", szlid);

			int nread = tag.value(data, m*dataSize);
			if (nread == dataSize)
			{
				switch (dataSize)
				{
				case FE_DOUBLE:	map.setValue(n, data[0]); break;
				case FE_VEC2D :	map.setValue(n, vec2d(data[0], data[1])); break;
				case FE_VEC3D :	map.setValue(n, vec3d(data[0], data[1], data[2])); break;
				default:
					assert(false);
				}
			}
			else if (nread == m*dataSize)
			{
				double* pd = data;
				for (int i = 0; i < m; ++i, pd += dataSize)
				{
					switch (dataSize)
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
	}
	while (!tag.isend());
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
		if (tag == "elem")
		{
			// get the local element number
			const char* szlid = tag.AttributeValue("lid");
			int n = atoi(szlid)-1;

			// make sure the number is valid
			if ((n<0) || (n>=nelems)) throw XMLReader::InvalidAttributeValue(tag, "lid", szlid);

			ELEMENT_DATA& data = values[n];
			data.nval = tag.value(data.val, nvalues);
		}
		else throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());
}
