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
		if (tag == "ElementData")
		{
			ParseElementDataSection(tag);
		}
		else if (tag == "SurfaceData")
		{
			ParseSurfaceDataSection(tag);
		}
		else if (tag == "EdgeData")
		{
			ParseEdgeDataSection(tag);
		}
		else if (tag == "NodeData")
		{
			ParseNodeDataSection(tag);
		}
		else throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());
}

//-----------------------------------------------------------------------------
void FEBioMeshDataSection3::ParseNodeDataSection(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	const char* szset = tag.AttributeValue("node_set");
	FENodeSet* nodeSet = mesh.FindNodeSet(szset);
	if (nodeSet == 0) throw XMLReader::InvalidAttributeValue(tag, "node_set", szset);

	const char* sztype = tag.AttributeValue("data_type", true);
	if (sztype == 0) sztype = "scalar";

	int dataType = -1;
	if (strcmp(sztype, "scalar") == 0) dataType = FE_DOUBLE;
	else if (strcmp(sztype, "vec2") == 0) dataType = FE_VEC2D;
	else if (strcmp(sztype, "vec3") == 0) dataType = FE_VEC3D;
	if (dataType == -1) throw XMLReader::InvalidAttributeValue(tag, "data_type", sztype);

	const char* szname = tag.AttributeValue("name");

	FENodeDataMap* pdata = new FENodeDataMap(dataType);
	fem.AddDataArray(szname, pdata);

	pdata->Create(nodeSet->size());

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
				if (gen.Generate(*pdata, *nodeSet) == false) throw XMLReader::InvalidValue(tag);
			}
			else throw XMLReader::InvalidTag(tag);
			++tag;
		} while (!tag.isend());
	}
	else
	{
		FEBioImport* feb = GetFEBioImport();
		feb->ParseDataArray(tag, *pdata, "node");
	}
}

//-----------------------------------------------------------------------------
void FEBioMeshDataSection3::ParseEdgeDataSection(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

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
	feb->ParseDataArray(tag, *pdata, "edge");
	*/
}

//-----------------------------------------------------------------------------
void FEBioMeshDataSection3::ParseSurfaceDataSection(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	const char* szsurf = tag.AttributeValue("surface");
	FEFacetSet* psurf = mesh.FindFacetSet(szsurf);
	if (psurf == 0) throw XMLReader::InvalidAttributeValue(tag, "surface", szsurf);

	const char* sztype = tag.AttributeValue("data_type", true);
	if (sztype == 0) sztype = "scalar";

	int dataType = -1;
	if (strcmp(sztype, "scalar") == 0) dataType = FE_DOUBLE;
	else if (strcmp(sztype, "vec2") == 0) dataType = FE_VEC2D;
	else if (strcmp(sztype, "vec3") == 0) dataType = FE_VEC3D;
	if (dataType == -1) throw XMLReader::InvalidAttributeValue(tag, "data_type", sztype);

	const char* szname = tag.AttributeValue("name");
	FESurfaceMap* pdata = new FESurfaceMap(dataType);
	fem.AddDataArray(szname, pdata);

	pdata->SetName(szname);
	pdata->Create(psurf);

	const char* szgen = tag.AttributeValue("generator", true);
	if (szgen)
	{
		// data will be generated
		FEDataGenerator* gen = fecore_new<FEDataGenerator>(FEDATAGENERATOR_ID, szgen, &fem);
		if (gen == 0) throw XMLReader::InvalidAttributeValue(tag, "generator", szgen);

		// read the parameters
		ReadParameterList(tag, gen);

		// give the generator a chance to validate itself
		if (gen->Init() == false) throw FEBioImport::DataGeneratorError();

		// generate the data
		if (gen->Generate(*pdata, *psurf) == false) throw FEBioImport::DataGeneratorError();
	}
	else
	{
		FEBioImport* feb = GetFEBioImport();
		feb->ParseDataArray(tag, *pdata, "face");
	}
}

//-----------------------------------------------------------------------------
void FEBioMeshDataSection3::ParseElementDataSection(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// get the element set
	const char* szset = tag.AttributeValue("elem_set");
	FEElementSet* part = mesh.FindElementSet(szset);
	if (part == 0) throw XMLReader::InvalidAttributeValue(tag, "elem_set", szset);

	// see if the data will be generated or tabulated
	const char* szgen = tag.AttributeValue("generator", true);

	// get the data type
	const char* sztype = tag.AttributeValue("data_type", true);
	if (sztype == 0) sztype = "scalar";

	int dataType = -1;
	if      (strcmp(sztype, "scalar") == 0) dataType = FE_DOUBLE;
	else if (strcmp(sztype, "vec2"  ) == 0) dataType = FE_VEC2D;
	else if (strcmp(sztype, "vec3"  ) == 0) dataType = FE_VEC3D;
	if (dataType == -1) throw XMLReader::InvalidAttributeValue(tag, "data_type", sztype);

	// see if it references a variable
	const char* szvar = tag.AttributeValue("var", true);

	// see if it defines a name
	const char* szname = tag.AttributeValue("name", true);

	// either var or name must be defined
	if ((szvar == 0) && (szname == 0)) throw XMLReader::MissingAttribute(tag, "[var|name]");
	if (szvar) szname = szvar;

	if (szgen == 0)
	{
		// data is tabulated and mapped directly to a variable.
		if      (strcmp(szvar, "shell thickness") == 0) ParseShellThickness(tag, *part);
		else if (strcmp(szvar, "fiber"    ) == 0) ParseMaterialFibers(tag, *part);
		else if (strcmp(szvar, "mat_axis" ) == 0) ParseMaterialAxes(tag, *part);
		else if (strstr(szvar, ".fiber"   ) != 0) ParseMaterialFiberProperty(tag, *part);
		else if (strstr(szvar, ".mat_axis") != 0) ParseMaterialAxesProperty(tag, *part);
		else ParseMaterialData(tag, *part, szvar);
	}
	else
	{
		// data will be generated
		FEDataGenerator* gen = fecore_new<FEDataGenerator>(FEDATAGENERATOR_ID, szgen, &fem);
		if (gen == 0) throw XMLReader::InvalidAttributeValue(tag, "generator", szgen);

		// read the parameters
		ReadParameterList(tag, gen);

		// initialize the generator
		if (gen->Init() == false) throw FEBioImport::DataGeneratorError();

		// create a new domain map
		FEDomainMap* map = new FEDomainMap(FE_DOUBLE);
		map->Create(part);
		map->SetName(szname);
		fem.AddDataArray(szname, map);

		// if a variable is defined, apply the map to the variable
		if (szvar)
		{
			string pname = szvar;

			// assume the first part is the name of the material
			size_t n = pname.find('.');
			if (n == string::npos) throw XMLReader::InvalidAttributeValue(tag, "var", pname.c_str());

			string matName = pname.substr(0, n);
			string paramName = pname.substr(n + 1);

			// get the material
			FEMaterial* mat = GetFEModel()->FindMaterial(matName);
			if (mat == 0) throw XMLReader::InvalidAttributeValue(tag, "var", pname.c_str());

			FEParameterList& PL = mat->GetParameterList();
			FEParam* param = PL.FindFromName(paramName.c_str());
			assert(param);
			if (param == 0) throw XMLReader::InvalidAttributeValue(tag, "var", pname.c_str());;

			if (param->type() != FE_PARAM_DOUBLE_MAPPED) throw XMLReader::InvalidAttributeValue(tag, "var", pname.c_str());

			FEParamDouble& mp = param->value<FEParamDouble>();
			mp.setValuator(new FEMappedValue(map));
		}

		// generate the data
		if (gen->Generate(*map, *part) == false) throw XMLReader::InvalidValue(tag);
	}
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
void FEBioMeshDataSection3::ParseMaterialData(XMLTag& tag, FEElementSet& set, const string& pname)
{
	vector<ELEMENT_DATA> data;
	ParseElementData(tag, set, data, 1);

	// assume the first part is the name of the material
	size_t n = pname.find('.');
	if (n == string::npos) throw XMLReader::InvalidAttributeValue(tag, "var", pname.c_str());

	string matName = pname.substr(0, n - 1);
	string paramName = pname.substr(n + 1);

	// get the material
	FEMaterial* mat = GetFEModel()->FindMaterial(matName);
	if (mat == 0) throw XMLReader::InvalidAttributeValue(tag, "var", pname.c_str());

	FEParameterList& PL = mat->GetParameterList();
	FEParam* param = PL.FindFromName(paramName.c_str());
	assert(param);
	if (param == 0) return;

	if (param->type() != FE_PARAM_DOUBLE_MAPPED) return;

	FEParamDouble& mp = param->value<FEParamDouble>();
	
	FEDomainMap* map = new FEDomainMap(FE_DOUBLE);
	map->Create(&set);
	
	assert(set.Elements() == (int)data.size());
	for (int i = 0; i < set.Elements(); ++i) map->setValue(i, data[i].val[0]);
	GetFEModel()->AddDataArray(pname.c_str(), map);

	mp.setValuator(new FEMappedValue(map));
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
