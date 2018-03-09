#include "stdafx.h"
#include "FEBioMeshDataSection.h"
#include "FECore/FEModel.h"
#include "FECore/DOFS.h"
#include <FECore/FEDataGenerator.h>
#include <FECore/FECoreKernel.h>
#include <FECore/FEDataMathGenerator.h>

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
	for (int i=0; i<mesh.Domains(); ++i)
	{
		FEDomain& d = mesh.Domain(i);
		for (int j=0; j<d.Elements(); ++j)
		{
			FEElement& el = d.ElementRef(j);
			if (el.GetID() > max_id) max_id = el.GetID();
		}
	}

	// create the pelem array
	m_pelem.assign(max_id, static_cast<FEElement*>(0));
	for (int nd=0; nd<mesh.Domains(); ++nd)
	{
		FEDomain& d = mesh.Domain(nd);
		for (int i=0; i<d.Elements(); ++i)
		{
			FEElement& el = d.ElementRef(i);
			assert(m_pelem[el.GetID()-1] == 0);
			m_pelem[el.GetID()-1] = &el;
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
				const char* szvar = tag.AttributeValue("var");
				if      (strcmp(szvar, "shell thickness") == 0) ParseShellThickness(tag, *part);
				else if (strcmp(szvar, "fiber"          ) == 0) ParseMaterialFibers(tag, *part);
				else if (strcmp(szvar, "mat_axis"       ) == 0) ParseMaterialAxes  (tag, *part);
				else if (strstr(szvar, ".fiber") != 0) ParseMaterialFiberProperty(tag, *part);
				else ParseMaterialData(tag, *part, szvar);
			}
			else
			{
				// data will be generated
				FEDataGenerator* gen = fecore_new<FEDataGenerator>(FEDATAGENERATOR_ID, szgen, &fem);
				if (gen == 0) throw XMLReader::InvalidAttributeValue(tag, "generator", szgen);

				// get the variable
				const char* szvar = tag.AttributeValue("var");

				// read the parameters
				ReadParameterList(tag, gen);

				// get the domain
				FEMesh& mesh = fem.GetMesh();
				FEDomain* dom = mesh.FindDomain(szset);
				if (dom == 0) throw XMLReader::InvalidAttributeValue(tag, "el_set", szset);

				// give the generator a chance to validate itself
				if (gen->Init() == false) throw FEBioImport::DataGeneratorError();

				// generate the data
				if (gen->Apply(dom, szvar) == false) throw FEBioImport::DataGeneratorError();
			}
		}
		else if (tag == "SurfaceData")
		{
			const char* szsurf = tag.AttributeValue("surface");
			FEFacetSet* psurf = mesh.FindFacetSet(szsurf);
			if (psurf == 0) throw XMLReader::InvalidAttributeValue(tag, "surface", szsurf);

			const char* sztype = tag.AttributeValue("data_type", true);
			if (sztype == 0) sztype = "scalar";

			int dataType = -1;
			if      (strcmp(sztype, "scalar") == 0) dataType = FE_DOUBLE;
			else if (strcmp(sztype, "vec2"  ) == 0) dataType = FE_VEC2D;
			else if (strcmp(sztype, "vec3"  ) == 0) dataType = FE_VEC3D;
			if (dataType == -1) throw XMLReader::InvalidAttributeValue(tag, "data_type", sztype);

			const char* szname = tag.AttributeValue("name");
			FESurfaceMap* pdata = new FESurfaceMap(dataType);
			fem.AddDataArray(szname, pdata);

			pdata->SetName(szname);
			pdata->Create(psurf);

			const char* szgen = tag.AttributeValue("generator", true);
			if (szgen)
			{
				if (dataType != FE_DOUBLE) throw XMLReader::InvalidAttributeValue(tag, "generator", szgen);

				++tag;
				do
				{
					if (tag == "math")
					{
						FEDataMathGenerator gen;
						gen.setExpression(tag.szvalue());
						if (gen.Generate(*pdata, *psurf) == false) throw XMLReader::InvalidValue(tag);
					}
					else throw XMLReader::InvalidTag(tag);
					++tag;
				} 
				while (!tag.isend());
			}
			else
				feb->ParseDataArray(tag, *pdata, "face");
		}
		else if (tag == "EdgeData")
		{
			const char* szedge = tag.AttributeValue("edge");
			FESegmentSet* pset = mesh.FindSegmentSet(szedge);
			if (pset == 0) throw XMLReader::InvalidAttributeValue(tag, "edge", szedge);

			const char* sztype = tag.AttributeValue("data_type", true);
			if (sztype == 0) sztype = "scalar";

			int dataType = -1;
			if      (strcmp(sztype, "scalar") == 0) dataType = FE_DOUBLE;
			else if (strcmp(sztype, "vec2"  ) == 0) dataType = FE_VEC2D;
			else if (strcmp(sztype, "vec3"  ) == 0) dataType = FE_VEC3D;
			if (dataType == -1) throw XMLReader::InvalidAttributeValue(tag, "data_type", sztype);
/*
			const char* szname = tag.AttributeValue("name");
			FEDataArray* pdata = new FEDataArray(dataType);
			fem.AddDataArray(szname, pdata);

			pdata->resize(pset->Segments());
			feb->ParseDataArray(tag, *pdata, "edge");
*/
		}
		else if (tag == "NodeData")
		{
			const char* szset = tag.AttributeValue("node_set");
			FENodeSet* nodeSet = mesh.FindNodeSet(szset);
			if (nodeSet == 0) throw XMLReader::InvalidAttributeValue(tag, "node_set", szset);

			const char* sztype = tag.AttributeValue("data_type", true);
			if (sztype == 0) sztype = "scalar";

			int dataType = -1;
			if      (strcmp(sztype, "scalar") == 0) dataType = FE_DOUBLE;
			else if (strcmp(sztype, "vec2"  ) == 0) dataType = FE_VEC2D;
			else if (strcmp(sztype, "vec3"  ) == 0) dataType = FE_VEC3D;
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
						FEDataMathGenerator gen;
						gen.setExpression(tag.szvalue());
						if (gen.Generate(*pdata, *nodeSet) == false) throw XMLReader::InvalidValue(tag);
					}
					else throw XMLReader::InvalidTag(tag);
					++tag;
				}
				while (!tag.isend());
			}
			else
			{
				feb->ParseDataArray(tag, *pdata, "node");
			}
		}
		else throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());
}

//-----------------------------------------------------------------------------
void FEBioMeshDataSection::ParseShellThickness(XMLTag& tag, FEElementSet& set)
{
	if (tag.isleaf())
	{
		FEMesh& mesh = GetFEModel()->GetMesh();
		double h[FEElement::MAX_NODES];
		int nval = tag.value(h, FEElement::MAX_NODES);

		for (int i=0; i<set.size(); ++i)
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
		for (size_t i=0; i<data.size(); ++i)
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
void set_element_mat_axis(FEElement& el, const vec3d& v1, const vec3d& v2);

//-----------------------------------------------------------------------------
void FEBioMeshDataSection::ParseMaterialFibers(XMLTag& tag, FEElementSet& set)
{
	vector<ELEMENT_DATA> data;
	ParseElementData(tag, set, data, 3);
	for (size_t i=0; i<data.size(); ++i)
	{
		ELEMENT_DATA& di = data[i];
		if (di.nval > 0)
		{
			FEElement& el = *m_pelem[set[i]-1];

			if (di.nval != 3) throw XMLReader::InvalidTag(tag);
			vec3d v(di.val[0], di.val[1], di.val[2]);

			set_element_fiber(el, v, 0);
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
	int n = atoi(ch);

	vector<ELEMENT_DATA> data;
	ParseElementData(tag, set, data, 3);
	for (size_t i = 0; i<data.size(); ++i)
	{
		ELEMENT_DATA& di = data[i];
		if (di.nval > 0)
		{
			FEElement& el = *m_pelem[set[i] - 1];

			if (di.nval != 3) throw XMLReader::InvalidTag(tag);
			vec3d v(di.val[0], di.val[1], di.val[2]);

			set_element_fiber(el, v, n);
		}
	}
}

//-----------------------------------------------------------------------------
void FEBioMeshDataSection::ParseMaterialAxes(XMLTag& tag, FEElementSet& set)
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
			if ((lid<0) || (lid >= set.size())) throw XMLReader::InvalidAttributeValue(tag, "lid", szlid);

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
			set_element_mat_axis(*el, v1, v2);
		}
		else throw XMLReader::InvalidTag(tag);
		++tag;	
	}
	while (!tag.isend());
}

//-----------------------------------------------------------------------------
void FEBioMeshDataSection::ParseMaterialData(XMLTag& tag, FEElementSet& set, const string& pname)
{
	vector<ELEMENT_DATA> data;
	ParseElementData(tag, set, data, 1);
	for (size_t i=0; i<data.size(); ++i)
	{
		ELEMENT_DATA& di = data[i];
		FEElement& el = *m_pelem[set[i]-1];

		for (int j=0; j<el.GaussPoints(); ++j)
		{
			FEMaterialPoint* pt = el.GetMaterialPoint(j);
			while (pt)
			{
				FEParameterList& pl = pt->GetParameterList();
				FEParam* p = pl.Find(pname.c_str());
				if (p) 
				{
					if ((p->dim() == 1) && (p->type() == FE_PARAM_DOUBLE))
					{
						p->value<double>() = di.val[0];
					}
					pt = 0;
				}
				else
				{
					pt = pt->Next();
					if (pt == 0) throw XMLReader::InvalidTag(tag);
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------
void FEBioMeshDataSection::ParseElementData(XMLTag& tag, FEElementSet& set, vector<ELEMENT_DATA>& values, int nvalues)
{
	// get the total nr of elements
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	int nelems = set.size();

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
