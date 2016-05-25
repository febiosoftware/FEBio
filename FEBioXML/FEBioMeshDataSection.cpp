#include "stdafx.h"
#include "FEBioMeshDataSection.h"
#include "FECore/FEModel.h"
#include "FECore/DOFS.h"

//-----------------------------------------------------------------------------
void FEBioMeshDataSection::Parse(XMLTag& tag)
{
	// make sure that the version is 2.5
	int nversion = m_pim->Version();
	if (nversion < 0x0205) throw XMLReader::InvalidTag(tag);

	// Make sure there is something in this tag
	if (tag.isleaf()) return;

	FEModel& fem = *GetFEModel();
	FEMesh& mesh = *m_pim->GetFEMesh();

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

	// loop over all mesh data
	++tag;
	do
	{
		if (tag == "ElementData")
		{
			const char* szset = tag.AttributeValue("elem_set");
			FEElementSet* part = mesh.FindElementSet(szset);
			if (part == 0) throw XMLReader::InvalidAttributeValue(tag, "elem_set", szset);

			const char* sztype = tag.AttributeValue("var");
			if      (strcmp(sztype, "shell thickness") == 0) ParseShellThickness(tag, *part);
			else if (strcmp(sztype, "fiber"          ) == 0) ParseMaterialFibers(tag, *part);
			else if (strcmp(sztype, "mat_axis"       ) == 0) ParseMaterialAxes  (tag, *part);
			else ParseMaterialData(tag, *part, sztype);
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
			m_pim->AddDataArray(szname, pdata);

			pdata->SetName(szname);
			pdata->Create(psurf);
			m_pim->ParseDataArray(tag, *pdata, "face");
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

			FEDataArray* pdata = new FEDataArray(dataType);
			m_pim->AddDataArray(szname, pdata);

			pdata->resize(nodeSet->size());
			m_pim->ParseDataArray(tag, *pdata, "node");
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
void set_element_fiber(FEElement& el, const vec3d& v);
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

			set_element_fiber(el, v);
		}
	}
}

//-----------------------------------------------------------------------------
void FEBioMeshDataSection::ParseMaterialAxes(XMLTag& tag, FEElementSet& set)
{
	// TODO: implement this (I can't use ParseElementData for this).
	assert(false);
	throw XMLReader::InvalidTag(tag);
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
	FEMesh& mesh = *m_pim->GetFEMesh();
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
