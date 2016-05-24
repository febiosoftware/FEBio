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

	// create the pelem array
	m_pelem.assign(nelems, static_cast<FEElement*>(0));
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
		if (tag == "elem_data")
		{
			const char* sztype = tag.AttributeValue("var");
			if      (strcmp(sztype, "shell thickness") == 0) ParseShellThickness(tag);
			else if (strcmp(sztype, "fiber"          ) == 0) ParseMaterialFibers(tag);
			else if (strcmp(sztype, "mat_axis"       ) == 0) ParseMaterialAxes  (tag);
			else ParseMaterialData(tag, sztype);
		}
		else if (tag == "surface_data")
		{
			const char* szsurf = tag.AttributeValue("surf");
			FEFacetSet* psurf = mesh.FindFacetSet(szsurf);
			if (psurf == 0) throw XMLReader::InvalidAttributeValue(tag, "surf", szsurf);

			const char* sztype = tag.AttributeValue("data_type", true);
			if (sztype == 0) sztype = "scalar";

			int dataType = -1;
			if      (strcmp(sztype, "scalar") == 0) dataType = FE_DOUBLE;
			else if (strcmp(sztype, "vec3"  ) == 0) dataType = FE_VEC3D;
			if (dataType == -1) throw XMLReader::InvalidAttributeValue(tag, "data_type", sztype);

			const char* szname = tag.AttributeValue("name");
			FESurfaceMap* pdata = new FESurfaceMap(dataType);
			m_pim->AddSurfaceMap(pdata);

			pdata->SetName(szname);
			pdata->Create(psurf);
			m_pim->ParseSurfaceMap(tag, *pdata);
		}
		else throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());
}

//-----------------------------------------------------------------------------
void FEBioMeshDataSection::ParseShellThickness(XMLTag& tag)
{
	vector<ELEMENT_DATA> data;
	ParseElementData(tag, data, FEElement::MAX_NODES);
	for (size_t i=0; i<data.size(); ++i)
	{
		ELEMENT_DATA& di = data[i];
		FEElement& el = *m_pelem[di.nid];

		if (el.Class() != FE_ELEM_SHELL) throw XMLReader::InvalidTag(tag);
		FEShellElement& shell = static_cast<FEShellElement&> (el);
		
		int ne = shell.Nodes();
		if (ne != di.nval) throw XMLReader::InvalidTag(tag);
		for (int j=0; j<ne; ++j) shell.m_h0[j] = di.val[j];
	}
}

//-----------------------------------------------------------------------------
// Defined in FEBioGeometrySection.cpp
void set_element_fiber(FEElement& el, const vec3d& v);
void set_element_mat_axis(FEElement& el, const vec3d& v1, const vec3d& v2);

//-----------------------------------------------------------------------------
void FEBioMeshDataSection::ParseMaterialFibers(XMLTag& tag)
{
	vector<ELEMENT_DATA> data;
	ParseElementData(tag, data, 3);
	for (size_t i=0; i<data.size(); ++i)
	{
		ELEMENT_DATA& di = data[i];
		FEElement& el = *m_pelem[di.nid];

		if (di.nval != 3) throw XMLReader::InvalidTag(tag);
		vec3d v(di.val[0], di.val[1], di.val[2]);

		set_element_fiber(el, v);
	}
}

//-----------------------------------------------------------------------------
void FEBioMeshDataSection::ParseMaterialAxes(XMLTag& tag)
{
	// TODO: implement this (I can't use ParseElementData for this).
	assert(false);
	throw XMLReader::InvalidTag(tag);
}

//-----------------------------------------------------------------------------
void FEBioMeshDataSection::ParseMaterialData(XMLTag& tag, const string& pname)
{
	vector<ELEMENT_DATA> data;
	ParseElementData(tag, data, 1);
	for (size_t i=0; i<data.size(); ++i)
	{
		ELEMENT_DATA& di = data[i];
		FEElement& el = *m_pelem[di.nid];

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
void FEBioMeshDataSection::ParseElementData(XMLTag& tag, vector<ELEMENT_DATA>& values, int nvalues)
{
	// get the total nr of elements
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = *m_pim->GetFEMesh();
	int nelems = mesh.Elements();

	// make sure we start with an empty array
	values.clear();
	values.reserve(nelems);

	ELEMENT_DATA data;
	++tag;
	do
	{
		if (tag == "elem")
		{
			// get the element number
			const char* szid = tag.AttributeValue("id");
			int n = atoi(szid)-1;

			// make sure the number is valid
			if ((n<0) || (n>=nelems)) throw XMLReader::InvalidAttributeValue(tag, "id", szid);

			data.nid = n;
			data.nval = tag.value(data.val, nvalues);
			values.push_back(data);
		}
		else if (tag == "elem_set")
		{
			const char* szname = tag.AttributeValue("elset");
			// find domain with this name
			FEElementSet* pset = mesh.FindElementSet(szname);
			if (pset == 0) throw XMLReader::InvalidAttributeValue(tag, "elset", szname);

			double d[FEElement::MAX_NODES]; 
			int nval = tag.value(d, nvalues);

			int n = pset->size();
			for (int i=0; i<n; ++i)
			{
				// get a pointer to the element
				int nid = (*pset)[i] - 1;
				
				data.nid = nid;
				data.nval = nval;
				for (int j=0; j<nval; ++j) data.val[j] = d[j];
				values.push_back(data);
			}
		}
		else throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());
}
