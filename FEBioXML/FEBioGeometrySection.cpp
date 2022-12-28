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
#include "FEBioGeometrySection.h"
#include "FECore/FESolidDomain.h"
#include "FECore/FEShellDomain.h"
#include "FECore/FETrussDomain.h"
#include "FECore/FEDomain2D.h"
#include "FECore/FEModel.h"
#include "FECore/FEMaterial.h"
#include "FECore/FECoreKernel.h"
#include <FECore/FENodeNodeList.h>

//-----------------------------------------------------------------------------
bool FEBioGeometrySection::ReadElement(XMLTag &tag, FEElement& el, int nid)
{
	el.SetID(nid);
	int n[FEElement::MAX_NODES];
	int m = tag.value(n, el.Nodes());
	if (m != el.Nodes()) return false;
	GetBuilder()->GlobalToLocalID(n, el.Nodes(), el.m_node);
	return true;
}

//=============================================================================
// FEBioGeometrySection1x
//============================================================================= 

//-----------------------------------------------------------------------------
void FEBioGeometrySection1x::Parse(XMLTag& tag)
{
	FEModelBuilder* feb = GetBuilder();
	feb->m_maxid = 0;

	++tag;
	do
	{
		if      (tag == "Nodes"      ) ParseNodeSection(tag);
		else if (tag == "Elements"   ) ParseElementSection(tag);
		else if (tag == "ElementData") ParseElementDataSection(tag);
		else throw XMLReader::InvalidTag(tag);
		++tag;
	} 
	while (!tag.isend());

	// At this point the mesh is completely read in.
	// Now we can allocate the degrees of freedom.
	// NOTE: We do this here since the mesh no longer automatically allocates the dofs.
	//       At some point I want to be able to read the mesh before deciding any physics.
	//       When that happens I'll have to move this elsewhere.
	FEModel& fem = *GetFEModel();
	int MAX_DOFS = fem.GetDOFS().GetTotalDOFS();
	fem.GetMesh().SetDOFS(MAX_DOFS);
}

//-----------------------------------------------------------------------------
//! Reads the Nodes section of the FEBio input file
void FEBioGeometrySection1x::ParseNodeSection(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	int N0 = mesh.Nodes();

	// get the largest nodal ID
	// (It is assumed that nodes are sorted by ID so the last node should have the largest ID)
	int max_id = 0;
	if (N0 > 0) max_id = mesh.Node(N0 - 1).GetID();

	// first we need to figure out how many nodes there are
	XMLTag t(tag);
	int nodes = tag.children();

	// see if this list defines a set
	const char* szl = tag.AttributeValue("set", true);
	FENodeSet* ps = 0;
	if (szl)
	{
		ps = new FENodeSet(&fem);

		ps->SetName(szl);
		mesh.AddNodeSet(ps);
	}

	// resize node's array
	mesh.AddNodes(nodes);

	// read nodal coordinates
	++tag;
	for (int i = 0; i<nodes; ++i)
	{
		FENode& node = mesh.Node(N0 + i);
		value(tag, node.m_r0);
		node.m_rt = node.m_r0;

		// get the nodal ID
		int nid = -1;
		tag.AttributeValue("id", nid);

		// Make sure it is valid
		if (nid <= max_id) throw XMLReader::InvalidAttributeValue(tag, "id");

		// set the ID
		node.SetID(nid);
		max_id = nid;

		// go on to the next node
		++tag;
	}

	// If a node set is defined add these nodes to the node-set
	if (ps)
	{
		for (int i = 0; i<nodes; ++i) ps->Add(N0 + i);
	}

	// tell the file reader to rebuild the node ID table
	GetBuilder()->BuildNodeList();
}

//-----------------------------------------------------------------------------
//! This function reads the Element section from the FEBio input file. It also
//! creates the domain classes which store the element data. A domain is defined
//! by the module (structural, poro, heat, etc), the element type (solid, shell,
//! etc.) and the material. 
//!
void FEBioGeometrySection1x::ParseElementSection(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// first we need to figure out how many elements 
	// and how many domains there are
	vector<FEDOMAIN> dom;
	vector<int>	ED; ED.reserve(1000);
	XMLTag t(tag); ++t;
	int elems = 0, i;
	while (!t.isend())
	{
		// get the material ID
		const char* szmat = t.AttributeValue("mat");
		int nmat = atoi(szmat) - 1;
		if ((nmat < 0) || (nmat >= fem.Materials())) throw FEBioImport::InvalidMaterial(elems + 1);

		// get the element type
		FE_Element_Spec spec = GetBuilder()->ElementSpec(t.Name());
		if (FEElementLibrary::IsValid(spec) == false) throw FEBioImport::InvalidElementType();

		// get the element ID
		int nid = -1;
		t.AttributeValue("id", nid);

		// keep track of the largest element ID
		if (nid > GetBuilder()->m_maxid) GetBuilder()->m_maxid = nid;

		// find a domain for this element
		int ndom = -1;
		for (i = 0; i<(int)dom.size(); ++i)
		{
			FEDOMAIN& d = dom[i];
			if ((d.mat == nmat) && (d.elem.eshape == spec.eshape))
			{
				ndom = i;
				d.nel++;
				break;
			}
		}
		if (ndom == -1)
		{
			FEDOMAIN d;
			d.mat = nmat;
			d.elem = spec;
			d.nel = 1;
			ndom = (int)dom.size();
			dom.push_back(d);
		}

		ED.push_back(ndom);
		elems++;
		++t;
	}

	FECoreKernel& febio = FECoreKernel::GetInstance();

	// create the domains
	for (i = 0; i<(int)dom.size(); ++i)
	{
		FEDOMAIN& d = dom[i];

		// get material class
		FEMaterial* pmat = fem.GetMaterial(d.mat);

		// create the new domain
		FEDomain* pdom = GetBuilder()->CreateDomain(d.elem, pmat);
		if (pdom == 0) throw FEBioImport::FailedCreatingDomain();

		// add it to the mesh
		assert(d.nel);
		pdom->Create(d.nel, d.elem);
		mesh.AddDomain(pdom);

		// we reset the nr of elements since we'll be using 
		// that variable is a counter in the next loop
		d.nel = 0;
	}

	// read element data
	++tag;
	int nid = 1;
	for (i = 0; i<elems; ++i, ++nid)
	{
		int nd = ED[i];
		int ne = dom[nd].nel++;

		// get the domain to which this element belongs
		FEDomain& domi = mesh.Domain(nd);

		// get the material ID
		int nmat = atoi(tag.AttributeValue("mat")) - 1;
		domi.SetMatID(nmat);

		// get the material class
		FEMaterial* pmat = fem.GetMaterial(nmat);
		assert(pmat == domi.GetMaterial());

		// determine element shape
		FE_Element_Spec espec = GetBuilder()->ElementSpec(tag.Name());
		if (FEElementLibrary::IsValid(espec) == false) throw FEBioImport::InvalidElementType();

		if (espec.etype != dom[ED[i]].elem.etype) throw XMLReader::InvalidTag(tag);
		if (ReadElement(tag, domi.ElementRef(ne), nid) == false) throw XMLReader::InvalidValue(tag);

		// go to next tag
		++tag;
	}

	// assign material point data
	for (i = 0; i<mesh.Domains(); ++i)
	{
		FEDomain& d = mesh.Domain(i);
		d.CreateMaterialPointData();
	}
}

//-----------------------------------------------------------------------------
void set_element_fiber(FEElement& el, const vec3d& v, int ncomp)
{
	// normalize fiber
	vec3d a = v;
	a.unit();

	// set up a orthonormal coordinate system
	vec3d b(0, 1, 0);
	if (fabs(fabs(a*b) - 1) < 1e-7) b = vec3d(0, 0, 1);
	vec3d c = a^b;
	b = c^a;

	// make sure they are unit vectors
	b.unit();
	c.unit();

	for (int i = 0; i<el.GaussPoints(); ++i)
	{
		FEMaterialPoint* mp = (ncomp == -1) ? el.GetMaterialPoint(i) : el.GetMaterialPoint(i)->GetPointData(ncomp);
/*
		FEElasticMaterialPoint& pt = *mp->ExtractData<FEElasticMaterialPoint>();
		mat3d& m = pt.m_Q;
		m.zero();
		m[0][0] = a.x; m[0][1] = b.x; m[0][2] = c.x;
		m[1][0] = a.y; m[1][1] = b.y; m[1][2] = c.y;
		m[2][0] = a.z; m[2][1] = b.z; m[2][2] = c.z;
*/
	}
}

//-----------------------------------------------------------------------------
void set_element_mat_axis(FEElement& el, const vec3d& v1, const vec3d& v2, int ncomp)
{
	vec3d a = v1;
	vec3d d = v2;

	vec3d c = a^d;
	vec3d b = c^a;

	// normalize
	a.unit();
	b.unit();
	c.unit();

	// assign to element
	for (int i = 0; i<el.GaussPoints(); ++i)
	{
        FEMaterialPoint* mp = (ncomp == -1) ? el.GetMaterialPoint(i) : el.GetMaterialPoint(i)->GetPointData(ncomp);

//		FEElasticMaterialPoint& pt = *mp->ExtractData<FEElasticMaterialPoint>();
//		pt.m_Q = mat3d(a, b, c);
	}
}

//-----------------------------------------------------------------------------
void FEBioGeometrySection1x::ParseElementData(FEElement& el, XMLTag& tag)
{
	vec3d a, d;
	if (tag == "fiber")
	{
		// read the fiber direction
		value(tag, a);
		set_element_fiber(el, a, 0);
	}
	else if (tag == "mat_axis")
	{
		++tag;
		do
		{
			if (tag == "a") value(tag, a);
			else if (tag == "d") value(tag, d);
			else throw XMLReader::InvalidTag(tag);

			++tag;
		} while (!tag.isend());
		set_element_mat_axis(el, a, d, -1);
	}
	else if (tag == "thickness")
	{
		if (el.Class() != FE_ELEM_SHELL) throw XMLReader::InvalidTag(tag);
		FEShellElement& shell = static_cast<FEShellElement&> (el);

		// read shell thickness
		tag.value(&shell.m_h0[0], shell.Nodes());
	}
	else if (tag == "area")
	{
		if (el.Class() != FE_ELEM_TRUSS) throw XMLReader::InvalidTag(tag);
		FETrussElement& truss = static_cast<FETrussElement&>(el);

		// read truss area
		value(tag, truss.m_a0);
	}
	else
	{
		for (int i = 0; i<el.GaussPoints(); ++i)
		{
			FEMaterialPoint* pt = el.GetMaterialPoint(i);
			// TODO: Material point parameters are no longer supported so I need to reimplement this
/*			while (pt)
			{
				FEParameterList& pl = pt->GetParameterList();
				if (ReadParameter(tag, pl)) break;

				bool tagFound = false;
				for (int i = 0; i<pt->Components(); ++i)
				{
					FEParameterList& pl = pt->GetPointData(i)->GetParameterList();
					if (ReadParameter(tag, pl))
					{
						tagFound = true;
						break;
					}
					
				}
				if (tagFound) break;

				pt = pt->Next();
				if (pt == 0) throw XMLReader::InvalidTag(tag);
			}
*/		}
	}
}

//-----------------------------------------------------------------------------
//! Reads the ElementData section from the FEBio input file

void FEBioGeometrySection1x::ParseElementDataSection(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// get the total nr of elements
	int nelems = mesh.Elements();

	//make sure we've read the element section
	if (nelems == 0) throw XMLReader::InvalidTag(tag);

	// create the pelem array
	vector<FEElement*> pelem;
	pelem.assign(nelems, static_cast<FEElement*>(0));

	for (int nd = 0; nd<mesh.Domains(); ++nd)
	{
		FEDomain& d = mesh.Domain(nd);
		for (int i = 0; i<d.Elements(); ++i)
		{
			FEElement& el = d.ElementRef(i);
			assert(pelem[el.GetID() - 1] == 0);
			pelem[el.GetID() - 1] = &el;
		}
	}

	// read additional element data
	++tag;
	do
	{
		// make sure this is an "element" tag
		if (tag == "element")
		{
			// get the element number
			const char* szid = tag.AttributeValue("id");
			int n = atoi(szid) - 1;

			// make sure the number is valid
			if ((n<0) || (n >= nelems)) throw XMLReader::InvalidAttributeValue(tag, "id", szid);

			// get a pointer to the element
			FEElement* pe = pelem[n];

			++tag;
			do
			{
				ParseElementData(*pe, tag);
				++tag;
			} while (!tag.isend());
		}
		else if (tag == "elset")
		{
			const char* szname = tag.AttributeValue("set");
			// find domain with this name
			FEElementSet* pset = mesh.FindElementSet(szname);
			if (pset == 0) throw XMLReader::InvalidAttributeValue(tag, "set", szname);

			++tag;
			do
			{
				int n = pset->Elements();
				for (int i = 0; i<n; ++i)
				{
					// get a pointer to the element
					int nid = (*pset)[i] - 1;
					FEElement* pe = pelem[nid];
					ParseElementData(*pe, tag);
				}
				++tag;
			} while (!tag.isend());
		}
		else throw XMLReader::InvalidTag(tag);

		++tag;
	} while (!tag.isend());
}

//-----------------------------------------------------------------------------
void FEBioGeometrySection1x::ParseNodeSetSection(XMLTag& tag)
{
	// read the node set
	FENodeSet* pns = GetFEBioImport()->ParseNodeSet(tag, "set");
	if (pns == 0) throw XMLReader::InvalidTag(tag);
}

//=============================================================================
// FEBioGeometrySection2
//============================================================================= 

//-----------------------------------------------------------------------------
void FEBioGeometrySection2::Parse(XMLTag& tag)
{
	FEModelBuilder* feb = GetBuilder();
	feb->m_maxid = 0;

	++tag;
	do
	{
		if      (tag == "Nodes"      ) ParseNodeSection       (tag);
		else if (tag == "Elements"   ) ParseElementSection    (tag);
		else if (tag == "NodeSet"    ) ParseNodeSetSection    (tag);
		else if (tag == "Surface"    ) ParseSurfaceSection    (tag);
		else if (tag == "Edge"       ) ParseEdgeSection       (tag);
		else if (tag == "ElementSet" ) ParseElementSetSection (tag);
		else if (tag == "ElementData") ParseElementDataSection(tag);
		else throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());

	// At this point the mesh is completely read in.
	// Now we can allocate the degrees of freedom.
	// NOTE: We do this here since the mesh no longer automatically allocates the dofs.
	//       At some point I want to be able to read the mesh before deciding any physics.
	//       When that happens I'll have to move this elsewhere.
	FEModel& fem = *GetFEModel();
	int MAX_DOFS = fem.GetDOFS().GetTotalDOFS();
	fem.GetMesh().SetDOFS(MAX_DOFS);
}

//-----------------------------------------------------------------------------
//! Reads the Nodes section of the FEBio input file
void FEBioGeometrySection2::ParseNodeSection(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	int N0 = mesh.Nodes();

	// get the largest nodal ID
	// (It is assumed that nodes are sorted by ID so the last node should have the largest ID)
	int max_id = 0;
	if (N0 > 0) max_id = mesh.Node(N0 - 1).GetID();

	// first we need to figure out how many nodes there are
	XMLTag t(tag);
	int nodes = tag.children();

	// see if this list defines a set
	const char* szl = tag.AttributeValue("set", true);
	FENodeSet* ps = 0;
	if (szl)
	{
		ps = new FENodeSet(&fem);

		ps->SetName(szl);
		mesh.AddNodeSet(ps);
	}

	// resize node's array
	mesh.AddNodes(nodes);

	// read nodal coordinates
	++tag;
	for (int i = 0; i<nodes; ++i)
	{
		FENode& node = mesh.Node(N0 + i);
		value(tag, node.m_r0);
		node.m_rt = node.m_r0;

		// get the nodal ID
		int nid = -1;
		tag.AttributeValue("id", nid);

		// Make sure it is valid
		if (nid <= max_id) throw XMLReader::InvalidAttributeValue(tag, "id");

		// set the ID
		node.SetID(nid);
		max_id = nid;

		// go on to the next node
		++tag;
	}

	// If a node set is defined add these nodes to the node-set
	if (ps)
	{
		for (int i = 0; i<nodes; ++i) ps->Add(N0 + i);
	}

	// tell the file reader to rebuild the node ID table
	GetBuilder()->BuildNodeList();
}

//-----------------------------------------------------------------------------
//! This function reads the Element section from the FEBio input file. It also
//! creates the domain classes which store the element data. A domain is defined
//! by the module (structural, poro, heat, etc), the element type (solid, shell,
//! etc.) and the material. 
//!
void FEBioGeometrySection2::ParseElementSection(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// get the material ID
	const char* szmat = tag.AttributeValue("mat");

	// as of 2.5, we allow materials to be referenced by name, so try this first
	FEMaterial* pmat = fem.FindMaterial(szmat);
	if (pmat == 0)
	{
		// if we didn't find the material we assume that the material tag uses an ID
		int nmat = atoi(szmat) - 1;
		if ((nmat < 0) || (nmat >= fem.Materials())) throw FEBioImport::InvalidDomainMaterial();

		// get the domain's material class
		pmat = fem.GetMaterial(nmat);
	}

	// get the name
	const char* szname = tag.AttributeValue("elset", true);
	if (szname == 0) szname = "_unnamed";

	// get the element type
	const char* sztype = tag.AttributeValue("type");
	FE_Element_Spec espec = GetBuilder()->ElementSpec(sztype);
	if (FEElementLibrary::IsValid(espec) == false) throw FEBioImport::InvalidElementType();

	// create the new domain
	FEDomain* pdom = GetBuilder()->CreateDomain(espec, pmat);
	if (pdom == 0) throw FEBioImport::FailedCreatingDomain();
	FEDomain& dom = *pdom;
	dom.SetName(szname);

	// count elements
	int elems = tag.children();
	assert(elems);

	// add domain it to the mesh
	pdom->Create(elems, espec);
	pdom->SetMatID(pmat->GetID() - 1);
	mesh.AddDomain(pdom);

	// for named domains, we'll also create an element set
	FEElementSet* pg = 0;
	if (szname)
	{
		pg = new FEElementSet(&fem);
		pg->SetName(szname);
		mesh.AddElementSet(pg);
	}

	// read element data
	++tag;
	for (int i = 0; i<elems; ++i)
	{
		if ((tag == "elem") == false) throw XMLReader::InvalidTag(tag);

		// get the element ID
		int nid;
		tag.AttributeValue("id", nid);

		// Make sure element IDs increase
		//		if (nid <= m_pim->m_maxid) throw XMLReader::InvalidAttributeValue(tag, "id");

		// keep track of the largest element ID
		// (which by assumption is the ID that was just read in)
		GetBuilder()->m_maxid = nid;

		// read the element data
		if (ReadElement(tag, dom.ElementRef(i), nid) == false) throw XMLReader::InvalidValue(tag);

		// go to next tag
		++tag;
	}

	// create the element set
	if (pg) pg->Create(pdom);

	// assign material point data
	dom.CreateMaterialPointData();
}

//-----------------------------------------------------------------------------
void FEBioGeometrySection2::ParseElementData(FEElement& el, XMLTag& tag)
{
	vec3d a, d;
	if (tag == "fiber")
	{
		// read the fiber direction
		value(tag, a);
		set_element_fiber(el, a, -1);
	}
	else if (tag == "mat_axis")
	{
		++tag;
		do
		{
			if (tag == "a") value(tag, a);
			else if (tag == "d") value(tag, d);
			else throw XMLReader::InvalidTag(tag);

			++tag;
		} while (!tag.isend());
		set_element_mat_axis(el, a, d, -1);
	}
	else if (tag == "thickness")
	{
		if (el.Class() != FE_ELEM_SHELL) throw XMLReader::InvalidTag(tag);
		FEShellElement& shell = static_cast<FEShellElement&> (el);

		// read shell thickness
		tag.value(&shell.m_h0[0], shell.Nodes());
	}
	else if (tag == "area")
	{
		if (el.Class() != FE_ELEM_TRUSS) throw XMLReader::InvalidTag(tag);
		FETrussElement& truss = static_cast<FETrussElement&>(el);

		// read truss area
		value(tag, truss.m_a0);
	}
	else
	{
		for (int i = 0; i<el.GaussPoints(); ++i)
		{
			FEMaterialPoint* pt = el.GetMaterialPoint(i);
			// TODO: material point parameters are no longer supported so I need to reimplement this.
/*			while (pt)
			{
				FEParameterList& pl = pt->GetParameterList();
				if (ReadParameter(tag, pl)) break;

				bool tagFound = false;
				if (pt->Components() > 1)
				{
					for (int i = 0; i<pt->Components(); ++i)
					{
						FEParameterList& pl = pt->GetPointData(i)->GetParameterList();
						if (ReadParameter(tag, pl))
						{
							tagFound = true;
							break;
						}
					}
				}

				if (tagFound) break;

				pt = pt->Next();
				if (pt == 0) throw XMLReader::InvalidTag(tag);
			}
*/		}
	}
}

//-----------------------------------------------------------------------------
//! Reads the ElementData section from the FEBio input file

void FEBioGeometrySection2::ParseElementDataSection(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// get the total nr of elements
	int nelems = mesh.Elements();

	//make sure we've read the element section
	if (nelems == 0) throw XMLReader::InvalidTag(tag);

	// create the pelem array
	vector<FEElement*> pelem;
	pelem.assign(nelems, static_cast<FEElement*>(0));

	for (int nd = 0; nd<mesh.Domains(); ++nd)
	{
		FEDomain& d = mesh.Domain(nd);
		for (int i = 0; i<d.Elements(); ++i)
		{
			FEElement& el = d.ElementRef(i);
			assert(pelem[el.GetID() - 1] == 0);
			pelem[el.GetID() - 1] = &el;
		}
	}

	// read additional element data
	++tag;
	do
	{
		// make sure this is an "element" tag
		if (tag == "element")
		{
			// get the element number
			const char* szid = tag.AttributeValue("id");
			int n = atoi(szid) - 1;

			// make sure the number is valid
			if ((n<0) || (n >= nelems)) throw XMLReader::InvalidAttributeValue(tag, "id", szid);

			// get a pointer to the element
			FEElement* pe = pelem[n];

			++tag;
			do
			{
				ParseElementData(*pe, tag);
				++tag;
			} while (!tag.isend());
		}
		else if (tag == "elset")
		{
			const char* szname = tag.AttributeValue("set");
			// find domain with this name
			FEElementSet* pset = mesh.FindElementSet(szname);
			if (pset == 0) throw XMLReader::InvalidAttributeValue(tag, "set", szname);

			++tag;
			do
			{
				int n = pset->Elements();
				for (int i = 0; i<n; ++i)
				{
					// get a pointer to the element
					int nid = (*pset)[i] - 1;
					FEElement* pe = pelem[nid];
					ParseElementData(*pe, tag);
				}
				++tag;
			} while (!tag.isend());
		}
		else throw XMLReader::InvalidTag(tag);

		++tag;
	} while (!tag.isend());
}

//-----------------------------------------------------------------------------
void FEBioGeometrySection2::ParseNodeSetSection(XMLTag& tag)
{
	// read the node set
	FENodeSet* pns = GetFEBioImport()->ParseNodeSet(tag, "set");
	if (pns == 0) throw XMLReader::InvalidTag(tag);
}

//-----------------------------------------------------------------------------
//! Reads a Geometry\Edge section.
void FEBioGeometrySection2::ParseEdgeSection(XMLTag& tag)
{
	// get the mesh
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// get the number of nodes
	// (we use this for checking the node indices of the facets)
	int NN = mesh.Nodes();

	// get the required name attribute
	const char* szname = tag.AttributeValue("name");

	// count nr of segments
	int nsegs = tag.children();

	// allocate storage for segments
	FESegmentSet* ps = new FESegmentSet(&fem);
	ps->Create(nsegs);
	ps->SetName(szname);

	// add it to the mesh
	mesh.AddSegmentSet(ps);

	// read segments
	++tag;
	int nf[FEElement::MAX_NODES];
	for (int i = 0; i<nsegs; ++i)
	{
		FESegmentSet::SEGMENT& line = ps->Segment(i);

		// set the facet type
		if (tag == "line2") line.ntype = 2;
		else throw XMLReader::InvalidTag(tag);

		// we assume that the segment type also defines the number of nodes
		int N = line.ntype;
		tag.value(nf, N);
		for (int j = 0; j<N; ++j)
		{
			int nid = nf[j] - 1;
			if ((nid<0) || (nid >= NN)) throw XMLReader::InvalidValue(tag);
			line.node[j] = nid;
		}

		++tag;
	}
}

//-----------------------------------------------------------------------------
//! Reads a Geometry\Surface section.
void FEBioGeometrySection2::ParseSurfaceSection(XMLTag& tag)
{
	// get the mesh
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// get the number of nodes
	// (we use this for checking the node indices of the facets)
	int NN = mesh.Nodes();

	// get the required name attribute
	const char* szname = tag.AttributeValue("name");

	// count nr of faces
	int faces = tag.children();

	// allocate storage for faces
	FEFacetSet* ps = new FEFacetSet(&fem);
	ps->Create(faces);
	ps->SetName(szname);

	// add it to the mesh
	mesh.AddFacetSet(ps);

	// read faces
	++tag;
	int nf[FEElement::MAX_NODES];
	for (int i = 0; i<faces; ++i)
	{
		FEFacetSet::FACET& face = ps->Face(i);

		// set the facet type
		if (tag == "quad4") face.ntype = 4;
		else if (tag == "tri3") face.ntype = 3;
		else if (tag == "tri6") face.ntype = 6;
		else if (tag == "tri7") face.ntype = 7;
		else if (tag == "quad8") face.ntype = 8;
		else if (tag == "quad9") face.ntype = 9;
		else if (tag == "tri10") face.ntype = 10;
		else throw XMLReader::InvalidTag(tag);

		// we assume that the facet type also defines the number of nodes
		int N = face.ntype;
		tag.value(nf, N);
		for (int j = 0; j<N; ++j)
		{
			int nid = nf[j] - 1;
			if ((nid<0) || (nid >= NN)) throw XMLReader::InvalidValue(tag);
			face.node[j] = nid;
		}

		++tag;
	}
}

//-----------------------------------------------------------------------------
void FEBioGeometrySection2::ParseElementSetSection(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// get the name attribute
	const char* szname = tag.AttributeValue("name");

	// create a new element set
	FEElementSet* pg = new FEElementSet(&fem);
	pg->SetName(szname);

	vector<int> l;
	++tag;
	do
	{
		if (tag == "elem")
		{
			int nid = -1;
			tag.AttributeValue("id", nid);
			l.push_back(nid);
		}
		else { delete pg; throw XMLReader::InvalidTag(tag); }
		++tag;
	} while (!tag.isend());

	// only add non-empty element sets
	if (l.empty() == false)
	{
		// assign indices to element set
		pg->Create(l);

		// add the element set to the mesh
		mesh.AddElementSet(pg);
	}
	else delete pg;
}

//=============================================================================
// FEBioGeometrySection25
//============================================================================= 

//-----------------------------------------------------------------------------
void FEBioGeometrySection25::Parse(XMLTag& tag)
{
	FEModelBuilder* feb = GetBuilder();
	feb->m_maxid = 0;

	// read all sections
	++tag;
	do
	{
		if      (tag == "Nodes"      ) ParseNodeSection       (tag);
		else if (tag == "Elements"   ) ParseElementSection    (tag);
		else if (tag == "NodeSet"    ) ParseNodeSetSection    (tag);
		else if (tag == "Surface"    ) ParseSurfaceSection    (tag);
		else if (tag == "Edge"       ) ParseEdgeSection       (tag);
		else if (tag == "ElementSet" ) ParseElementSetSection (tag);
		else if (tag == "DiscreteSet") ParseDiscreteSetSection(tag);
		else if (tag == "SurfacePair") ParseSurfacePairSection(tag);
		else if (tag == "NodeSetPair") ParseNodeSetPairSection(tag);
		else if (tag == "NodeSetSet" ) ParseNodeSetSetSection (tag);
		else if (tag == "Part"       ) ParsePartSection       (tag);
		else if (tag == "Instance"   ) ParseInstanceSection   (tag);
		else throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());

	// At this point the mesh is completely read in.
	// Now we can allocate the degrees of freedom.
	// NOTE: We do this here since the mesh no longer automatically allocates the dofs.
	//       At some point I want to be able to read the mesh before deciding any physics.
	//       When that happens I'll have to move this elsewhere.
	FEModel& fem = *GetFEModel();
	int MAX_DOFS = fem.GetDOFS().GetTotalDOFS();
	fem.GetMesh().SetDOFS(MAX_DOFS);
}

//-----------------------------------------------------------------------------
void FEBioGeometrySection25::ParsePartSection(XMLTag& tag)
{
	// get the part name
	const char* szname = tag.AttributeValue("name");

	// Create a new part
	FEBModel::Part* part = m_feb.AddPart(szname);

	// get the type attribute
	const char* sztype = tag.AttributeValue("type", true);
	if (sztype == 0) sztype = "template";

	// check the value
	bool binstance = false;
	if      (strcmp(sztype, "instance") == 0) binstance = true;
	else if (strcmp(sztype, "template") == 0) binstance = false;
	else throw XMLReader::InvalidAttributeValue(tag, "type", sztype);

	// see if the from attribute is defined
	const char* szfrom = tag.AttributeValue("from", true);
	if (szfrom)
	{
		// make sure this is an empty leaf
		if (tag.isempty() == false) throw XMLReader::InvalidValue(tag);

		// redirect input to another file
		char xpath[256] = {0};
		sprintf(xpath, "febio_spec/Geometry/Part[@name=%s]", szname);
		XMLReader xml;
		if (xml.Open(szfrom) == false) throw XMLReader::InvalidAttributeValue(tag, "from", szfrom);
		XMLTag tag2;
		if (xml.FindTag(xpath, tag2) == false) throw XMLReader::InvalidAttributeValue(tag, "from", szfrom);

		// read the part
		ParsePart(tag2, part);
	}
	else ParsePart(tag, part);

	// instantiate the part
	if (binstance) 
	{
		if (m_feb.BuildPart(*GetFEModel(), *part) == false) throw FEBioImport::FailedBuildingPart(part->Name());
	}
}

//-----------------------------------------------------------------------------
void FEBioGeometrySection25::ParsePart(XMLTag& tag, FEBModel::Part* part)
{
	// read the part's sections
	++tag;
	do
	{
		if      (tag == "Nodes"      ) ParsePartNodeSection   (tag, part);
		else if (tag == "Elements"   ) ParsePartElementSection(tag, part);
		else if (tag == "NodeSet"    ) ParsePartNodeSetSection(tag, part);
		else if (tag == "Surface"    ) ParsePartSurfaceSection(tag, part);
		else throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());
}

//-----------------------------------------------------------------------------
void FEBioGeometrySection25::ParseInstanceSection(XMLTag& tag)
{
	// get the name and part tags
	const char* szname = tag.AttributeValue("name", true);
	const char* szpart = tag.AttributeValue("part");
	if (szname == 0) szname = szpart;

	// find the part
	FEBModel::Part* oldPart = m_feb.FindPart(szpart);
	if (oldPart == 0) throw XMLReader::InvalidAttributeValue(tag, "part", szpart);

	// copy the part
	FEBModel::Part* newPart = new FEBModel::Part(*oldPart);
	m_feb.AddPart(newPart);

	// rename the part
	newPart->SetName(szname);

	// parse any child tags
	Transform transform;
	if (tag.isleaf() == false)
	{
		++tag;
		do
		{
			if (tag == "translate")
			{
				double r[3];
				tag.value(r, 3);
				transform.SetPosition(vec3d(r[0], r[1], r[2]));
			}
			else if (tag == "rotate")
			{
				const char* sztype = tag.AttributeValue("type", true);
				if (sztype == 0) sztype = "quaternion";

				if (strcmp(sztype, "quaternion") == 0)
				{
					double v[4];
					tag.value(v, 4);
					quatd q(v[0], v[1], v[2], v[3]);
					transform.SetRotation(q);
				}
				else if (strcmp(sztype, "vector") == 0)
				{
					double v[3];
					tag.value(v, 3);
					transform.SetRotation(vec3d(v[0], v[1], v[2]));
				}
				else if (strcmp(sztype, "Euler") == 0)
				{
					double v[3];
					tag.value(v, 3);
					transform.SetRotation(v[0], v[1], v[2]);
				}
				else throw XMLReader::InvalidAttributeValue(tag, "type", sztype);
			}
			else if (tag == "scale")
			{
				double s[3];
				tag.value(s, 3);
				transform.SetScale(s[0], s[1], s[2]);
			}
			else if (tag == "Elements")
			{
				const char* szdom = tag.AttributeValue("name");
				const char* szmat = tag.AttributeValue("mat");
				
				FEBModel::Domain* dom = newPart->FindDomain(szdom);
				if (dom == 0) throw XMLReader::InvalidAttributeValue(tag, "name", szdom);

				// make sure the material is defined
				FEModel& fem = *GetFEModel();
				FEMaterial* mat = fem.FindMaterial(szmat);
				if (mat == 0) throw FEBioImport::InvalidDomainMaterial();

				dom->SetMaterialName(szmat);
			}
			else throw XMLReader::InvalidTag(tag);
			++tag;
		}
		while (!tag.isend());
	}

	// build this part
	if (m_feb.BuildPart(*GetFEModel(), *newPart, true, transform) == false) throw FEBioImport::FailedBuildingPart(newPart->Name());
}

//-----------------------------------------------------------------------------
//! Reads the Nodes section of the FEBio input file
void FEBioGeometrySection25::ParseNodeSection(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	int N0 = mesh.Nodes();

	// get the largest nodal ID
	// (It is assumed that nodes are sorted by ID so the last node should have the largest ID)
	int max_id = 0;
	if (N0 > 0) max_id = mesh.Node(N0 - 1).GetID();

	// first we need to figure out how many nodes there are
	XMLTag t(tag);
	int nodes = tag.children();

	// see if this list defines a set
	const char* szl = tag.AttributeValue("name", true);
	FENodeSet* ps = 0;
	if (szl)
	{
		ps = new FENodeSet(&fem);

		ps->SetName(szl);
		mesh.AddNodeSet(ps);
	}

	// resize node's array
	mesh.AddNodes(nodes);

	// read nodal coordinates
	++tag;
	for (int i = 0; i<nodes; ++i)
	{
		FENode& node = mesh.Node(N0 + i);
		value(tag, node.m_r0);
		node.m_rt = node.m_r0;

		// get the nodal ID
		int nid = -1;
		tag.AttributeValue("id", nid);

		// Make sure it is valid
		if (nid <= max_id) throw XMLReader::InvalidAttributeValue(tag, "id");

		// set the ID
		node.SetID(nid);
		max_id = nid;

		// go on to the next node
		++tag;
	}

	// If a node set is defined add these nodes to the node-set
	if (ps)
	{
		for (int i = 0; i<nodes; ++i) ps->Add(N0 + i);
	}

	// tell the file reader to rebuild the node ID table
	GetBuilder()->BuildNodeList();
}

//-----------------------------------------------------------------------------
//! Reads the Nodes section of the FEBio input file
void FEBioGeometrySection25::ParsePartNodeSection(XMLTag& tag, FEBModel::Part* part)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	int N0 = mesh.Nodes();

	// first we need to figure out how many nodes there are
	XMLTag t(tag);
	int nodes = tag.children();

	// see if this list defines a set
	const char* szname = tag.AttributeValue("name", true);
	FEBModel::NodeSet* ps = 0;
	if (szname)
	{
		ps = new FEBModel::NodeSet(szname);
		part->AddNodeSet(ps);
	}

	// allocate node
	vector<FEBModel::NODE> node(nodes);
	vector<int> nodeList(nodes);

	// read nodal coordinates
	++tag;
	for (int i = 0; i<nodes; ++i)
	{
		FEBModel::NODE& nd = node[i];
		value(tag, nd.r);

		// get the nodal ID
		tag.AttributeValue("id", nd.id);
		nodeList[i] = nd.id;

		// go on to the next node
		++tag;
	}

	// add nodes to the part
	part->AddNodes(node);

	// If a node set is defined add these nodes to the node-set
	if (ps) ps->SetNodeList(nodeList);
}





//-----------------------------------------------------------------------------
//! This function reads the Element section from the FEBio input file. It also
//! creates the domain classes which store the element data. A domain is defined
//! by the module (structural, poro, heat, etc), the element type (solid, shell,
//! etc.) and the material. 
//!
void FEBioGeometrySection25::ParseElementSection(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// get the material ID
	const char* szmat = tag.AttributeValue("mat");

	// as of 2.5, we allow materials to be referenced by name, so try this first
	FEMaterial* pmat = fem.FindMaterial(szmat);
	if (pmat == 0)
	{
		// if we didn't find the material we assume that the material tag uses an ID
		int nmat = atoi(szmat) - 1;
		if ((nmat < 0) || (nmat >= fem.Materials())) throw FEBioImport::InvalidDomainMaterial();

		// get the domain's material class
		pmat = fem.GetMaterial(nmat);
	}

	// get the name
	const char* szname = tag.AttributeValue("name", true);
	if (szname == 0) szname = "_unnamed";

	// get the element type
	const char* sztype = tag.AttributeValue("type");
	FE_Element_Spec espec = GetBuilder()->ElementSpec(sztype);
	if (FEElementLibrary::IsValid(espec) == false) throw FEBioImport::InvalidElementType();

	// create the new domain
	FEDomain* pdom = GetBuilder()->CreateDomain(espec, pmat);
	if (pdom == 0) throw FEBioImport::FailedCreatingDomain();
	FEDomain& dom = *pdom;
	dom.SetName(szname);

	// active flag
	const char* szactive = tag.AttributeValue("active", true);
	if (szactive)
	{
		if (strcmp(szactive, "false") == 0) pdom->SetActive(false);
	}

	// count elements
	vector<FEModelBuilder::ELEMENT> elemList; elemList.reserve(512000);
	++tag;
	do
	{
		if ((tag == "elem") == false) throw XMLReader::InvalidTag(tag);

		// get the element ID
		FEModelBuilder::ELEMENT el;
		tag.AttributeValue("id", el.nid);

		el.nodes = tag.value(el.node, FEElement::MAX_NODES);
		elemList.push_back(el);
		++tag;
	}
	while (!tag.isend());

	int elems = (int) elemList.size();
	assert(elems);

	// add domain it to the mesh
	pdom->Create(elems, espec);
	pdom->SetMatID(pmat->GetID() - 1);
	mesh.AddDomain(pdom);

	// for named domains, we'll also create an element set
	FEElementSet* pg = 0;
	if (szname)
	{
		pg = new FEElementSet(&fem);
		pg->SetName(szname);
		mesh.AddElementSet(pg);
	}

	// read element data
	for (int i = 0; i<elems; ++i)
	{
		FEModelBuilder::ELEMENT& elem = elemList[i];
		
		// get the element ID
		int nid = elem.nid;

		// Make sure element IDs increase
		//		if (nid <= m_pim->m_maxid) throw XMLReader::InvalidAttributeValue(tag, "id");

		// keep track of the largest element ID
		// (which by assumption is the ID that was just read in)
		GetBuilder()->m_maxid = nid;

		// process the element data
		FEElement& el = dom.ElementRef(i);
		el.SetID(nid);
		if (elem.nodes != el.Nodes()) throw XMLReader::InvalidTag(tag);
		GetBuilder()->GlobalToLocalID(elem.node, el.Nodes(), el.m_node);
	}

	// create the element set
	if (pg) pg->Create(pdom);

	// assign material point data
	dom.CreateMaterialPointData();
}


//-----------------------------------------------------------------------------
//! This function reads the Element section from the FEBio input file. It also
//! creates the domain classes which store the element data. A domain is defined
//! by the module (structural, poro, heat, etc), the element type (solid, shell,
//! etc.) and the material. 
//!
void FEBioGeometrySection25::ParsePartElementSection(XMLTag& tag, FEBModel::Part* part)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// get the material ID
	const char* szmat = tag.AttributeValue("mat", true);

	// get the (optional) name
	const char* szname = tag.AttributeValue("name", true);

	// get the element type
	const char* sztype = tag.AttributeValue("type");
	FE_Element_Spec espec = GetBuilder()->ElementSpec(sztype);
	if (FEElementLibrary::IsValid(espec) == false) throw FEBioImport::InvalidElementType();

	// create the new domain
	FEBModel::Domain* dom = new FEBModel::Domain(espec);
	if (szname) dom->SetName(szname);
	if (szmat) dom->SetMaterialName(szmat);

	// count elements
	int elems = tag.children();
	assert(elems);

	// add domain it to the mesh
	dom->Create(elems);
	part->AddDomain(dom);

	// for named domains, we'll also create an element set
	FEBModel::ElementSet* pg = 0;
	if (szname)
	{
		FEBModel::ElementSet* pg = new FEBModel::ElementSet(szname);
		part->AddElementSet(pg);
	}

	vector<int> elemList(elems);

	// read element data
	++tag;
	for (int i = 0; i<elems; ++i)
	{
		if ((tag == "elem") == false) throw XMLReader::InvalidTag(tag);

		FEBModel::ELEMENT& el = dom->GetElement(i);

		// get the element ID
		tag.AttributeValue("id", el.id);
		elemList[i] = el.id;

		// read the element data
		tag.value(el.node, FEElement::MAX_NODES);

		// go to next tag
		++tag;
	}

	// set the element list
	if (pg) pg->SetElementList(elemList);
}

//-----------------------------------------------------------------------------
void FEBioGeometrySection25::ParseNodeSetSection(XMLTag& tag)
{
	// read the node set
	FENodeSet* pns = GetFEBioImport()->ParseNodeSet(tag, "node_set");
	if (pns == 0) throw XMLReader::InvalidTag(tag);
}

//-----------------------------------------------------------------------------
//! Reads the Geometry::Groups section of the FEBio input file
void FEBioGeometrySection25::ParsePartNodeSetSection(XMLTag& tag, FEBModel::Part* part)
{
	const char* szname = tag.AttributeValue("name");

	// create the node set
	FEBModel::NodeSet* set = new FEBModel::NodeSet(szname);
	part->AddNodeSet(set);

	int nodes = tag.children();
	vector<int> nodeList(nodes);

	++tag;
	for (int i=0; i<nodes; ++i)
	{
		if (tag == "node")
		{
			// get the ID
			int nid;
			tag.AttributeValue("id", nid);
			nodeList[i] = nid;
		}
		else throw XMLReader::InvalidTag(tag);
		
		++tag;
	}

	// add nodes to the list
	set->SetNodeList(nodeList);
}

//-----------------------------------------------------------------------------
void FEBioGeometrySection25::ParseDiscreteSetSection(XMLTag& tag)
{
	// get the mesh
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// get the name
	const char* szname = tag.AttributeValue("name");

	// create the discrete element set
	FEDiscreteSet* ps = new FEDiscreteSet(&mesh);
	ps->SetName(szname);
	mesh.AddDiscreteSet(ps);

	// see if the generate attribute was defined
	const char* szgen = tag.AttributeValue("generate", true);
	if (szgen == 0)
	{
		// read the node pairs
		++tag;
		do
		{
			if (tag == "delem")
			{
				int n[2];
				tag.value(n, 2);
				n[0] -= 1; n[1] -= 1;
				ps->add(n[0], n[1]);
			}
			else throw XMLReader::InvalidTag(tag);
			++tag;
		}
		while (!tag.isend());
	}
	else
	{
		assert(tag.isempty());

		// generate all the springs from the mesh
		FENodeNodeList NNL;
		NNL.Create(mesh);

		// generate the springs
		for (int i = 0; i<NNL.Size(); ++i)
		{
			int nval = NNL.Valence(i);
			for (int j = 0; j<nval; ++j)
			{
				int n0 = i;
				int nj = NNL.NodeList(i)[j];
				if (n0 < nj)
				{
					ps->add(n0, nj);
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! Reads a Geometry\Edge section.
void FEBioGeometrySection25::ParseEdgeSection(XMLTag& tag)
{
	// get the mesh
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// get the number of nodes
	// (we use this for checking the node indices of the facets)
	int NN = mesh.Nodes();

	// get the required name attribute
	const char* szname = tag.AttributeValue("name");

	// count nr of segments
	int nsegs = tag.children();

	// allocate storage for segments
	FESegmentSet* ps = new FESegmentSet(&fem);
	ps->Create(nsegs);
	ps->SetName(szname);

	// add it to the mesh
	mesh.AddSegmentSet(ps);

	// read segments
	++tag;
	int nf[FEElement::MAX_NODES];
	for (int i = 0; i<nsegs; ++i)
	{
		FESegmentSet::SEGMENT& line = ps->Segment(i);

		// set the facet type
		if (tag == "line2") line.ntype = 2;
		else throw XMLReader::InvalidTag(tag);

		// we assume that the segment type also defines the number of nodes
		int N = line.ntype;
		tag.value(nf, N);
		for (int j = 0; j<N; ++j)
		{
			int nid = nf[j] - 1;
			if ((nid<0) || (nid >= NN)) throw XMLReader::InvalidValue(tag);
			line.node[j] = nid;
		}

		++tag;
	}
}

//-----------------------------------------------------------------------------
void FEBioGeometrySection25::ParseSurfacePairSection(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// create new surface pair
	FESurfacePair* p = new FESurfacePair(&mesh);

	// set its name
	const char* szname = tag.AttributeValue("name");
	p->SetName(szname);

	++tag;
	do
	{
		if (tag == "master")
		{
			const char* sz = tag.AttributeValue("surface");
			p->SetSecondarySurface(mesh.FindFacetSet(sz));
			if (p->GetSecondarySurface() == 0) throw XMLReader::InvalidAttributeValue(tag, "surface", sz);
		}
		else if (tag == "slave")
		{
			const char* sz = tag.AttributeValue("surface");
			p->SetPrimarySurface(mesh.FindFacetSet(sz));
			if (p->GetPrimarySurface() == 0) throw XMLReader::InvalidAttributeValue(tag, "surface", sz);
		}
		else throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());

	// add it to the mesh
	mesh.AddSurfacePair(p);
}

//-----------------------------------------------------------------------------
void FEBioGeometrySection25::ParseNodeSetPairSection(XMLTag& tag)
{
	FEModelBuilder::NodeSetPair p;
	const char* szname = tag.AttributeValue("name");
	strcpy(p.szname, szname);

	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	++tag;
	do
	{
		if (tag == "master")
		{
			const char* sz = tag.AttributeValue("node_set");
			p.set2 = mesh.FindNodeSet(sz);
			if (p.set2 == 0) throw XMLReader::InvalidAttributeValue(tag, "node_set", sz);
		}
		else if (tag == "slave")
		{
			const char* sz = tag.AttributeValue("node_set");
			p.set1 = mesh.FindNodeSet(sz);
			if (p.set1 == 0) throw XMLReader::InvalidAttributeValue(tag, "node_set", sz);
		}
		else throw XMLReader::InvalidTag(tag);
		++tag;
	} while (!tag.isend());

	GetBuilder()->AddNodeSetPair(p);
}

//-----------------------------------------------------------------------------
void FEBioGeometrySection25::ParseNodeSetSetSection(XMLTag& tag)
{
	FEModelBuilder::NodeSetSet p;
	const char* szname = tag.AttributeValue("name");
	strcpy(p.szname, szname);

	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	++tag;
	do
	{
		if (tag == "node_set")
		{
			const char* sz = tag.AttributeValue("node_set");
			FENodeSet* ns = mesh.FindNodeSet(sz);
			if (ns == 0) throw XMLReader::InvalidAttributeValue(tag, "node_set", sz);
			p.add(ns);
		}
		else throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());

	GetBuilder()->AddNodeSetSet(p);
}

//-----------------------------------------------------------------------------
//! Reads a Geometry\Surface section.
void FEBioGeometrySection25::ParseSurfaceSection(XMLTag& tag)
{
	// get the mesh
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// get the number of nodes
	// (we use this for checking the node indices of the facets)
	int NN = mesh.Nodes();

	// get the required name attribute
	const char* szname = tag.AttributeValue("name");

	// if parts are defined we use the new format
	if (m_feb.Parts() > 0)
	{
		FEFacetSet* ps = new FEFacetSet(&fem);
		ps->SetName(szname);

		// add it to the mesh
		mesh.AddFacetSet(ps);

		// read the child surfaces
		++tag;
		do
		{
			if (tag == "Surface")
			{
				const char* szatt = tag.AttributeValue("surface");
				FEFacetSet* pf = mesh.FindFacetSet(szatt);
				if (pf == 0) throw XMLReader::InvalidAttributeValue(tag, "surface", szatt);

				ps->Add(pf);
			}
			else throw XMLReader::InvalidTag(tag);
			++tag;
		} while (!tag.isend());

	}
	else
	{
		// count nr of faces
		int faces = tag.children();

		// allocate storage for faces
		FEFacetSet* ps = new FEFacetSet(&fem);
		ps->Create(faces);
		ps->SetName(szname);

		// add it to the mesh
		mesh.AddFacetSet(ps);

		// read faces
		++tag;
		int nf[FEElement::MAX_NODES];
		for (int i = 0; i<faces; ++i)
		{
			FEFacetSet::FACET& face = ps->Face(i);

			// set the facet type
			if (tag == "quad4") face.ntype = 4;
			else if (tag == "tri3") face.ntype = 3;
			else if (tag == "tri6") face.ntype = 6;
			else if (tag == "tri7") face.ntype = 7;
			else if (tag == "quad8") face.ntype = 8;
			else if (tag == "quad9") face.ntype = 9;
			else if (tag == "tri10") face.ntype = 10;
			else throw XMLReader::InvalidTag(tag);

			// we assume that the facet type also defines the number of nodes
			int N = face.ntype;
			int nread = tag.value(nf, N);
			if (nread != face.ntype) throw XMLReader::InvalidValue(tag);
			for (int j = 0; j<N; ++j)
			{
				int nid = nf[j] - 1;
				if ((nid<0) || (nid >= NN)) throw XMLReader::InvalidValue(tag);
				face.node[j] = nid;
			}

			++tag;
		}
	}
}

//-----------------------------------------------------------------------------
//! Reads a Geometry\Surface section.
void FEBioGeometrySection25::ParsePartSurfaceSection(XMLTag& tag, FEBModel::Part* part)
{
	// get the mesh
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// get the required name attribute
	const char* szname = tag.AttributeValue("name");

	// count nr of faces
	int faces = tag.children();

	// allocate storage for faces
	FEBModel::Surface* ps = new FEBModel::Surface(szname);
	part->AddSurface(ps);
	ps->Create(faces);

	// read faces
	++tag;
	for (int i = 0; i<faces; ++i)
	{
		FEBModel::FACET& face = ps->GetFacet(i);

		// get the ID (although we don't really use this)
		tag.AttributeValue("id", face.id);

		// set the facet type
		if (tag == "quad4") face.ntype = 4;
		else if (tag == "tri3") face.ntype = 3;
		else if (tag == "tri6") face.ntype = 6;
		else if (tag == "tri7") face.ntype = 7;
		else if (tag == "quad8") face.ntype = 8;
		else if (tag == "quad9") face.ntype = 9;
		else throw XMLReader::InvalidTag(tag);

		// we assume that the facet type also defines the number of nodes
		int N = face.ntype;
		tag.value(face.node, N);

		++tag;
	}
}

//-----------------------------------------------------------------------------
void FEBioGeometrySection25::ParseElementSetSection(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// get the name attribute
	const char* szname = tag.AttributeValue("name");

	// create a new element set
	FEElementSet* pg = new FEElementSet(&fem);
	pg->SetName(szname);

	vector<int> l;
	++tag;
	do
	{
		if (tag == "elem")
		{
			int nid = -1;
			tag.AttributeValue("id", nid);
			l.push_back(nid);
		}
		else { delete pg; throw XMLReader::InvalidTag(tag); }
		++tag;
	} while (!tag.isend());

	// only add non-empty element sets
	if (l.empty() == false)
	{
		// see if all elements belong to the same domain
		bool oneDomain = true;
		FEElement* el = mesh.FindElementFromID(l[0]); assert(el);
		FEDomain* dom = dynamic_cast<FEDomain*>(el->GetMeshPartition());
		for (int i = 1; i < l.size(); ++i)
		{
			FEElement* el_i = mesh.FindElementFromID(l[i]); assert(el);
			FEDomain* dom_i = dynamic_cast<FEDomain*>(el_i->GetMeshPartition());

			if (dom != dom_i)
			{
				oneDomain = false;
				break;
			}
		}

		// assign indices to element set
		if (oneDomain)
			pg->Create(dom, l);
		else
			pg->Create(l);

		// add the element set to the mesh
		mesh.AddElementSet(pg);
	}
	else delete pg;
}
