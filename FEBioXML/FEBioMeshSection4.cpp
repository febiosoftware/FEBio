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
#include "FEBioMeshSection4.h"
#include <FECore/FESolidDomain.h>
#include <FECore/FEShellDomain.h>
#include <FECore/FETrussDomain.h>
#include <FECore/FEDomain2D.h>
#include <FECore/FEModel.h>
#include <FECore/FEMaterial.h>
#include <FECore/FECoreKernel.h>
#include <FECore/FENodeNodeList.h>
#include <sstream>

//-----------------------------------------------------------------------------
FEBioMeshSection4::FEBioMeshSection4(FEBioImport* pim) : FEBioFileSection(pim) {}

//-----------------------------------------------------------------------------
void FEBioMeshSection4::Parse(XMLTag& tag)
{
	FEModelBuilder* builder = GetBuilder();
	builder->m_maxid = 0;

	// create a default part
	// NOTE: Do not specify a name for the part, otherwise
	//       all lists will be given the name: partname.listname
	FEBModel& feb = builder->GetFEBModel();
	assert(feb.Parts() == 0);
	FEBModel::Part* part = feb.AddPart("");

	// read all sections
	++tag;
	do
	{
		if      (tag == "Nodes"      ) ParseNodeSection       (tag, part);
		else if (tag == "Elements"   ) ParseElementSection    (tag, part);
		else if (tag == "NodeSet"    ) ParseNodeSetSection    (tag, part);
		else if (tag == "Surface"    ) ParseSurfaceSection    (tag, part);
		else if (tag == "Edge"       ) ParseEdgeSection       (tag, part);
		else if (tag == "ElementSet" ) ParseElementSetSection (tag, part);
		else if (tag == "SurfacePair") ParseSurfacePairSection(tag, part);
		else if (tag == "DiscreteSet") ParseDiscreteSetSection(tag, part);
		else throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());
}

//-----------------------------------------------------------------------------
//! Reads the Nodes section of the FEBio input file
void FEBioMeshSection4::ParseNodeSection(XMLTag& tag, FEBModel::Part* part)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	int N0 = mesh.Nodes();

	// see if this list defines a set
	const char* szname = tag.AttributeValue("name", true);
	FEBModel::NodeSet* ps = 0;
	if (szname)
	{
		ps = new FEBModel::NodeSet(szname);
		part->AddNodeSet(ps);
	}

	// allocate node

	vector<FEBModel::NODE> node; node.reserve(10000);
	vector<int> nodeList; nodeList.reserve(10000);

	// read nodal coordinates
	++tag;
	do {
		// nodal coordinates
		FEBModel::NODE nd;
		value(tag, nd.r);

		// get the nodal ID
		tag.AttributeValue("id", nd.id);

		// add it to the pile
		node.push_back(nd);
		nodeList.push_back(nd.id);

		// go on to the next node
		++tag;
	} while (!tag.isend());

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
void FEBioMeshSection4::ParseElementSection(XMLTag& tag, FEBModel::Part* part)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// get the (required!) name
	const char* szname = tag.AttributeValue("name");

	// get the element spec
	const char* sztype = tag.AttributeValue("type");
	FE_Element_Spec espec = GetBuilder()->ElementSpec(sztype);
	if (FEElementLibrary::IsValid(espec) == false) throw FEBioImport::InvalidElementType();

	// make sure the domain does not exist yet
	FEBModel::Domain* dom = part->FindDomain(szname);
	if (dom)
	{
		stringstream ss;
		ss << "Duplicate part name found : " << szname;
		throw std::runtime_error(ss.str());
	}

	// create the new domain
	dom = new FEBModel::Domain(espec);
	if (szname) dom->SetName(szname);

	// add domain it to the mesh
	part->AddDomain(dom);

	// for named domains, we'll also create an element set
	FEBModel::ElementSet* pg = 0;
	if (szname)
	{
		pg = new FEBModel::ElementSet(szname);
		part->AddElementSet(pg);
	}

	dom->Reserve(10000);
	vector<int> elemList; elemList.reserve(10000);

	// read element data
	++tag;
	do
	{
		FEBModel::ELEMENT el;

		// get the element ID
		tag.AttributeValue("id", el.id);

		// read the element data
		tag.value(el.node, FEElement::MAX_NODES);

		dom->AddElement(el);
		elemList.push_back(el.id);

		// go to next tag
		++tag;
	} while (!tag.isend());

	// set the element list
	if (pg) pg->SetElementList(elemList);
}

//-----------------------------------------------------------------------------
//! Reads the Geometry::Groups section of the FEBio input file
void FEBioMeshSection4::ParseNodeSetSection(XMLTag& tag, FEBModel::Part* part)
{
	// get the required name attribute
	const char* szname = tag.AttributeValue("name");

	// make sure it's a leaf
	if (tag.isleaf() == false) throw XMLReader::InvalidValue(tag);

	// see if a nodeset with this name already exists
	FEBModel::NodeSet* set = part->FindNodeSet(szname);
	if (set) throw FEBioImport::RepeatedNodeSet(szname);

	// read the node IDs
	vector<int> nodeList;
	tag.value(nodeList);

	// create the node set
	set = new FEBModel::NodeSet(szname);
	part->AddNodeSet(set);

	// add nodes to the list
	set->SetNodeList(nodeList);
}

//-----------------------------------------------------------------------------
void FEBioMeshSection4::ParseSurfaceSection(XMLTag& tag, FEBModel::Part* part)
{
	// get the mesh
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// get the required name attribute
	const char* szname = tag.AttributeValue("name");

	// see if this surface was already defined
	FEBModel::Surface* ps = part->FindSurface(szname);
	if (ps) throw FEBioImport::RepeatedSurface(szname);

	// count nr of faces
	int faces = tag.children();

	// allocate storage for faces
	ps = new FEBModel::Surface(szname);
	part->AddSurface(ps);
	ps->Create(faces);

	// read faces
	++tag;
	for (int i = 0; i < faces; ++i)
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
void FEBioMeshSection4::ParseElementSetSection(XMLTag& tag, FEBModel::Part* part)
{
	// get the mesh
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// get the required name attribute
	const char* szname = tag.AttributeValue("name");

	// see if this element set was already defined
	FEBModel::ElementSet* ps = part->FindElementSet(szname);
	if (ps) throw FEBioImport::RepeatedElementSet(szname);

	// allocate storage for faces
	ps = new FEBModel::ElementSet(szname);
	part->AddElementSet(ps);

	// read elements
	vector<int> elemList;
	tag.value(elemList);

	if (elemList.empty()) throw XMLReader::InvalidTag(tag);

	ps->SetElementList(elemList);
}

//-----------------------------------------------------------------------------
void FEBioMeshSection4::ParseEdgeSection(XMLTag& tag, FEBModel::Part* part)
{
	// get the mesh
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// get the required name attribute
	const char* szname = tag.AttributeValue("name");

	// see if this edge was already defined
	FEBModel::EdgeSet* ps = part->FindEdgeSet(szname);
	if (ps) throw FEBioImport::RepeatedEdgeSet(szname);

	// count nr of edges
	int edges = tag.children();

	// allocate storage for edges
	ps = new FEBModel::EdgeSet(szname);
	part->AddEdgeSet(ps);
	

	// read edges
	vector<FEBModel::EDGE> edgeSet(edges);
	++tag;
	for (int i = 0; i < edges; ++i)
	{
		FEBModel::EDGE& edge = edgeSet[i];

		// get the ID (although we don't really use this)
		tag.AttributeValue("id", edge.id);

		// set the facet type
		if      (tag == "line2") edge.ntype = 2;
		else if (tag == "line3") edge.ntype = 3;
		else throw XMLReader::InvalidTag(tag);

		// we assume that the facet type also defines the number of nodes
		int N = edge.ntype;
		tag.value(edge.node, N);

		++tag;
	}

	ps->SetEdgeList(edgeSet);
}

//-----------------------------------------------------------------------------
void FEBioMeshSection4::ParseSurfacePairSection(XMLTag& tag, FEBModel::Part* part)
{
	const char* szname = tag.AttributeValue("name");
	FEBModel::SurfacePair* surfPair = new FEBModel::SurfacePair;
	surfPair->m_name = szname;
	part->AddSurfacePair(surfPair);

	++tag;
	do
	{
		if (tag == "primary")
		{
			const char* sz = tag.szvalue();
			FEBModel::Surface* surf = part->FindSurface(sz);
			if (surf == nullptr) throw XMLReader::InvalidValue(tag);
			surfPair->m_primary = sz;
		}
		else if (tag == "secondary")
		{
			const char* sz = tag.szvalue();
			FEBModel::Surface* surf = part->FindSurface(sz);
			if (surf == nullptr) throw XMLReader::InvalidValue(tag);
			surfPair->m_secondary = sz;
		}
		else throw XMLReader::InvalidTag(tag);
		++tag;
	} 
	while (!tag.isend());
}

//-----------------------------------------------------------------------------
void FEBioMeshSection4::ParseDiscreteSetSection(XMLTag& tag, FEBModel::Part* part)
{
	const char* szname = tag.AttributeValue("name");
	FEBModel::DiscreteSet* dset = new FEBModel::DiscreteSet;
	dset->SetName(szname);
	part->AddDiscreteSet(dset);

	int n[2];
	++tag;
	do
	{
		if (tag == "delem")
		{
			tag.value(n, 2);
			dset->AddElement(n[0], n[1]);
		}
		else throw XMLReader::InvalidTag(tag);
		++tag;
	} while (!tag.isend());
}
