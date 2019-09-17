/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, The Trustees of Columbia University in 
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
#include "FEBioLoadsSection.h"
#include "FEBioMech/FEPointBodyForce.h"
#include "FECore/FEModel.h"
#include "FECore/FECoreKernel.h"
#include <FECore/FENodalLoad.h>
#include <FECore/FESurfaceLoad.h>
#include <FECore/FEEdgeLoad.h>
#include <FECore/FEEdge.h>

//-----------------------------------------------------------------------------
void FEBioLoadsSection3::Parse(XMLTag& tag)
{
	// make sure this tag has children
	if (tag.isleaf()) return;

	++tag;
	do
	{
		if      (tag == "nodal_load"  ) ParseNodalLoad  (tag);
		else if (tag == "surface_load") ParseSurfaceLoad(tag);
		else if (tag == "edge_load"   ) ParseEdgeLoad   (tag);
		else if (tag == "body_load"   ) ParseBodyLoad   (tag);
		else throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());
}

//-----------------------------------------------------------------------------
void FEBioLoadsSection3::ParseBodyLoad(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	const char* sztype = tag.AttributeValue("type");
	FEBodyLoad* pbl = fecore_new<FEBodyLoad>(sztype, &fem);
	if (pbl == 0) throw XMLReader::InvalidAttributeValue(tag, "type", sztype);

	// see if a name was defined
	const char* szname = tag.AttributeValue("name", true);
	if (szname) pbl->SetName(szname);

	// see if a specific domain was referenced
	const char* szpart = tag.AttributeValue("elem_set", true);
	if (szpart)
	{
		FEMesh& mesh = fem.GetMesh();
		FEElementSet* elset = mesh.FindElementSet(szpart);
		if (elset == 0) throw XMLReader::InvalidAttributeValue(tag, "elem_set", szpart);
		pbl->SetDomainList(elset);
	}

	ReadParameterList(tag, pbl);
	fem.AddBodyLoad(pbl);
}

//-----------------------------------------------------------------------------
void FEBioLoadsSection3::ParseNodalLoad(XMLTag &tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// get the type
	const char* sztype = tag.AttributeValue("type");

	// create nodal load
	FENodalLoad* pfc = fecore_new<FENodalLoad>(sztype, GetFEModel());

	// read (optional) name attribute
	const char* szname = tag.AttributeValue("name", true);
	if (szname) pfc->SetName(szname);

	// get the node set
	const char* szset = tag.AttributeValue("node_set");
	FENodeSet* nodeSet = GetBuilder()->FindNodeSet(szset);
	if (nodeSet == 0) throw XMLReader::InvalidAttributeValue(tag, "node_set", szset);

	// assign the node set
	pfc->SetNodeSet(nodeSet);

	// add it to the model
	GetBuilder()->AddNodalLoad(pfc);

	// read parameters
	ReadParameterList(tag, pfc);
}

//-----------------------------------------------------------------------------
void FEBioLoadsSection3::ParseSurfaceLoad(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// create surface load
	const char* sztype = tag.AttributeValue("type");
	FESurfaceLoad* psl = fecore_new<FESurfaceLoad>(sztype, &fem);
	if (psl == 0) throw XMLReader::InvalidTag(tag);

	// read (optional) name attribute
	const char* szname = tag.AttributeValue("name", true);
	if (szname) psl->SetName(szname);

	// read required surface attribute
	const char* surfaceName = tag.AttributeValue("surface");
	FEFacetSet* pface = mesh.FindFacetSet(surfaceName);
	if (pface == 0) throw XMLReader::InvalidAttributeValue(tag, "surface", surfaceName);

	// create a surface from this facet set
	FESurface* psurf = new FESurface(&fem);
	GetBuilder()->BuildSurface(*psurf, *pface);

	// assign it
	mesh.AddSurface(psurf);
	psl->SetSurface(psurf);

	// add it to the model
	GetBuilder()->AddSurfaceLoad(psl);

	// read the parameters
	ReadParameterList(tag, psl);
}

//-----------------------------------------------------------------------------
void FEBioLoadsSection3::ParseEdgeLoad(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// create edge load
	const char* sztype = tag.AttributeValue("type");
	FEEdgeLoad* pel = fecore_new<FEEdgeLoad>(sztype, &fem);
	if (pel == 0) throw XMLReader::InvalidTag(tag);

	// create a new edge
	FEEdge* pedge = new FEEdge(&fem);
	mesh.AddEdge(pedge);
	pel->SetEdge(pedge);

	// get the segment set
	const char* szedge = tag.AttributeValue("edge");
	FESegmentSet* pset = mesh.FindSegmentSet(szedge);
	if (pset == 0) throw XMLReader::InvalidAttributeValue(tag, "edge", szedge);
	if (GetBuilder()->BuildEdge(*pedge, *pset) == false) throw XMLReader::InvalidTag(tag);

	// read the parameters
	ReadParameterList(tag, pel);

	// add edge load to model
	GetBuilder()->AddEdgeLoad(pel);
}
