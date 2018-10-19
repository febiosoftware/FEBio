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
	FEBodyLoad* pbl = fecore_new<FEBodyLoad>(FEBODYLOAD_ID, sztype, &fem);
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
	DOFS& dofs = fem.GetDOFS();

	// get the bc
	const char* sz = tag.AttributeValue("bc");
	int bc = dofs.GetDOF(sz);
	if (bc == -1) throw XMLReader::InvalidAttributeValue(tag, "bc", sz);

	// get the node set
	const char* szset = tag.AttributeValue("node_set");
	FENodeSet* nodeSet = mesh.FindNodeSet(szset);
	if (nodeSet == 0) throw XMLReader::InvalidAttributeValue(tag, "node_set", szset);

	// create nodal load
	FENodalLoad* pfc = dynamic_cast<FENodalLoad*>(fecore_new<FEBoundaryCondition>(FEBC_ID, "nodal load", &fem));
	pfc->SetDOF(bc);
	pfc->AddNodes(*nodeSet);

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
	FESurfaceLoad* psl = fecore_new<FESurfaceLoad>(FESURFACELOAD_ID, sztype, &fem);
	if (psl == 0) throw XMLReader::InvalidTag(tag);

	// read name attribute
	const char* szname = tag.AttributeValue("name", true);
	if (szname) psl->SetName(szname);

	// get the surface
	const char* szset = tag.AttributeValue("surface");
	FEFacetSet* pface = mesh.FindFacetSet(szset);
	if (pface == 0) throw XMLReader::InvalidAttributeValue(tag, "surface", szset);

	FESurface* psurf = new FESurface(&fem, pface);
	GetBuilder()->BuildSurface(*psurf, *pface);
	
	mesh.AddSurface(psurf);
	psl->SetSurface(psurf);

	// read the parameters
	ReadParameterList(tag, psl);

	// add it to the model
	GetBuilder()->AddSurfaceLoad(psl);
}

//-----------------------------------------------------------------------------
void FEBioLoadsSection3::ParseEdgeLoad(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// create edge load
	const char* sztype = tag.AttributeValue("type");
	FEEdgeLoad* pel = fecore_new<FEEdgeLoad>(FEEDGELOAD_ID, sztype, &fem);
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
