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
#include <FEBioMech/FEPointBodyForce.h>
#include <FECore/FEModel.h>
#include <FECore/FECoreKernel.h>
#include <FECore/FENodalLoad.h>
#include <FECore/FESurfaceLoad.h>
#include <FECore/FEEdgeLoad.h>
#include <FECore/FEEdge.h>

//=============================================================================
// FEBioLoadsSection1x
//=============================================================================

//-----------------------------------------------------------------------------
void FEBioLoadsSection1x::Parse(XMLTag& tag)
{
	// make sure this tag has children
	if (tag.isleaf()) return;

	++tag;
	do
	{
		if      (tag == "force"      ) ParseNodalLoad(tag);
		else if (tag == "body_force" ) ParseBodyForce(tag);
		else if (tag == "heat_source") ParseBodyLoad (tag);
		else ParseSurfaceLoad(tag);
		++tag;
	}
	while (!tag.isend());
}


//-----------------------------------------------------------------------------
// NOTE: note that this section used to be in the Globals section (version 1.1)
void FEBioLoadsSection1x::ParseBodyForce(XMLTag &tag)
{
	FEModel& fem = *GetFEModel();

	const char* szt = tag.AttributeValue("type", true);
	if (szt == 0) szt = "const";

	if (strcmp(szt, "point") == 0)
	{
		FEPointBodyForce* pf = fecore_alloc(FEPointBodyForce, &fem);
		FEParameterList& pl = pf->GetParameterList();
		++tag;
		do
		{
			if (tag == "a")
			{
				const char* szlc = tag.AttributeValue("lc");
				//						pf->lc[0] = pf->lc[1] = pf->lc[2] = atoi(szlc);
				value(tag, pf->m_a);
			}
			else if (tag == "node")
			{
				value(tag, pf->m_inode);
				pf->m_inode -= 1;
			}
			else if (ReadParameter(tag, pl) == false) throw XMLReader::InvalidTag(tag);
			++tag;
		} while (!tag.isend());

		fem.AddBodyLoad(pf);
	}
	else
	{
		// create the body force
		FEBodyLoad* pf = fecore_new<FEBodyLoad>(szt, &fem);
		if (pf == 0) throw XMLReader::InvalidAttributeValue(tag, "type", szt);

		// add it to the model
		fem.AddBodyLoad(pf);

		// read the parameter list
		ReadParameterList(tag, pf);
	}
}

//-----------------------------------------------------------------------------
void FEBioLoadsSection1x::ParseBodyLoad(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	FEBodyLoad* pbl = fecore_new<FEBodyLoad>(tag.Name(), &fem);
	if (pbl == 0) throw XMLReader::InvalidTag(tag);
	ReadParameterList(tag, pbl);
	fem.AddBodyLoad(pbl);
}

//-----------------------------------------------------------------------------
void FEBioLoadsSection1x::ParseNodalLoad(XMLTag &tag)
{
	// No longer supported
	throw XMLReader::InvalidTag(tag);
/*
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	DOFS& dofs = fem.GetDOFS();

	// count how many nodal forces there are
	int ncnf = tag.children();

	// read the prescribed data
	++tag;
	for (int i = 0; i<ncnf; ++i)
	{
		int n = atoi(tag.AttributeValue("id")) - 1, bc;
		const char* sz = tag.AttributeValue("bc");

		bc = dofs.GetDOF(sz);
		if (bc == -1) throw XMLReader::InvalidAttributeValue(tag, "bc", sz);

		sz = tag.AttributeValue("lc");
		int lc = atoi(sz) - 1;

		double scale = 0.0;
		value(tag, scale);

		FENodalLoad* pfc = dynamic_cast<FENodalLoad*>(fecore_new<FEBoundaryCondition>("nodal_load", &fem));
		pfc->SetDOF(bc);
		pfc->SetLoad(scale);
		pfc->AddNode(n);

		if (lc >= 0)
		{
			FEParam* p = pfc->GetParameter("scale");
			if (p == nullptr) throw XMLReader::InvalidTag(tag);
			fem.AttachLoadController(p, lc);
		}

		// add it to the model
		GetBuilder()->AddNodalLoad(pfc);

		++tag;
	}
*/
}

//-----------------------------------------------------------------------------
void FEBioLoadsSection1x::ParseSurfaceLoad(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();

	// count how many pressure cards there are
	int npr = tag.children();

	// create a new surface
	FESurface* psurf = fecore_alloc(FESurface, &fem);
	psurf->Create(npr);
	fem.GetMesh().AddSurface(psurf);

	// create surface load
	FESurfaceLoad* ps = fecore_new<FESurfaceLoad>(tag.Name(), &fem);
	if (ps == 0) throw XMLReader::InvalidTag(tag);

	FEModelBuilder* feb = GetBuilder();

	// read the pressure data
	++tag;
	int nf[FEElement::MAX_NODES], N;
	for (int i = 0; i<npr; ++i)
	{
		FESurfaceElement& el = psurf->Element(i);

		if (tag == "quad4") el.SetType(FE_QUAD4G4);
		else if (tag == "tri3") el.SetType(feb->m_ntri3);
		else if (tag == "tri6") el.SetType(feb->m_ntri6);
		else if (tag == "tri7") el.SetType(feb->m_ntri7);
		else if (tag == "tri10") el.SetType(feb->m_ntri10);
		else if (tag == "quad8") el.SetType(FE_QUAD8G9);
		else if (tag == "quad9") el.SetType(FE_QUAD9G9);
		else throw XMLReader::InvalidTag(tag);

		N = el.Nodes();
		tag.value(nf, N);
		for (int j = 0; j<N; ++j) el.m_node[j] = nf[j] - 1;

		++tag;
	}

	ps->SetSurface(psurf);

	// add it to the model
	GetBuilder()->AddSurfaceLoad(ps);
}

//=============================================================================
// FEBioLoadsSection2
//=============================================================================

//-----------------------------------------------------------------------------
void FEBioLoadsSection2::Parse(XMLTag& tag)
{
	// make sure this tag has children
	if (tag.isleaf()) return;

	++tag;
	do
	{
		if      (tag == "nodal_load"  ) ParseNodalLoad(tag);
		else if (tag == "surface_load") ParseSurfaceLoad(tag);
		else if (tag == "edge_load"   ) ParseEdgeLoad(tag);
		else if (tag == "body_load"   ) ParseBodyLoad(tag);
		else throw XMLReader::InvalidTag(tag);
		++tag;
	} while (!tag.isend());
}

//-----------------------------------------------------------------------------
void FEBioLoadsSection2::ParseBodyLoad(XMLTag& tag)
{
	const char* sztype = tag.AttributeValue("type");
	FEModel& fem = *GetFEModel();
	FEBodyLoad* pbl = fecore_new<FEBodyLoad>(sztype, &fem);
	if (pbl == 0) throw XMLReader::InvalidAttributeValue(tag, "type", sztype);
	ReadParameterList(tag, pbl);
	fem.AddBodyLoad(pbl);
}


//-----------------------------------------------------------------------------
void FEBioLoadsSection2::ParseNodalLoad(XMLTag &tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	DOFS& dofs = fem.GetDOFS();

	// count how many nodal forces there are
	int ncnf = tag.children();

	// get the bc
	const char* sz = tag.AttributeValue("bc");
	int bc = dofs.GetDOF(sz);
	if (bc == -1) throw XMLReader::InvalidAttributeValue(tag, "bc", sz);

	// get the load curve
	sz = tag.AttributeValue("lc");
	int lc = atoi(sz) - 1;

	// see if there is a set defined
	const char* szset = tag.AttributeValue("set", true);
	if (szset)
	{
		// make sure this is a leaf tag
		if (tag.isleaf() == false) throw XMLReader::InvalidValue(tag);

		// find the node set
		FENodeSet* pns = mesh.FindNodeSet(szset);
		if (pns == 0) throw XMLReader::InvalidAttributeValue(tag, "set", szset);

		// see if the scale attribute is defined
		double scale = 1.0;
		tag.AttributeValue("scale", scale, true);

		// create new nodal force
		FENodalDOFLoad* pfc = fecore_alloc(FENodalDOFLoad, &fem);
		pfc->SetDOF(bc);
		pfc->SetLoad(scale);

		// add it to the model
		GetBuilder()->AddNodalLoad(pfc);
	}
	else
	{
		// NOTE: no longer supported!
		throw XMLReader::InvalidTag(tag);
/*
		// read the prescribed data
		++tag;
		for (int i = 0; i<ncnf; ++i)
		{
			// get the nodal ID
			int n = atoi(tag.AttributeValue("id")) - 1;

			// get the load scale factor
			double scale = 0.0;
			value(tag, scale);

			// create new nodal force
			FENodalDOFLoad* pfc = fecore_alloc(FENodalDOFLoad, &fem);
			pfc->SetDOF(bc);
			pfc->SetLoad(scale);
			pfc->AddNode(n);

			if (lc >= 0)
			{
				FEParam* p = pfc->GetParameter("scale");
				if (p == nullptr) throw XMLReader::InvalidTag(tag);
				fem.AttachLoadController(p, lc);
			}

			// add it to the model
			GetBuilder()->AddNodalLoad(pfc);

			++tag;
		}
*/
	}
}

//-----------------------------------------------------------------------------
void FEBioLoadsSection2::ParseSurfaceLoad(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();

	// create surface load
	const char* sztype = tag.AttributeValue("type");
	FESurfaceLoad* psl = fecore_new<FESurfaceLoad>(sztype, &fem);
	if (psl == 0) throw XMLReader::InvalidTag(tag);

	// read name attribute
	const char* szname = tag.AttributeValue("name", true);
	if (szname) psl->SetName(szname);

	// create a new surface
	FESurface* psurf = fecore_alloc(FESurface, &fem);
	fem.GetMesh().AddSurface(psurf);

	// we need to find the surface tag first
	ParseSurfaceLoadSurface(tag, psurf);
	psl->SetSurface(psurf);

	// read the parameters
	FEParameterList& pl = psl->GetParameterList();

	// read the pressure data
	++tag;
	do
	{
		if (ReadParameter(tag, pl) == false)
		{
			if (tag == "surface")
			{
				// skip it, we already processed it
				tag.m_preader->SkipTag(tag);
			}
			else throw XMLReader::InvalidTag(tag);
		}
		else ++tag;
	} 
	while (!tag.isend());

	// add it to the model
	GetBuilder()->AddSurfaceLoad(psl);
}

//-----------------------------------------------------------------------------
void FEBioLoadsSection2::ParseSurfaceLoadSurface(XMLTag& tag, FESurface* psurf)
{
	FEModelBuilder* feb = GetBuilder();
	XMLTag tag2(tag);
	++tag2;
	do
	{
		if (tag2 == "surface")
		{
			// see if the surface is referenced by a set of defined explicitly
			const char* szset = tag2.AttributeValue("set", true);
			if (szset)
			{
				// make sure this tag does not have any children
				if (!tag2.isleaf()) throw XMLReader::InvalidTag(tag2);

				// see if we can find the facet set
				FEMesh& m = GetFEModel()->GetMesh();
				FEFacetSet* ps = m.FindFacetSet(szset);

				// create a surface from the facet set
				if (ps)
				{
					if (feb->BuildSurface(*psurf, *ps) == false) throw XMLReader::InvalidTag(tag2);
				}
				else throw XMLReader::InvalidAttributeValue(tag2, "set", szset);
			}
			else
			{
				// count how many pressure cards there are
				int npr = tag2.children();
				psurf->Create(npr);

				// read surface
				++tag2;
				int nf[FEElement::MAX_NODES], N;
				for (int i = 0; i<npr; ++i)
				{
					FESurfaceElement& el = psurf->Element(i);

					if      (tag2 == "quad4") el.SetType(FE_QUAD4G4);
					else if (tag2 == "tri3" ) el.SetType(feb->m_ntri3);
					else if (tag2 == "tri6" ) el.SetType(feb->m_ntri6);
					else if (tag2 == "tri7" ) el.SetType(feb->m_ntri7);
					else if (tag2 == "tri10") el.SetType(feb->m_ntri10);
					else if (tag2 == "quad8") el.SetType(FE_QUAD8G9);
					else if (tag2 == "quad9") el.SetType(FE_QUAD9G9);
					else throw XMLReader::InvalidTag(tag2);

					N = el.Nodes();
					tag2.value(nf, N);
					for (int j = 0; j<N; ++j) el.m_node[j] = nf[j] - 1;

					++tag2;
				}
			}

			++tag2;
		}
		else tag2.m_preader->SkipTag(tag2);
	}
	while (!tag2.isend());
}

//-----------------------------------------------------------------------------
void FEBioLoadsSection2::ParseEdgeLoad(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();

	// create edge load
	const char* sztype = tag.AttributeValue("type");
	FEEdgeLoad* pel = fecore_new<FEEdgeLoad>(sztype, &fem);
	if (pel == 0) throw XMLReader::InvalidTag(tag);

	// create a new edge
	FEEdge* pedge = new FEEdge(&fem);
	fem.GetMesh().AddEdge(pedge);
	pel->SetEdge(pedge);

	// read the parameters
	FEParameterList& pl = pel->GetParameterList();

	// read the pressure data
	++tag;
	do
	{
		if (ReadParameter(tag, pl) == false)
		{
			if (tag == "edge")
			{
				// see if the surface is referenced by a set of defined explicitly
				const char* szset = tag.AttributeValue("set", true);
				if (szset)
				{
					// make sure this tag does not have any children
					if (!tag.isleaf()) throw XMLReader::InvalidTag(tag);

					// see if we can find the facet set
					FEMesh& m = GetFEModel()->GetMesh();
					FESegmentSet* ps = m.FindSegmentSet(szset);

					// create an edge from the segment set
					if (ps)
					{
						if (GetBuilder()->BuildEdge(*pedge, *ps) == false) throw XMLReader::InvalidTag(tag);
						pel->Create(pedge->Elements());
					}
					else throw XMLReader::InvalidAttributeValue(tag, "set", szset);
				}
				else
				{
					// count how many load cards there are
					int npr = tag.children();
					pedge->Create(npr);
					pel->Create(npr);

					++tag;
					int nf[FEElement::MAX_NODES], N;
					for (int i = 0; i<npr; ++i)
					{
						FELineElement& el = pedge->Element(i);

						if (tag == "line2") el.SetType(FE_LINE2G1);
						else throw XMLReader::InvalidTag(tag);

						N = el.Nodes();
						tag.value(nf, N);
						for (int j = 0; j<N; ++j) el.m_node[j] = nf[j] - 1;

						++tag;
					}
				}
			}
			else throw XMLReader::InvalidTag(tag);
		}
		++tag;
	} while (!tag.isend());

	// add edge load to model
	GetBuilder()->AddEdgeLoad(pel);
}

//=============================================================================
// FEBioLoadsSection25
//=============================================================================

//-----------------------------------------------------------------------------
void FEBioLoadsSection25::Parse(XMLTag& tag)
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
void FEBioLoadsSection25::ParseBodyLoad(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
    FEMesh& mesh = fem.GetMesh();
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
		FEElementSet* elset = mesh.FindElementSet(szpart);
		if (elset == 0) throw XMLReader::InvalidAttributeValue(tag, "elem_set", szpart);

		pbl->SetDomainList(elset);
	}

	ReadParameterList(tag, pbl);
	fem.AddBodyLoad(pbl);
}

//-----------------------------------------------------------------------------
void FEBioLoadsSection25::ParseNodalLoad(XMLTag &tag)
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
	FENodalDOFLoad* pfc = fecore_alloc(FENodalDOFLoad, &fem);
	pfc->SetDOF(bc);
	pfc->SetNodeSet(nodeSet);

	// add it to the model
	GetBuilder()->AddNodalLoad(pfc);

	// read parameters
	ReadParameterList(tag, pfc);
}

//-----------------------------------------------------------------------------
void FEBioLoadsSection25::ParseSurfaceLoad(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// create surface load
	const char* sztype = tag.AttributeValue("type");
	FESurfaceLoad* psl = fecore_new<FESurfaceLoad>(sztype, &fem);
	if (psl == 0)
	{
		// there are some obsolete loads that have been moved elsewhere,
		// so let's check for those first;
		ParseObsoleteLoad(tag);
	}
	else
	{

		// read name attribute
		const char* szname = tag.AttributeValue("name", true);
		if (szname) psl->SetName(szname);

		// get the surface
		const char* szset = tag.AttributeValue("surface");
		FEFacetSet* pface = mesh.FindFacetSet(szset);
		if (pface == 0) throw XMLReader::InvalidAttributeValue(tag, "surface", szset);

		FESurface* psurf = fecore_alloc(FESurface, &fem);
		GetBuilder()->BuildSurface(*psurf, *pface);

		mesh.AddSurface(psurf);
		psl->SetSurface(psurf);

		// read the parameters
		if (!tag.isleaf()) ReadParameterList(tag, psl);

		// add it to the model
		GetBuilder()->AddSurfaceLoad(psl);
	}
}

//-----------------------------------------------------------------------------
void FEBioLoadsSection25::ParseEdgeLoad(XMLTag& tag)
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

//-----------------------------------------------------------------------------
void FEBioLoadsSection25::ParseObsoleteLoad(XMLTag& tag)
{
	FEModel* fem = GetFEModel();
	FEMesh& mesh = fem->GetMesh();
	const char* sztype = tag.AttributeValue("type");

	// this "load" was moved to the boundary section
	if (strcmp(sztype, "fluid rotational velocity") == 0)
	{
		FEBoundaryCondition* pbc = fecore_new<FEBoundaryCondition>(sztype, fem);
		if (pbc == nullptr) throw XMLReader::InvalidAttributeValue(tag, "type", sztype);

		// get the node_set property
		FEProperty* prop = pbc->FindProperty("node_set");
		if (prop == nullptr) throw XMLReader::InvalidAttributeValue(tag, "type", sztype);

		// get the surface
		const char* szset = tag.AttributeValue("surface");
		FEFacetSet* pface = mesh.FindFacetSet(szset);
		if (pface == 0) throw XMLReader::InvalidAttributeValue(tag, "surface", szset);

		// extract the nodeset from the surface
		FENodeList nodeList = pface->GetNodeList();
		FENodeSet* nodeSet = fecore_alloc(FENodeSet, fem);
		nodeSet->SetName(pface->GetName());
		nodeSet->Add(nodeList);
		mesh.AddNodeSet(nodeSet);

		// assign the surface
		prop->SetProperty(nodeSet);

		// Read the parameter list
		ReadParameterList(tag, pbc);

		// Add the boundary condition
		GetBuilder()->AddBC(pbc);
	}
	else throw XMLReader::InvalidAttributeValue(tag, "type", sztype);
}
