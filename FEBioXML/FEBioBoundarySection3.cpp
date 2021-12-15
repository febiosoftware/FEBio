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
#include "FEBioBoundarySection3.h"
#include <FECore/FEBoundaryCondition.h>
#include <FECore/FEModel.h>
#include <FECore/FELinearConstraintManager.h>

//-----------------------------------------------------------------------------
void FEBioBoundarySection3::Parse(XMLTag& tag)
{
	if (tag.isleaf()) return;

	// build the node set map for faster lookup
	BuildNodeSetMap();

	++tag;
	do
	{
		if (tag == "bc") ParseBC(tag);
		else throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());
}

//-----------------------------------------------------------------------------
void FEBioBoundarySection3::ParseBC(XMLTag& tag)
{
	FEModel* fem = GetFEModel();
	FEMesh& mesh = fem->GetMesh();

	// get the type string
	const char* sztype = tag.AttributeValue("type");

	// handle "rigid" bc separately
	if (strcmp(sztype, "rigid") == 0)
	{
		ParseBCRigid(tag);
		return;
	}

	// handle linar constraints
	if (strcmp(sztype, "linear constraint") == 0)
	{
		ParseLinearConstraint(tag);
		return;
	}

	// create the boundary condition
	FEBoundaryCondition* pbc = fecore_new<FEBoundaryCondition>(sztype, fem);
	if (pbc == 0) throw XMLReader::InvalidTag(tag);

	// read the (optional) name 
	const char* szname = tag.AttributeValue("name", true);
	if (szname) pbc->SetName(szname);

	// get the node set
	FEProperty* pn = pbc->FindProperty("node_set");
	if (pn)
	{
		// read required node_set attribute
		const char* szset = tag.AttributeValue("node_set");
		FENodeSet* nodeSet = GetBuilder()->FindNodeSet(szset);
		if (nodeSet == nullptr) throw XMLReader::InvalidAttributeValue(tag, "node_set", szset);
		pn->SetProperty(nodeSet);
	}

	// get the surface
	FEProperty* ps = pbc->FindProperty("surface");
	if (ps)
	{
		// read required surface attribute
		const char* surfaceName = tag.AttributeValue("surface");
		FEFacetSet* pface = mesh.FindFacetSet(surfaceName);
		if (pface == 0) throw XMLReader::InvalidAttributeValue(tag, "surface", surfaceName);

		// create a surface from this facet set
		FESurface* psurf = fecore_alloc(FESurface, fem);
		GetBuilder()->BuildSurface(*psurf, *pface);

		// assign it
		mesh.AddSurface(psurf);
		ps->SetProperty(psurf);
	}

	// add this boundary condition to the current step
	GetBuilder()->AddBC(pbc);

	// Read the parameter list
	if (tag.isleaf() == false)
		ReadParameterList(tag, pbc);
}

//-----------------------------------------------------------------------------
// Rigid node sets are defined in the Boundary section since version 2.5
// (Used to be defined in the Contact section)
void FEBioBoundarySection3::ParseBCRigid(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	FEModelBuilder* feb = GetBuilder();
	int NMAT = fem.Materials();

	// get the nodeset
	const char* szset = tag.AttributeValue("node_set");
	FENodeSet* nodeSet = GetBuilder()->FindNodeSet(szset);
	if (nodeSet == 0) throw XMLReader::InvalidAttributeValue(tag, "node_set", szset);

	// create new rigid node set
	FEBoundaryCondition* prn = fecore_new_class<FEBoundaryCondition>("FERigidNodeSet", &fem);

	prn->SetNodeSet(nodeSet);

	// the default shell bc depends on the shell formulation
	// hinged shell = 0
	// clamped shell = 1
	prn->SetParameter("clamp_shells", feb->m_default_shell == OLD_SHELL ? 0 : 1);

	feb->AddBC(prn);

	// read the parameter list
	ReadParameterList(tag, prn);
}


//-----------------------------------------------------------------------------
//! Parse the linear constraints section of the xml input file
//! This section is a subsection of the Boundary section

void FEBioBoundarySection3::ParseLinearConstraint(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	DOFS& dofs = fem.GetDOFS();

	// make sure there is a constraint defined
	if (tag.isleaf()) return;

	FEModelBuilder* feb = GetBuilder();

	FELinearConstraint* lc = new FELinearConstraint(&fem);

	++tag;
	do
	{
		if (tag == "node")
		{
			int nodeId;
			tag.value(nodeId);
			lc->SetParentNode(feb->FindNodeFromID(nodeId));
		}
		else if (tag == "dof")
		{
			lc->SetParentDof(dofs.GetDOF(tag.szvalue()));
		}
		else if (tag == "child_dof")
		{
			FELinearConstraintDOF* dof = new FELinearConstraintDOF(&fem);
			++tag;
			do
			{
				if (tag == "node")
				{
					int nodeId;
					tag.value(nodeId);
					dof->node = feb->FindNodeFromID(nodeId);
				}
				else if (tag == "dof")
				{
					dof->dof = dofs.GetDOF(tag.szvalue());
				}
				else if (tag == "value")
				{
					double v;
					tag.value(v);
					dof->val = v;

					const char* szlc = tag.AttributeValue("lc", true);
					if (szlc)
					{
						int lc = atoi(szlc) - 1;
						FEParam* pp = dof->FindParameter("value"); assert(pp);
						if (pp) fem.AttachLoadController(pp, lc);
					}
				}
				else throw XMLReader::InvalidTag(tag);
				++tag;
			} 
			while (!tag.isend());

			lc->AddChildDof(dof);
		}
		else throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());

	// add the linear constraint to the system
	fem.GetLinearConstraintManager().AddLinearConstraint(lc);
	GetBuilder()->AddComponent(lc);
}