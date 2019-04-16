/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, Columbia University, and others.

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
	} while (!tag.isend());
}

//-----------------------------------------------------------------------------
void FEBioBoundarySection3::ParseBC(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// get the type string
	const char* sztype = tag.AttributeValue("type");

	// create the boundary condition
	FEBoundaryCondition* pbc = fecore_new<FEBoundaryCondition>(sztype, &fem);
	if (pbc == 0) throw XMLReader::InvalidTag(tag);

	// get the node set
	const char* szset = tag.AttributeValue("node_set", true);
	if (szset)
	{
		FENodeSet* nodeSet = mesh.FindNodeSet(szset);
		if (nodeSet == 0) throw XMLReader::InvalidAttributeValue(tag, "node_set", szset);

		// add the nodes to the BC
		pbc->AddNodes(*nodeSet);

		// Read the parameter list
		FEParameterList& pl = pbc->GetParameterList();
		ReadParameterList(tag, pl);
	}
	else
	{
		// if a node set is not defined, see if a surface is defined
		szset = tag.AttributeValue("surface");
		FEFacetSet* set = mesh.FindFacetSet(szset);

		// Read the parameter list (before setting the surface)
		FEParameterList& pl = pbc->GetParameterList();
		ReadParameterList(tag, pl);

		// add the surface nodes
		pbc->AddNodes(*set);
	}

	// add this boundary condition to the current step
	GetBuilder()->AddBC(pbc);
}
