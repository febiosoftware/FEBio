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
#include "FEBioBoundarySection3.h"
#include <FECore/FEBoundaryCondition.h>
#include <FEBioMech/FEMechModel.h>
#include <FEBioMech/FERigidSystem.h>

//-----------------------------------------------------------------------------
void FEBioBoundarySection3::Parse(XMLTag& tag)
{
	if (tag.isleaf()) return;

	// build the node set map for faster lookup
	BuildNodeSetMap();

	FEModel* fem = GetFEModel();

	++tag;
	do
	{
		if (tag == "bc")
		{
			// get the type string
			const char* sztype = tag.AttributeValue("type");

			// create the boundary condition
			FEBoundaryCondition* pbc = fecore_new<FEBoundaryCondition>(sztype, fem);
			if (pbc == 0) throw XMLReader::InvalidTag(tag);

			// add this boundary condition to the current step
			GetBuilder()->AddBC(pbc);

			// Read the parameter list
			ReadParameterList(tag, pbc);
		}
		else if (tag == "rigid") ParseBCRigid(tag);
		else throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());
}

//-----------------------------------------------------------------------------
// Rigid node sets are defined in the Boundary section since version 2.5
// (Used to be defined in the Contact section)
void FEBioBoundarySection3::ParseBCRigid(XMLTag& tag)
{
	FEMechModel& fem = static_cast<FEMechModel&>(*GetFEModel());
	FEMesh& mesh = fem.GetMesh();
	FERigidSystem& rigid = *fem.GetRigidSystem();
	FEModelBuilder* feb = GetBuilder();
	int NMAT = fem.Materials();

	// get the rigid body material ID
	int rb = -1;
	tag.AttributeValue("rb", rb);
	rb -= 1;

	// make sure we have a valid rigid body reference
	if ((rb < 0) || (rb >= NMAT)) throw XMLReader::InvalidAttributeValue(tag, "rb", tag.AttributeValue("rb"));

	// get the nodeset
	const char* szset = tag.AttributeValue("node_set");
	FENodeSet* nodeSet = mesh.FindNodeSet(szset);
	if (nodeSet == 0) throw XMLReader::InvalidAttributeValue(tag, "node_set", szset);

	// create new rigid node set
	FERigidNodeSet* prn = new FERigidNodeSet(&fem);

	// the default shell bc depends on the shell formulation
	prn->SetShellBC(feb->m_default_shell == OLD_SHELL ? FERigidNodeSet::HINGED_SHELL : FERigidNodeSet::CLAMPED_SHELL);

	prn->SetRigidID(rb);
	prn->SetNodeSet(*nodeSet);

	rigid.AddRigidNodeSet(prn);

	// add it to the current step
	GetBuilder()->AddComponent(prn);

	// read the parameter list
	ReadParameterList(tag, prn);
}
