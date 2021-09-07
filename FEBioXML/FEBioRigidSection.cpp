/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2020 University of Utah, The Trustees of Columbia University in 
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
#include "FEBioRigidSection.h"
#include <FECore/FEModel.h>
#include <FECore/FECoreKernel.h>
#include <FECore/FEModelComponent.h>
#include <FECore/FEModelLoad.h>

void FEBioRigidSection::Parse(XMLTag& tag)
{
	if (tag.isleaf()) return;

	++tag;
	do
	{
		if      (tag == "rigid_constraint") ParseRigidBC(tag);
		else if (tag == "rigid_connector" ) ParseRigidConnector(tag);
		else throw XMLReader::InvalidTag(tag);
		++tag;
	} 
	while (!tag.isend());
}

void FEBioRigidSection::ParseRigidBC(XMLTag& tag)
{
	FEModel* fem = GetFEModel();
	FEModelBuilder& feb = *GetBuilder();

	// get the type
	const char* sztype = tag.AttributeValue("type");

	if (strcmp(sztype, "fix") == 0)
	{
		// create the fixed dof
		FEModelComponent* pBC = fecore_new_class<FEModelComponent>("FERigidBodyFixedBC", fem);
		feb.AddRigidBC(pBC);
		ReadParameterList(tag, pBC);
	}
	else if (strcmp(sztype, "prescribe") == 0)
	{
		// create the rigid displacement constraint
		FEModelComponent* pDC = fecore_new_class<FEModelComponent>("FERigidBodyDisplacement", fem);
		feb.AddRigidBC(pDC);
		ReadParameterList(tag, pDC);
	}
	else if (strcmp(sztype, "force") == 0)
	{
		// create the rigid body force
		FEModelLoad* pFC = fecore_new<FEModelLoad>(FERIGIDLOAD_ID, "rigid_force", fem);
		feb.AddRigidBC(pFC);
		ReadParameterList(tag, pFC);
	}
	else if (strcmp(sztype, "initial_rigid_velocity") == 0)
	{
		FEModelComponent* pic = fecore_new_class<FEModelComponent>("FERigidBodyVelocity", fem);
		feb.AddRigidBC(pic);
		ReadParameterList(tag, pic);
	}
	else if (strcmp(sztype, "initial_rigid_angular_velocity") == 0)
	{
		FEModelComponent* pic = fecore_new_class<FEModelComponent>("FERigidBodyAngularVelocity", fem);
		feb.AddRigidBC(pic);
		ReadParameterList(tag, pic);
	}
    else if (strcmp(sztype, "follower force") == 0)
    {
        FEModelLoad* rc = fecore_new_class<FEModelLoad>("FERigidFollowerForce", fem);
        feb.AddModelLoad(rc);
        ReadParameterList(tag, rc);
    }
    else if (strcmp(sztype, "follower moment") == 0)
    {
		FEModelLoad* rc = fecore_new_class<FEModelLoad>("FERigidFollowerMoment", fem);
        feb.AddModelLoad(rc);
        ReadParameterList(tag, rc);
    }
	else
	{
		// create the rigid constraint
		FEModelComponent* pBC = fecore_new<FEModelComponent>(FERIGIDBC_ID, sztype, fem);
		if (pBC == nullptr)  throw XMLReader::InvalidAttributeValue(tag, "type", sztype);
		feb.AddRigidBC(pBC);
		ReadParameterList(tag, pBC);
	}
}

void FEBioRigidSection::ParseRigidConnector(XMLTag& tag)
{
	const char* sztype = tag.AttributeValue("type");

	FENLConstraint* plc = fecore_new<FENLConstraint>(sztype, GetFEModel());
	if (plc == 0) throw XMLReader::InvalidAttributeValue(tag, "type", sztype);

	const char* szname = tag.AttributeValue("name", true);
	if (szname) plc->SetName(szname);

	// read the parameter list
	ReadParameterList(tag, plc);

	// add this constraint to the current step
	GetBuilder()->AddNonlinearConstraint(plc);
}
