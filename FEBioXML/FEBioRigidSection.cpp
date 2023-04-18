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
#include "FEBioRigidSection.h"
#include <FECore/FEModel.h>
#include <FECore/FECoreKernel.h>
#include <FECore/FEModelComponent.h>
#include <FECore/FEModelLoad.h>
#include <FECore/FENLConstraint.h>
#include <FECore/FEBoundaryCondition.h>

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
		FEBoundaryCondition* pBC = fecore_new_class<FEBoundaryCondition>("FERigidFixedBCOld", fem);
		feb.AddRigidComponent(pBC);
		ReadParameterList(tag, pBC);
	}
	else if (strcmp(sztype, "prescribe") == 0)
	{
		// create the rigid displacement constraint
		FEBoundaryCondition* pDC = fecore_new_class<FEBoundaryCondition>("FERigidPrescribedOld", fem);
		feb.AddRigidComponent(pDC);
		ReadParameterList(tag, pDC);
	}
	else if (strcmp(sztype, "force") == 0)
	{
		string name;
		const char* szname = tag.AttributeValue("name", true);
		if (szname) name = szname;

		// we need to decide whether we want to apply a force or a moment, since these
		// are now two separate classes. Unfortunately, this means we first need to 
		// read all parameters, before we can allocate the correct class. 
		int rb = -1;
		int ntype = 0;
		int bc = -1;
		double val = 0.;
		bool brel = false;
		int lc = -1;
		++tag;
		do
		{
			if      (tag == "rb"       ) tag.value(rb);
			else if (tag == "value")
			{
				tag.value(val);
				const char* szlc = tag.AttributeValue("lc", true);
				if (szlc) lc = atoi(szlc) - 1;
			}
			else if (tag == "load_type") tag.value(ntype);
			else if (tag == "relative" ) tag.value(brel);
			else if (tag == "dof"      )
			{
				const char* sz = tag.szvalue();
				if ((strcmp(sz, "Rx") == 0) || (strcmp(sz, "0") == 0)) bc = 0;
				if ((strcmp(sz, "Ry") == 0) || (strcmp(sz, "1") == 0)) bc = 1;
				if ((strcmp(sz, "Rz") == 0) || (strcmp(sz, "2") == 0)) bc = 2;
				if ((strcmp(sz, "Ru") == 0) || (strcmp(sz, "3") == 0)) bc = 3;
				if ((strcmp(sz, "Rv") == 0) || (strcmp(sz, "4") == 0)) bc = 4;
				if ((strcmp(sz, "Rw") == 0) || (strcmp(sz, "5") == 0)) bc = 5;

				if (bc < 0)
				{
					throw XMLReader::InvalidValue(tag);
				}
			}
			else throw XMLReader::InvalidTag(tag);

			++tag;
		} while (!tag.isend());

		FEModelLoad* pFC = nullptr;
		if (bc < 3)
		{
			pFC = fecore_new_class<FEModelLoad>("FERigidBodyForce", fem);
			pFC->SetParameter("load_type", ntype);
			pFC->SetParameter("rb", rb);
			pFC->SetParameter("dof", bc);
			pFC->SetParameter("value", val);
			pFC->SetParameter("relative", brel);
		}
		else
		{
			pFC = fecore_new_class<FEModelLoad>("FERigidBodyMoment", fem);
			pFC->SetParameter("rb", rb);
			pFC->SetParameter("dof", bc - 3);
			pFC->SetParameter("value", val);
			pFC->SetParameter("relative", brel);
		}

		if (lc >= 0)
		{
			FEParam* p = pFC->GetParameter("value");
			if (p == nullptr) throw XMLReader::InvalidTag(tag);
			GetFEModel()->AttachLoadController(p, lc);
		}

		if (name.empty() == false)
		{
			pFC->SetName(name);
		}

		feb.AddModelLoad(pFC);
	}
	else if (strcmp(sztype, "initial_rigid_velocity") == 0)
	{
		FEBoundaryCondition* pic = fecore_new_class<FEBoundaryCondition>("FERigidBodyVelocity", fem);
		feb.AddRigidComponent(pic);
		ReadParameterList(tag, pic);
	}
	else if (strcmp(sztype, "initial_rigid_angular_velocity") == 0)
	{
		FEBoundaryCondition* pic = fecore_new_class<FEBoundaryCondition>("FERigidBodyAngularVelocity", fem);
		feb.AddRigidComponent(pic);
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
		FEBoundaryCondition* pBC = fecore_new<FEBoundaryCondition>(sztype, fem);
		if (pBC == nullptr)  throw XMLReader::InvalidAttributeValue(tag, "type", sztype);
		feb.AddRigidComponent(pBC);
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
