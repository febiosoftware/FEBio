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
#include "FEBioRigidSection.h"
#include <FECore/FEModel.h>
#include <FECore/FECoreKernel.h>
#include <FEBioMech/FERigidSurface.h>
#include <FEBioMech/FEMechModel.h>
#include <FEBioMech/FERigidMaterial.h>
#include <FEBioMech/FERigidForce.h>

void FEBioRigidSection::Parse(XMLTag& tag)
{
	if (tag.isleaf()) return;

	FEMechModel& fem = static_cast<FEMechModel&>(*GetFEModel());

	++tag;
	do
	{
		if (tag == "rigid_constraint") ParseRigidBC(tag);
		else if (tag == "initial_velocity")
		{
			const char* szm = tag.AttributeValue("mat");
			assert(szm);

			// get the material ID
			int nmat = atoi(szm);
			if ((nmat <= 0) || (nmat > fem.Materials())) throw XMLReader::InvalidAttributeValue(tag, "mat", szm);

			// get the initial velocity
			vec3d v;
			value(tag, v);

			// create the initial condition
			FERigidBodyVelocity* pic = new FERigidBodyVelocity(&fem);
			pic->m_rid = nmat;
			pic->m_vel = v;
			fem.AddRigidInitialCondition(pic);

			// add this initial condition to the current step
			GetBuilder()->AddComponent(pic);
		}
		else if (tag == "initial_angular_velocity")
		{
			const char* szm = tag.AttributeValue("mat");
			assert(szm);

			// get the material ID
			int nmat = atoi(szm);
			if ((nmat <= 0) || (nmat > fem.Materials())) throw XMLReader::InvalidAttributeValue(tag, "mat", szm);

			// get the initial angular velocity
			vec3d w;
			value(tag, w);

			// create the initial condition
			FERigidBodyAngularVelocity* pic = new FERigidBodyAngularVelocity(&fem);
			pic->m_rid = nmat;
			pic->m_w = w;
			fem.AddRigidInitialCondition(pic);

			// add this initial condition to the current step
			GetBuilder()->AddComponent(pic);
		}
		else throw XMLReader::InvalidTag(tag);
		++tag;
	} 
	while (!tag.isend());
}

bool parseRigidDofs(const char* sz, vector<int>& bc)
{
	bc.clear();
	int l = strlen(sz);
	if (l < 1) return false;

	char* buf = new char[l + 1];
	strcpy(buf, sz);
	buf[l] = 0;

	char *ch = buf;
	bool ok = true;
	do
	{
		char* ch2 = strchr(ch, ',');
		if (ch2) *ch2++ = 0;

		if      (strcmp(ch, "Rx") == 0) bc.push_back(0);
		else if (strcmp(ch, "Ry") == 0) bc.push_back(1);
		else if (strcmp(ch, "Rz") == 0) bc.push_back(2);
		else if (strcmp(ch, "Ru") == 0) bc.push_back(3);
		else if (strcmp(ch, "Rv") == 0) bc.push_back(4);
		else if (strcmp(ch, "Rw") == 0) bc.push_back(5);
		else
		{
			ok = false;
			break;
		}

		ch = ch2;
	}
	while (ch);
	delete buf;
	return ok;
}

void FEBioRigidSection::ParseRigidBC(XMLTag& tag)
{
	FEMechModel& fem = static_cast<FEMechModel&>(*GetFEModel());

	// get the type
	const char* sztype = tag.AttributeValue("type");

	if (strcmp(sztype, "fix") == 0)
	{
		// create the fixed dof
		FERigidBodyFixedBC* pBC = static_cast<FERigidBodyFixedBC*>(fecore_new<FERigidBC>("rigid_fixed", &fem));
		fem.AddRigidFixedBC(pBC);

		// add this boundary condition to the current step
		GetBuilder()->AddComponent(pBC);

		// read parameters
		ReadParameterList(tag, pBC);
	}
	else if (strcmp(sztype, "prescribe") == 0)
	{
		// create the rigid displacement constraint
		FERigidBodyDisplacement* pDC = static_cast<FERigidBodyDisplacement*>(fecore_new<FERigidBC>("rigid_prescribed", &fem));
		fem.AddRigidPrescribedBC(pDC);

		GetBuilder()->AddComponent(pDC);

		// read parameters
		ReadParameterList(tag, pDC);
	}
	else if (strcmp(sztype, "force") == 0)
	{
		// create the rigid body force
		FERigidBodyForce* pFC = static_cast<FERigidBodyForce*>(fecore_new<FEModelLoad>(FEBC_ID, "rigid_force", &fem));

		// add it to the model
		GetBuilder()->AddModelLoad(pFC);

		// read the parameterlist
		ReadParameterList(tag, pFC);
	}
	else if (tag == "rigid_velocity")
	{
		FERigidBodyVelocity* rc = fecore_alloc(FERigidBodyVelocity, &fem);
		GetBuilder()->AddRigidIC(rc);

		ReadParameterList(tag, rc);
	}
	else if (tag == "rigid_angular_velocity")
	{
		FERigidBodyAngularVelocity* rc = fecore_alloc(FERigidBodyAngularVelocity, &fem);
		GetBuilder()->AddRigidIC(rc);

		ReadParameterList(tag, rc);
	}
	else throw XMLReader::InvalidAttributeValue(tag, "type", sztype);
}
