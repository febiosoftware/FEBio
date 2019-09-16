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
#include <FEBioMech/FERigidSystem.h>
#include <FECore/FECoreKernel.h>
#include <FEBioMech/FERigidSurface.h>
#include <FEBioMech/FEMechModel.h>
#include <FEBioMech/FERigidMaterial.h>
#include <FEBioMech/FERigidForce.h>

void FEBioRigidSection::Parse(XMLTag& tag)
{
	if (tag.isleaf()) return;

	FEMechModel& fem = static_cast<FEMechModel&>(*GetFEModel());
	FERigidSystem& rigid = *fem.GetRigidSystem();

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
			rigid.AddInitialVelocity(pic);

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
			rigid.AddInitialAngularVelocity(pic);

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

		if      (strcmp(ch, "x" ) == 0) bc.push_back(0);
		else if (strcmp(ch, "y" ) == 0) bc.push_back(1);
		else if (strcmp(ch, "z" ) == 0) bc.push_back(2);
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
	FERigidSystem& rigid = *fem.GetRigidSystem();

	const char* szm = tag.AttributeValue("mat");
	assert(szm);

	// get the material ID
	int nmat = atoi(szm);
	if ((nmat <= 0) || (nmat > fem.Materials())) throw XMLReader::InvalidAttributeValue(tag, "mat", szm);

	// make sure this is a valid rigid material
	FERigidMaterial* pm = dynamic_cast<FERigidMaterial*>(fem.GetMaterial(nmat - 1));
	if (pm == 0) throw XMLReader::InvalidAttributeValue(tag, "mat", szm);

	const char* sztype = tag.AttributeValue("type");

	if (strcmp(sztype, "fix") == 0)
	{
		vector<int> dofList;
		// get the dof list
		++tag;
		do
		{
			if (tag == "dofs")
			{
				const char* sz = tag.szvalue();
				if (parseRigidDofs(sz, dofList) == false) throw XMLReader::InvalidValue(tag);
			}
			else throw XMLReader::InvalidTag(tag);
			++tag;
		}
		while (!tag.isend());

		// create a separate rigid BC for each dof
		for (int i = 0; i < dofList.size(); ++i)
		{
			int bc = dofList[i];

			// create the fixed dof
			FERigidBodyFixedBC* pBC = static_cast<FERigidBodyFixedBC*>(fecore_new<FERigidBC>("rigid_fixed", &fem));
			pBC->id = nmat;
			pBC->bc = bc;
			rigid.AddFixedBC(pBC);

			// add this boundary condition to the current step
			GetBuilder()->AddComponent(pBC);
		}
	}
	else if (strcmp(sztype, "prescribe") == 0)
	{
		// create the rigid displacement constraint
		FERigidBodyDisplacement* pDC = static_cast<FERigidBodyDisplacement*>(fecore_new<FERigidBC>("rigid_prescribed", &fem));
		rigid.AddPrescribedBC(pDC);
		pDC->SetID(nmat);

		GetBuilder()->AddComponent(pDC);

		++tag;
		do
		{
			if (tag == "dof")
			{
				int bc = -1;
				const char* szbc = tag.szvalue();
				if      (strcmp(szbc, "x") == 0) bc = 0;
				else if (strcmp(szbc, "y") == 0) bc = 1;
				else if (strcmp(szbc, "z") == 0) bc = 2;
				else if (strcmp(szbc, "Ru") == 0) bc = 3;
				else if (strcmp(szbc, "Rv") == 0) bc = 4;
				else if (strcmp(szbc, "Rw") == 0) bc = 5;
				else throw XMLReader::InvalidValue(tag);
				pDC->SetBC(bc);
			}
			else if (tag == "value")
			{
				// get the loadcurve
				const char* szlc = tag.AttributeValue("lc");
				int lc = atoi(szlc) - 1;

				double val = 0.0;
				value(tag, val);
				pDC->SetValue(val);

				// assign a load curve
				if (lc >= 0)
				{
					FEParam* p = pDC->GetParameter("value");
					if (p == nullptr) throw XMLReader::InvalidTag(tag);
					GetFEModel()->AttachLoadController(p, lc);
				}
			}
			else if (tag == "relative")
			{
				bool b = true;
				tag.value(b);
				pDC->SetRelativeFlag(b);

			}
			++tag;
		}
		while (!tag.isend());
	}
	else if (strcmp(sztype, "force") == 0)
	{
		// create the rigid body force
		FERigidBodyForce* pFC = static_cast<FERigidBodyForce*>(fecore_new<FEModelLoad>(FEBC_ID, "rigid_force", &fem));
		pFC->SetID(nmat);

		// add it to the model
		GetBuilder()->AddModelLoad(pFC);

		++tag;
		do
		{
			if (tag == "dof")
			{
				const char* szbc = tag.szvalue();
				// get the dof
				int bc = -1;
				if      (strcmp(szbc, "x") == 0) bc = 0;
				else if (strcmp(szbc, "y") == 0) bc = 1;
				else if (strcmp(szbc, "z") == 0) bc = 2;
				else if (strcmp(szbc, "Ru") == 0) bc = 3;
				else if (strcmp(szbc, "Rv") == 0) bc = 4;
				else if (strcmp(szbc, "Rw") == 0) bc = 5;
				else throw XMLReader::InvalidAttributeValue(tag, "bc", szbc);
				pFC->SetBC(bc);
			}
			else if (tag == "value")
			{ 
				// get the loadcurve
				const char* szlc = tag.AttributeValue("lc", true);
				int lc = -1;
				if (szlc) lc = atoi(szlc) - 1;
				if (lc >= 0)
				{
					FEParam* p = pFC->GetParameter("force");
					if (p == nullptr) throw XMLReader::InvalidTag(tag);
					GetFEModel()->AttachLoadController(p, lc);
				}

				double val = 0.0;
				value(tag, val);
				pFC->SetForce(val);
			}
			else if (tag == "follow")
			{
				bool b = true;
				value(tag, b);
				pFC->SetFollowFlag(b);
			}
			// TODO: handle the former type attribute, which decides if the force is a ramp or not
			else throw XMLReader::InvalidTag(tag);

			++tag;
		}
		while (!tag.isend());

		// make sure there is a loadcurve for type=0 forces
//		if ((ntype == 0) && (lc == -1)) throw XMLReader::MissingAttribute(tag, "lc");
	}
	else if (tag == "rigid_velocity")
	{
		FERigidBodyVelocity* rc = fecore_alloc(FERigidBodyVelocity, &fem);
		rc->m_rid = nmat;
		GetBuilder()->AddRigidBodyVelocity(rc);

		ReadParameterList(tag, rc->GetParameterList());
	}
	else if (tag == "rigid_angular_velocity")
	{
		FERigidBodyAngularVelocity* rc = fecore_alloc(FERigidBodyAngularVelocity, &fem);
		rc->m_rid = nmat;
		GetBuilder()->AddRigidBodyAngularVelocity(rc);

		ReadParameterList(tag, rc->GetParameterList());
	}
	else throw XMLReader::InvalidAttributeValue(tag, "type", sztype);
}
