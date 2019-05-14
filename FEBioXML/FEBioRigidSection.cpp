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
	++tag;
	do
	{
		if (tag == "rigid_body") ParseRigidBody(tag);
		else throw XMLReader::InvalidTag(tag);
		++tag;
	} 
	while (!tag.isend());
}

void FEBioRigidSection::ParseRigidBody(XMLTag& tag)
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

	++tag;
	do
	{
		if (tag == "prescribed")
		{
			// get the dof
			int bc = -1;
			const char* szbc = tag.AttributeValue("bc");
			if (strcmp(szbc, "x") == 0) bc = 0;
			else if (strcmp(szbc, "y") == 0) bc = 1;
			else if (strcmp(szbc, "z") == 0) bc = 2;
			else if (strcmp(szbc, "Rx") == 0) bc = 3;
			else if (strcmp(szbc, "Ry") == 0) bc = 4;
			else if (strcmp(szbc, "Rz") == 0) bc = 5;
			else throw XMLReader::InvalidAttributeValue(tag, "bc", szbc);

			// get the loadcurve
			const char* szlc = tag.AttributeValue("lc");
			int lc = atoi(szlc) - 1;

			// get the (optional) type attribute
			bool brel = false;
			const char* szrel = tag.AttributeValue("type", true);
			if (szrel)
			{
				if (strcmp(szrel, "relative") == 0) brel = true;
				else if (strcmp(szrel, "absolute") == 0) brel = false;
				else throw XMLReader::InvalidAttributeValue(tag, "type", szrel);
			}

			// create the rigid displacement constraint
			FERigidBodyDisplacement* pDC = static_cast<FERigidBodyDisplacement*>(fecore_new<FERigidBC>("rigid_prescribed", &fem));
			rigid.AddPrescribedBC(pDC);

			pDC->SetID(nmat);
			pDC->SetBC(bc);
			pDC->SetRelativeFlag(brel);

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

			// add this boundary condition to the current step
			GetBuilder()->AddComponent(pDC);
		}
		else if (tag == "force")
		{
			// get the dof
			int bc = -1;
			const char* szbc = tag.AttributeValue("bc");
			if (strcmp(szbc, "x") == 0) bc = 0;
			else if (strcmp(szbc, "y") == 0) bc = 1;
			else if (strcmp(szbc, "z") == 0) bc = 2;
			else if (strcmp(szbc, "Rx") == 0) bc = 3;
			else if (strcmp(szbc, "Ry") == 0) bc = 4;
			else if (strcmp(szbc, "Rz") == 0) bc = 5;
			else throw XMLReader::InvalidAttributeValue(tag, "bc", szbc);

			// get the type
			int ntype = 0;
			bool bfollow = false;
			const char* sztype = tag.AttributeValue("type", true);
			if (sztype)
			{
				if (strcmp(sztype, "ramp") == 0) ntype = 1;
				else if (strcmp(sztype, "follow") == 0) bfollow = true;
				else throw XMLReader::InvalidAttributeValue(tag, "type", sztype);
			}

			// get the loadcurve
			const char* szlc = tag.AttributeValue("lc", true);
			int lc = -1;
			if (szlc) lc = atoi(szlc) - 1;

			// make sure there is a loadcurve for type=0 forces
			if ((ntype == 0) && (lc == -1)) throw XMLReader::MissingAttribute(tag, "lc");

			// create the rigid body force
			FERigidBodyForce* pFC = static_cast<FERigidBodyForce*>(fecore_new<FEModelLoad>(FEBC_ID, "rigid_force", &fem));
			pFC->SetType(ntype);
			pFC->SetID(nmat);
			pFC->SetBC(bc);
			pFC->SetFollowFlag(bfollow);

			double val = 0.0;
			value(tag, val);
			pFC->SetForce(val);

			if (lc >= 0)
			{
				FEParam* p = pFC->GetParameter("force");
				if (p == nullptr) throw XMLReader::InvalidTag(tag);
				GetFEModel()->AttachLoadController(p, lc);
			}

			// add it to the model
			GetBuilder()->AddModelLoad(pFC);
		}
		else if (tag == "fixed")
		{
			// get the dof
			int bc = -1;
			const char* szbc = tag.AttributeValue("bc");
			if (strcmp(szbc, "x") == 0) bc = 0;
			else if (strcmp(szbc, "y") == 0) bc = 1;
			else if (strcmp(szbc, "z") == 0) bc = 2;
			else if (strcmp(szbc, "Rx") == 0) bc = 3;
			else if (strcmp(szbc, "Ry") == 0) bc = 4;
			else if (strcmp(szbc, "Rz") == 0) bc = 5;
			else throw XMLReader::InvalidAttributeValue(tag, "bc", szbc);

			// create the fixed dof
			FERigidBodyFixedBC* pBC = static_cast<FERigidBodyFixedBC*>(fecore_new<FERigidBC>("rigid_fixed", &fem));
			pBC->id = nmat;
			pBC->bc = bc;
			rigid.AddFixedBC(pBC);

			// add this boundary condition to the current step
			GetBuilder()->AddComponent(pBC);
		}
		else if (tag == "initial_velocity")
		{
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
	} while (!tag.isend());
}
