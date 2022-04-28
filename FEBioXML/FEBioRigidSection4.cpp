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
#include "FEBioRigidSection4.h"
#include <FECore/FEModel.h>
#include <FECore/FECoreKernel.h>
#include <FECore/FEModelComponent.h>
#include <FECore/FEModelLoad.h>
#include <FECore/FENLConstraint.h>
#include <FECore/FEBoundaryCondition.h>
#include <FECore/FEInitialCondition.h>

void FEBioRigidSection4::Parse(XMLTag& tag)
{
	if (tag.isleaf()) return;

	++tag;
	do
	{
		if      (tag == "rigid_bc"       ) ParseRigidBC(tag);
		else if (tag == "rigid_ic"       ) ParseRigidIC(tag);
		else if (tag == "rigid_load"     ) ParseRigidLoad(tag);
		else if (tag == "rigid_connector") ParseRigidConnector(tag);
		else throw XMLReader::InvalidTag(tag);
		++tag;
	} while (!tag.isend());
}

void FEBioRigidSection4::ParseRigidBC(XMLTag& tag)
{
	FEModel* fem = GetFEModel();
	FEModelBuilder& feb = *GetBuilder();

	// get the name
	const char* szname = tag.AttributeValue("name", true);

	// get the type
	const char* sztype = tag.AttributeValue("type");

	FEBoundaryCondition* bc = fecore_new<FEBoundaryCondition>(sztype, fem);
	if (szname) bc->SetName(szname);

	GetBuilder()->AddRigidComponent(bc);
	ReadParameterList(tag, bc);
}

void FEBioRigidSection4::ParseRigidIC(XMLTag& tag)
{
	FEModel* fem = GetFEModel();

	// get the name
	const char* szname = tag.AttributeValue("name", true);

	// get the type
	const char* sztype = tag.AttributeValue("type");

	FEInitialCondition* ic = fecore_new<FEInitialCondition>(sztype, fem);
	if (szname) ic->SetName(szname);

	GetBuilder()->AddRigidComponent(ic);
	ReadParameterList(tag, ic);
}

void FEBioRigidSection4::ParseRigidLoad(XMLTag& tag)
{
	// get the name
	const char* szname = tag.AttributeValue("name", true);

	// get the type
	const char* sztype = tag.AttributeValue("type");

	FEModelLoad* ml = fecore_new<FEModelLoad>(sztype, GetFEModel());
	if (szname) ml->SetName(szname);

	GetBuilder()->AddRigidComponent(ml);
	ReadParameterList(tag, ml);
}

void FEBioRigidSection4::ParseRigidConnector(XMLTag& tag)
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
