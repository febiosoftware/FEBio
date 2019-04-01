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
#include "FEBioLoadDataSection.h"
#include <FECore/FEModel.h>
#include <FECore/FELoadController.h>
#include <FECore/FEPointFunction.h>

//-----------------------------------------------------------------------------
FEBioLoadDataSection3::FEBioLoadDataSection3(FEFileImport* pim) : FEFileSection(pim) 
{
}

//-----------------------------------------------------------------------------
//!  This function reads the load data section from the xml file
//!
void FEBioLoadDataSection3::Parse(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	
	++tag;
	do
	{
		if (tag == "load_controller")
		{
			// load curve ID
			int nid;
			tag.AttributeValue("id", nid);

			// get the type
			const char* sztype = tag.AttributeValue("type", true);
			if (sztype == 0) sztype = "loadcurve";

			// get the number of load curves
			int nlc = fem.LoadControllers();

			// check that the ID is one more than the number of load curves defined
			// This is to make sure that the ID's are in numerical order and no values are skipped.
			if (nid != nlc + 1) throw XMLReader::InvalidAttributeValue(tag, "id");

			// create the controller
			FELoadController* plc = fecore_new<FELoadController>(sztype, &fem); assert(plc);
			plc->SetID(nid - 1);

			// add the controller
			fem.AddLoadController(plc);

			// read the parameter list
			ReadParameterList(tag, plc);
		}
		else throw XMLReader::InvalidTag(tag);

		++tag;
	}
	while (!tag.isend());
}
