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
#include "FEBioLoadDataSection.h"
#include <FECore/FEModel.h>
#include <FECore/FELoadController.h>
#include <FECore/FEPointFunction.h>

//-----------------------------------------------------------------------------
FEBioLoadDataSection3::FEBioLoadDataSection3(FEFileImport* pim) : FEFileSection(pim) 
{
	m_redefineCurves = false;
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
			const char* sztype = tag.AttributeValue("type");

			// get the number of load curves
			int nlc = fem.LoadControllers();

			// create the controller
			FELoadController* plc = fecore_new<FELoadController>(sztype, &fem);
			if (plc == nullptr) throw XMLReader::InvalidAttributeValue(tag, "type");

			// set the ID
			plc->SetID(nid - 1);

			// see if this refers to a valid curve
			if (m_redefineCurves && ((nid > 0) && (nid <= nlc)))
			{
				fem.ReplaceLoadController(nid - 1, plc);
			}
			else
			{
				// check that the ID is one more than the number of load curves defined
				// This is to make sure that the ID's are in numerical order and no values are skipped.
				if (nid != nlc + 1)
				{
					delete plc;
					throw XMLReader::InvalidAttributeValue(tag, "id");
				}

				// add the controller
				fem.AddLoadController(plc);
			}

			// read the parameter list
			ReadParameterList(tag, plc);

			if (m_redefineCurves)
			{
				// We only get here during restart, in which case the new load controllers
				// will not get a chance to be initialized, so we'll do it here.
				plc->Init();
			}
		}
		else throw XMLReader::InvalidTag(tag);

		++tag;
	}
	while (!tag.isend());
}
