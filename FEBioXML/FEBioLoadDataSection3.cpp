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
