#include "stdafx.h"
#include "FEBioLoadDataSection.h"
#include <FECore/FEModel.h>
#include <FECore/LoadCurve.h>
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
		if (tag == "loadcurve")
		{
			// load curve ID
			int nid;
			tag.AttributeValue("id", nid);

			// get the number of load curves
			int nlc = fem.LoadCurves();

			// find or create the load curve
			FEFunction1D* pfnc = 0;

			// see if this refers to a valid curve
			if (m_redefineCurves)
			{
				if ((nid > 0) && (nid <= nlc)) 
				{
					FELoadCurve* plc = fem.GetLoadCurve(nid - 1);
					pfnc = plc->GetFunction();

					// clear the curve since we're about to read in new data points
					pfnc->Clear();
				}
			}

			// if the ID does not refer to an existing curve, make sure it defines the next curve
			if (pfnc == 0)
			{
				// check that the ID is one more than the number of load curves defined
				// This is to make sure that the ID's are in numerical order and no values are skipped.
				if (nid != nlc + 1) throw XMLReader::InvalidAttributeValue(tag, "id");

				// create the loadcurve
				pfnc = fecore_new<FEPointFunction>("point", &fem); assert(pfnc);

				// add the loadcurve
				fem.AddLoadCurve(new FELoadCurve(pfnc));
			}

			// read the parameter list
			ReadParameterList(tag, pfnc);
		}
		else throw XMLReader::InvalidTag(tag);

		++tag;
	}
	while (!tag.isend());
}
