#include "stdafx.h"
#include "FEBioLoadDataSection.h"
#include <FECore/FEModel.h>
#include <FECore/LoadCurve.h>
#include <FECore/FEPointFunction.h>

//-----------------------------------------------------------------------------
FEBioLoadDataSection::FEBioLoadDataSection(FEFileImport* pim) : FEFileSection(pim) 
{
	m_redefineCurves = false;
}

//-----------------------------------------------------------------------------
//!  This function reads the load data section from the xml file
//!
void FEBioLoadDataSection::Parse(XMLTag& tag)
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

			// default type and extend mode
			FEPointFunction::INTFUNC ntype = FEPointFunction::LINEAR;
			FEPointFunction::EXTMODE nextm = FEPointFunction::CONSTANT;

			// get the (optional) type
			XMLAtt* patt = tag.Attribute("type", true);
			if (patt)
			{
				XMLAtt& type = *patt;
				if      (type == "step"  ) ntype = FEPointFunction::STEP;
				else if (type == "linear") ntype = FEPointFunction::LINEAR;
				else if (type == "smooth") ntype = FEPointFunction::SMOOTH;
				else throw XMLReader::InvalidAttributeValue(tag, "type", type.cvalue());
			}

			// get the optional extend mode
			patt = tag.Attribute("extend", true);
			if (patt)
			{
				XMLAtt& ext = *patt;
				if      (ext == "constant"     ) nextm = FEPointFunction::CONSTANT;
				else if (ext == "extrapolate"  ) nextm = FEPointFunction::EXTRAPOLATE;
				else if (ext == "repeat"       ) nextm = FEPointFunction::REPEAT;
				else if (ext == "repeat offset") nextm = FEPointFunction::REPEAT_OFFSET;
				else throw XMLReader::InvalidAttributeValue(tag, "extend", ext.cvalue());
			}

			// get the number of load curves
			int nlc = fem.LoadCurves();

			// find or create the load curve
			FEPointFunction* pfnc = 0;

			// see if this refers to a valid curve
			if (m_redefineCurves)
			{
				if ((nid > 0) && (nid <= nlc)) 
				{
					FELoadCurve* plc = fem.GetLoadCurve(nid - 1);
					pfnc = dynamic_cast<FEPointFunction*>(plc);
					assert(pfnc);

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

			// set the load curve attributes
			pfnc->SetInterpolation(ntype);
			pfnc->SetExtendMode(nextm);

			// read the data points
			double d[2];
			++tag;
			do
			{
				tag.value(d, 2);
				pfnc->Add(d[0], d[1]);
				++tag;
			}
			while (!tag.isend());
		}
		else throw XMLReader::InvalidTag(tag);

		++tag;
	}
	while (!tag.isend());
}
