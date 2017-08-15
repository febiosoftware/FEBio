#include "stdafx.h"
#include "FEBioLoadDataSection.h"
#include "FECore/FEModel.h"
#include "FECore/FEDataLoadCurve.h"

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
			FEDataLoadCurve::INTFUNC ntype = FEDataLoadCurve::LINEAR;
			FEDataLoadCurve::EXTMODE nextm = FEDataLoadCurve::CONSTANT;

			// get the (optional) type
			XMLAtt* patt = tag.Attribute("type", true);
			if (patt)
			{
				XMLAtt& type = *patt;
				if      (type == "step"  ) ntype = FEDataLoadCurve::STEP;
				else if (type == "linear") ntype = FEDataLoadCurve::LINEAR;
				else if (type == "smooth") ntype = FEDataLoadCurve::SMOOTH;
				else throw XMLReader::InvalidAttributeValue(tag, "type", type.cvalue());
			}

			// get the optional extend mode
			patt = tag.Attribute("extend", true);
			if (patt)
			{
				XMLAtt& ext = *patt;
				if      (ext == "constant"     ) nextm = FEDataLoadCurve::CONSTANT;
				else if (ext == "extrapolate"  ) nextm = FEDataLoadCurve::EXTRAPOLATE;
				else if (ext == "repeat"       ) nextm = FEDataLoadCurve::REPEAT;
				else if (ext == "repeat offset") nextm = FEDataLoadCurve::REPEAT_OFFSET;
				else throw XMLReader::InvalidAttributeValue(tag, "extend", ext.cvalue());
			}

			// get the number of load curves
			int nlc = fem.LoadCurves();

			// find or create the load curve
			FEDataLoadCurve* plc = 0;

			// see if this refers to a valid curve
			if (m_redefineCurves)
			{
				if ((nid > 0) && (nid <= nlc)) 
				{
					plc = dynamic_cast<FEDataLoadCurve*>(fem.GetLoadCurve(nid - 1));
					assert(plc);

					// clear the curve since we're about to read in new data points
					plc->Clear();
				}
			}

			// if the ID does not refer to an existing curve, make sure it defines the next curve
			if (plc == 0)
			{
				// check that the ID is one more than the number of load curves defined
				// This is to make sure that the ID's are in numerical order and no values are skipped.
				if (nid != nlc + 1) throw XMLReader::InvalidAttributeValue(tag, "id");

				// create the loadcurve
				plc = dynamic_cast<FEDataLoadCurve*>(fecore_new<FELoadCurve>(FELOADCURVE_ID, "loadcurve", &fem));

				// add the loadcurve
				fem.AddLoadCurve(plc);
			}

			// set the load curve attributes
			plc->SetInterpolation(ntype);
			plc->SetExtendMode(nextm);

			// read the data points
			double d[2];
			++tag;
			do
			{
				tag.value(d, 2);
				plc->Add(d[0], d[1]);
				++tag;
			}
			while (!tag.isend());
		}
		else throw XMLReader::InvalidTag(tag);

		++tag;
	}
	while (!tag.isend());
}
