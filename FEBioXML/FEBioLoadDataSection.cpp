#include "stdafx.h"
#include "FEBioLoadDataSection.h"
#include "FECore/FEModel.h"
#include "FECore/FEDataLoadCurve.h"

//-----------------------------------------------------------------------------
//!  This function reads the load data section from the xml file
//!
void FEBioLoadDataSection::Parse(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	int nmplc = -1;
	
	if (m_pim->m_nsteps > 0)
	{
		FEAnalysis* pstep = m_pim->GetStep();
		nmplc = pstep->m_nmplc+1;
	}

	++tag;
	do
	{
		if (tag == "loadcurve")
		{
			// load curve ID
			int nid;
			tag.AttributeValue("id", nid);

			// check that the ID is one more than the number of load curves defined
			// This is to make sure that the ID's are in numerical order and no values are skipped.
			int nlc = fem.LoadCurves();
			if (nid != nlc+1) throw XMLReader::InvalidAttributeValue(tag, "id");

			// default type and extend mode
			FEDataLoadCurve::INTFUNC ntype = FEDataLoadCurve::LINEAR;
			FEDataLoadCurve::EXTMODE nextm = FEDataLoadCurve::CONSTANT;

			// For backward compatibility we must set the default to STEP for the mustpoint loadcurve
			// Note that this assumes that the step has already been read in.
			if (nid == nmplc) ntype = FEDataLoadCurve::STEP;

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

			// count how many points we have
			int nlp = tag.children();

			// create the loadcurve
			FEDataLoadCurve* plc = dynamic_cast<FEDataLoadCurve*>(fecore_new<FELoadCurve>(FELOADCURVE_ID, "loadcurve", &fem));
			plc->SetInterpolation(ntype);
			plc->SetExtendMode(nextm);

			// read the points
			double d[2];
			++tag;
			for (int i=0; i<nlp; ++i)
			{
				tag.value(d, 2);
				plc->Add(d[0], d[1]);
				++tag;
			}

			// add the loadcurve
			fem.AddLoadCurve(plc);
		}
		else throw XMLReader::InvalidTag(tag);

		++tag;
	}
	while (!tag.isend());
}
