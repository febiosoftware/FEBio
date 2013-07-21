#include "stdafx.h"
#include "FEBioLoadDataSection.h"

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

			// default type and extend mode
			FELoadCurve::INTFUNC ntype = FELoadCurve::LINEAR;
			FELoadCurve::EXTMODE nextm = FELoadCurve::CONSTANT;

			// For backward compatibility we must set the default to STEP for the mustpoint loadcurve
			// Note that this assumes that the step has already been read in.
			if (nid == nmplc) ntype = FELoadCurve::STEP;

			// get the (optional) type
			XMLAtt* patt = tag.Attribute("type", true);
			if (patt)
			{
				XMLAtt& type = *patt;
				if      (type == "step"  ) ntype = FELoadCurve::STEP;
				else if (type == "linear") ntype = FELoadCurve::LINEAR;
				else if (type == "smooth") ntype = FELoadCurve::SMOOTH;
				else throw XMLReader::InvalidAttributeValue(tag, "type", type.cvalue());
			}

			// get the optional extend mode
			patt = tag.Attribute("extend", true);
			if (patt)
			{
				XMLAtt& ext = *patt;
				if      (ext == "constant"     ) nextm = FELoadCurve::CONSTANT;
				else if (ext == "extrapolate"  ) nextm = FELoadCurve::EXTRAPOLATE;
				else if (ext == "repeat"       ) nextm = FELoadCurve::REPEAT;
				else if (ext == "repeat offset") nextm = FELoadCurve::REPEAT_OFFSET;
				else throw XMLReader::InvalidAttributeValue(tag, "extend", ext.cvalue());
			}

			// count how many points we have
			int nlp = tag.children();

			// create the loadcurve
			FELoadCurve* plc = new FELoadCurve;
			plc->Create(nlp);
			plc->SetInterpolation(ntype);
			plc->SetExtendMode(nextm);

			// read the points
			double d[2];
			++tag;
			for (int i=0; i<nlp; ++i)
			{
				tag.value(d, 2);
				plc->LoadPoint(i).time  = d[0];
				plc->LoadPoint(i).value = d[1];

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
