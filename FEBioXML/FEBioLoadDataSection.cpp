#include "stdafx.h"
#include "FEBioLoadDataSection.h"

//-----------------------------------------------------------------------------
//!  This function reads the load data section from the xml file
//!
void FEBioLoadDataSection::Parse(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	int nmplc = -1;
	if (m_pim->Version() <= 0x0100)
	{
		FEAnalysisStep* pstep = GetStep();
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
			const char* szt = tag.AttributeValue("type", true);
			if (szt)
			{
				if      (strcmp(szt, "step"  ) == 0) ntype = FELoadCurve::STEP;
				else if (strcmp(szt, "linear") == 0) ntype = FELoadCurve::LINEAR;
				else if (strcmp(szt, "smooth") == 0) ntype = FELoadCurve::SMOOTH;
				else throw XMLReader::InvalidAttributeValue(tag, "type", szt);
			}

			// get the optional extend mode
			const char* szm = tag.AttributeValue("extend", true);
			if (szm)
			{
				if      (strcmp(szm, "constant"     ) == 0) nextm = FELoadCurve::CONSTANT;
				else if (strcmp(szm, "extrapolate"  ) == 0) nextm = FELoadCurve::EXTRAPOLATE;
				else if (strcmp(szm, "repeat"       ) == 0) nextm = FELoadCurve::REPEAT;
				else if (strcmp(szm, "repeat offset") == 0) nextm = FELoadCurve::REPEAT_OFFSET;
				else throw XMLReader::InvalidAttributeValue(tag, "extend", szt);
			}

			// count how many points we have
			XMLTag t(tag); ++t;
			int nlp = 0;
			while (!t.isend()) { ++nlp; ++t; }

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
