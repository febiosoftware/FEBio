// FERestartImport.cpp: implementation of the FERestartImport class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FERestartImport.h"
#include "FECore/FESolver.h"
#include "FECore/FEAnalysis.h"
#include "FECore/FEModel.h"
#include "FECore/DumpFile.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

FERestartImport::FERestartImport()
{
	m_pfem = 0;
}

FERestartImport::~FERestartImport()
{
	m_pfem = 0;
}

//-----------------------------------------------------------------------------

bool FERestartImport::Load(FEModel& fem, const char* szfile)
{
	// open the XML file
	if (m_xml.Open(szfile) == false) return errf("FATAL ERROR: Failed opening restart file %s\n", szfile);

	// keep a pointer to the fem object
	m_pfem = &fem;

	m_szdmp[0] = 0;

	// loop over child tags
	try
	{
		// find the root element
		XMLTag tag;
		if (m_xml.FindTag("febio_restart", tag) == false) return errf("FATAL ERROR: File does not contain restart data.\n");

		// check the version number
		if (strcmp(tag.m_att[0].m_szatv, "1.0") != 0) return errf("FATAL ERROR: Incorrect restart file version\n");

		// the next tag has to be the archive
		++tag;
		if (tag != "Archive") return errf("FATAL ERROR: The first element must be the archive name\n");
		char szar[256];
		tag.value(szar);

		// open the archive
		DumpFile ar(fem);
		if (ar.Open(szar) == false) return errf("FATAL ERROR: failed opening restart archive\n");

		// read the archive
		if (fem.Serialize(ar) == false) return errf("FATAL ERROR: failed reading restart data from archive %s\n", szar);

		// read the restart data
		++tag;
		while (!tag.isend())
		{
			if      (tag == "Control") ParseControlSection  (tag);
			else if (tag == "LoadData") ParseLoadSection    (tag);
			else throw XMLReader::InvalidTag(tag);

			++tag;
		}
	}
	catch (XMLReader::Error& e)
	{
		fprintf(stderr, "FATAL ERROR: %s (line %d)\n", e.GetErrorString(), m_xml.GetCurrentLine());
		return false;
	}
	catch (...)
	{
		fprintf(stderr, "FATAL ERROR: unrecoverable error (line %d)\n", m_xml.GetCurrentLine());
		return false;
	}

	// close the XML file
	m_xml.Close();

	// we're done!
	return true;
}

//-----------------------------------------------------------------------------
//! This function reads the load data section from the restart file.

bool FERestartImport::ParseLoadSection(XMLTag& tag)
{
	FEModel& fem = *m_pfem;

	++tag;
	do
	{
		if (tag == "loadcurve")
		{
			// get the ID from the curve
			const char* szid = tag.AttributeValue("id");
			int nid = atoi(szid);

			// find the loadcurve with this ID
			FELoadCurve* plc = fem.GetLoadCurve(nid-1);

			// count how many points we have
			XMLTag t(tag); ++t;
			int nlp = 0;
			while (!t.isend()) { ++nlp; ++t; }

			// create the loadcurve
			plc->Create(nlp);

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
		}
		else throw XMLReader::InvalidTag(tag);

		++tag;
	}
	while (!tag.isend());

	return true;
}

//-----------------------------------------------------------------------------
//!  This function parses the control section from the restart file

bool FERestartImport::ParseControlSection(XMLTag& tag)
{
	FEModel& fem = *m_pfem;
	FEAnalysis* pstep = fem.GetCurrentStep();

	++tag;
	do
	{
		if      (tag == "time_steps"        ) tag.value(pstep->m_ntime);
		if      (tag == "final_time"        ) tag.value(pstep->m_final_time);
		else if (tag == "step_size"         ) tag.value(pstep->m_dt0);
//		else if (tag == "dtol"              ) tag.value(pstep->m_psolver->m_Dtol);
//		else if (tag == "ptol"              ) tag.value(pstep->m_psolver->m_Ptol);
//		else if (tag == "ctol"              ) tag.value(pstep->m_psolver->m_Ctol);
//		else if (tag == "etol"              ) tag.value(pstep->m_psolver->m_Etol);
//		else if (tag == "rtol"              ) tag.value(pstep->m_psolver->m_Rtol);
//		else if (tag == "lstol"             ) tag.value(pstep->m_psolver->m_pbfgs->m_LStol);
//		else if (tag == "max_refs"          ) tag.value(pstep->m_psolver->m_pbfgs->m_maxref);
//		else if (tag == "max_ups"           ) tag.value(pstep->m_psolver->m_pbfgs->m_maxups);
//		else if (tag == "cmax"              ) tag.value(pstep->m_psolver->m_pbfgs->m_cmax);
		else if (tag == "pressure_stiffness") tag.value(pstep->m_istiffpr);
		else if (tag == "restart" ) 
		{
			const char* szf = tag.AttributeValue("file", true);
			if (szf) strcpy(m_szdmp, szf);
			char szval[256];
			tag.value(szval);
			if		(strcmp(szval, "DUMP_DEFAULT"    ) == 0) {} // don't change the restart level
			else if (strcmp(szval, "DUMP_NEVER"      ) == 0) pstep->SetDumpLevel(FE_DUMP_NEVER);
			else if (strcmp(szval, "DUMP_MAJOR_ITRS" ) == 0) pstep->SetDumpLevel(FE_DUMP_MAJOR_ITRS);
			else if (strcmp(szval, "DUMP_STEP"       ) == 0) pstep->SetDumpLevel(FE_DUMP_STEP);
			else if (strcmp(szval, "0" ) == 0) pstep->SetDumpLevel(FE_DUMP_NEVER);		// for backward compatibility only
			else if (strcmp(szval, "1" ) == 0) pstep->SetDumpLevel(FE_DUMP_MAJOR_ITRS); // for backward compatibility only
			else throw XMLReader::InvalidValue(tag);
		}
		else if (tag == "time_stepper")
		{
			pstep->m_bautostep = true;
			++tag;
			do
			{
				if      (tag == "max_retries") tag.value(pstep->m_maxretries);
				else if (tag == "opt_iter"   ) tag.value(pstep->m_iteopt);
				else if (tag == "dtmin"      ) tag.value(pstep->m_dtmin);
				else throw XMLReader::InvalidTag(tag);

				++tag;
			}
			while (!tag.isend());
		}
		else if (tag == "plot_level")
		{
			char szval[256];
			tag.value(szval);
			if      (strcmp(szval, "PLOT_NEVER"        ) == 0) pstep->SetPlotLevel(FE_PLOT_NEVER);
			else if (strcmp(szval, "PLOT_MAJOR_ITRS"   ) == 0) pstep->SetPlotLevel(FE_PLOT_MAJOR_ITRS);
			else if (strcmp(szval, "PLOT_MINOR_ITRS"   ) == 0) pstep->SetPlotLevel(FE_PLOT_MINOR_ITRS);
			else if (strcmp(szval, "PLOT_MUST_POINTS"  ) == 0) pstep->SetPlotLevel(FE_PLOT_MUST_POINTS);
			else if (strcmp(szval, "PLOT_FINAL"        ) == 0) pstep->SetPlotLevel(FE_PLOT_FINAL);
			else if (strcmp(szval, "PLOT_AUGMENTATIONS") == 0) pstep->SetPlotLevel(FE_PLOT_AUGMENTATIONS);
			else throw XMLReader::InvalidValue(tag);
		}
		else throw XMLReader::InvalidTag(tag);

		++tag;
	}
	while (!tag.isend());

	return true;
}
