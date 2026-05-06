/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2026 University of Utah, The Trustees of Columbia University in
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
#include "FECoreStudy.h"
#include "XMLReader.h"
#include "xmltool.h"
#include "log.h"
#include "FEBioReport.h"
#include "FEModel.h"
#include "FEAnalysis.h"
#include <filesystem>

FECoreStudy::FECoreStudy(FEModel* fem) : FECoreTask(fem) {}

FECoreStudy::~FECoreStudy() {}

bool FECoreStudy::Init(const char* szfile)
{
	if (szfile == nullptr) {
		feLogError("No control file specified.");
		return false;
	}

	// set the file name as the study's name
	m_optionsFile = szfile;
	std::string name = std::filesystem::path(szfile).stem().string();
	SetName(name);

	// read the control file
	if (ReadControlFile(szfile) == false) return false;

	return Init();
}

void FECoreStudy::DisableStandardOutput()
{
	// suppress all standard output from FEBio (we assume that the study will report its results in a different way)
	// TODO: how do I block output to log file? I can use BlockLog, but I think that blocks everything. 
	FEModel* fem = GetFEModel();
	for (int i = 0; i < fem->Steps(); ++i)
	{
		fem->GetStep(i)->SetPlotLevel(FE_PLOT_NEVER);
	}
}

bool FECoreStudy::ReadControlFile(const char* szfile)
{
	XMLReader xml;
	if (xml.Open(szfile) == false) return false;
	try {
		// find the root tag (should be the same as the type string)
		const char* sztype = GetTypeStr();
		if (sztype == nullptr) {
			feLogError("The study class does not have a type string.");
			return false;
		}
		XMLTag tag;
		if (!xml.FindTag("febio_study", tag)) {
			feLogError("Failed finding root tag \"febio_study\".");
			return false;
		}

		// make sure the type matches
		const char* szval = tag.AttributeValue("type", true);
		if (szval == nullptr) {
			feLogError("The root tag is missing the required \"type\" attribute.");
			return false;
		}

		if (strcmp(szval, sztype) != 0) {
			feLogError("The study type \"%s\" does not match the expected type \"%s\".", szval, sztype);
			return false;
		}

		// read the data
		if (fexml::readParameterList(tag, this) == false)
		{
			feLogError("Failed reading study parameters.");
			return false;
		}
	}
	catch (XMLReader::Error& e) {
		feLogError(e.what());
		return false;
	}
	catch (std::exception& e) {
		feLogError(e.what());
		return false;
	}
	catch (...) {
		feLogError("An unknown error occurred while loading the options file.");
		return false;
	}

	xml.Close();
	return true;
}

bool FECoreStudy::Run()
{
	// execute the study
	m_success = Execute();

	// create the report
	FEBioReport report;
	report.SetTitle(GetName());
	report.SetOptionsFile(m_optionsFile);
	report.SetStatus(m_success ? 1 : 0);

	// build the report (this is overridden by derived classes to add content to the report)
	BuildReport(report);

	// write the report to a file
	string name = GetName();
	std::string filename = name + ".febr";
	bool b = report.Write(filename);
	if (b) 
		feLog("Study report written to %s", filename.c_str());
	else
		feLogError("Failed writing study report.");

	// all done!
	return m_success;
}

void FECoreStudy::BuildReport(FEBioReport& report)
{
	report.AddSection("Summary").AddText(m_success ? "Study completed successfully." : "Study failed.");
}
