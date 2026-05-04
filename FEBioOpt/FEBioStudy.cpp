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
#include "FEBioStudy.h"
#include <FEBioXML/XMLReader.h>
#include <FEBioXML/xmltool.h>
#include <FECore/log.h>
#include "FEBioReport.h"
#include <filesystem>

FEBioStudy::FEBioStudy(FEModel* fem) : FECoreTask(fem) {}

FEBioStudy::~FEBioStudy() {}

bool FEBioStudy::Init(const char* szfile)
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

bool FEBioStudy::ReadControlFile(const char* szfile)
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
		if (!xml.FindTag(sztype, tag)) {
			feLogError("Failed finding root tag '%s'.", sztype);
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

bool FEBioStudy::Run()
{
	// execute the study
	m_success = Execute();

	// build the report
	FEBioReport report;
	report.SetTitle(GetName());
	report.SetOptionsFile(m_optionsFile);
	report.SetStatus(m_success ? 1 : 0);
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

void FEBioStudy::BuildReport(FEBioReport& report)
{
	report.AddSection("Summary").AddText(m_success ? "Study completed successfully." : "Study failed.");
}
