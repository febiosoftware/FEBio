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
#pragma once
#include "FECoreTask.h"
#include "FEBioReport.h"
#include "fecore_api.h"

class FEBioReport;

// A study is a task that has a assumed structure. 
// It requires a separate input file that can be processed with FEBioXML.
class FECORE_API FECoreStudy : public FECoreTask
{
	FECORE_BASE_CLASS(FECoreStudy)

public:
	FECoreStudy(FEModel* fem);
	~FECoreStudy();

	bool Init(const char* szfile) final;

	bool Run() final;

	using FECoreBase::Init;

	virtual void BuildReport(FEBioReport& report);

	std::string GetOptionsFile() const { return m_optionsFile; }

	// Disable the standard output from FEBio. 
	// This is useful for studies that want to report their results in a different way.
	// (Make sure to call this before model initialization!)
	void DisableStandardOutput();

public:
	virtual bool Execute() = 0;

private:
	bool ReadControlFile(const char* szfile);

private:
	std::string m_optionsFile;
	bool m_success = false;
};
