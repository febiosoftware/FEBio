/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
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
#include "FECore/FEModel.h"
#include "FEBioXML/FEBioImport.h"

class FEDiagnostic;

//-----------------------------------------------------------------------------
class FEDiagnosticScenario : public FEParamContainer
{
public:
	FEDiagnosticScenario(FEDiagnostic* pdia) : m_pdia(pdia) {};

	FEDiagnostic* GetDiagnostic() { return m_pdia; }

	virtual bool Init() { return true; }

private:
	FEDiagnostic* m_pdia;
};

//-----------------------------------------------------------------------------
//! The FEDiagnostic class is a base class that can be used to create
//! diagnostic classes to test FEBio's performance.

class FEDiagnostic : public FECoreClass
{
public:
	FECORE_BASE_CLASS(FEDiagnostic)

public:
	//! constructor
	FEDiagnostic(FEModel* fem);

	//! destructor
	virtual ~FEDiagnostic();

	//! initialization
	virtual bool Init() { return true; }

	//! run the diagnostic. Returns true on pass, false on failure
	virtual bool Run() = 0;

	//! load data from file
	virtual bool ParseSection(XMLTag& tag) { return false; }

	//! create a scenario class
	virtual FEDiagnosticScenario* CreateScenario(const std::string& sname) { return 0; }

	void SetFileName(const std::string& fileName);
	const std::string& GetFileName();

private:
	std::string	m_file;	//!< the input file used
};

//-----------------------------------------------------------------------------
// Control Section
class FEDiagnosticControlSection : public FEFileSection
{
public:
	FEDiagnosticControlSection(FEFileImport* pim) : FEFileSection(pim) {}
	void Parse(XMLTag& tag);
};

//-----------------------------------------------------------------------------
// Scenario Section parser
class FEDiagnosticScenarioSection : public FEFileSection
{
public:
	FEDiagnosticScenarioSection(FEFileImport* pim) : FEFileSection(pim){}
	void Parse(XMLTag& tag);
};

//-----------------------------------------------------------------------------
//! The FEDiagnosticImport class creates a specific diagnostic test. Currently
//! the only way to create a diagnostic is to load a diagnostic from file

class FEDiagnosticImport : public FEFileImport
{
public:
	FEDiagnostic* LoadFile(FEModel& fem, const char* szfile);

protected:
	FEDiagnostic* m_pdia;

	friend class FEDiagnosticScenarioSection;
};
