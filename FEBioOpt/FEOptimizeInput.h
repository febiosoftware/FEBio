/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2020 University of Utah, The Trustees of Columbia University in 
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
#include <FEBioXML/XMLReader.h>
#include <FEBioXML/FileImport.h>
#include <FECore/FEModel.h>

//-----------------------------------------------------------------------------
//! FEBio error terminated during the optimization
class FEErrorTermination{};

//-----------------------------------------------------------------------------
class FEOptimizeData;

//=============================================================================
//! Class that reads the optimization input file
class FEOptimizeInput : public FEFileImport
{
public:
	bool Input(const char* szfile, FEOptimizeData* pOpt);

private:
	FEOptimizeData*	m_opt;
};

//=============================================================================
class FEOptionsSection : public FEFileSection
{
public:
	FEOptionsSection(FEOptimizeData* opt, FEFileImport* im) : FEFileSection(im), m_opt(opt) {}
	void Parse(XMLTag& tag) override;

private:
	FEOptimizeData*	m_opt;
};

//=============================================================================
class FETaskSection : public FEFileSection
{
public:
	FETaskSection(FEOptimizeData* opt, FEFileImport* im) : FEFileSection(im), m_opt(opt) {}
	void Parse(XMLTag& tag) override;

private:
	FEOptimizeData*	m_opt;
};

//=============================================================================
class FEObjectiveSection : public FEFileSection
{
public:
	FEObjectiveSection(FEOptimizeData* opt, FEFileImport* im) : FEFileSection(im), m_opt(opt) {}
	void Parse(XMLTag& tag) override;

private:
	FEDataSource* ParseDataSource(XMLTag& tag, FEOptimizeData& opt);

private:
	FEOptimizeData*	m_opt;
};

//=============================================================================
class FEParametersSection : public FEFileSection
{
public:
	FEParametersSection(FEOptimizeData* opt, FEFileImport* im) : FEFileSection(im), m_opt(opt) {}
	void Parse(XMLTag& tag) override;

private:
	FEOptimizeData*	m_opt;
};

//=============================================================================
class FEConstraintsSection : public FEFileSection
{
public:
	FEConstraintsSection(FEOptimizeData* opt, FEFileImport* im) : FEFileSection(im), m_opt(opt) {}
	void Parse(XMLTag& tag) override;

private:
	FEOptimizeData*	m_opt;
};
