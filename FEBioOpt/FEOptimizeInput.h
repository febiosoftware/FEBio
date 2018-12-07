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
