#pragma once
#include "FEBioImport.h"

//-----------------------------------------------------------------------------
// Control Section
class FEBioControlSection : public FEBioFileSection
{
public:
	FEBioControlSection(FEBioImport* pim) : FEBioFileSection(pim) {}
	void Parse(XMLTag& tag);

protected:
	bool ParseCommonParams(XMLTag& tag);
	void ParseIntegrationRules(XMLTag& tag);
};

//-----------------------------------------------------------------------------
// Control Section for steps
class FEStepControlSection : public FEFileSection
{
public:
	FEStepControlSection(FEFileImport* pim) : FEFileSection(pim) {}
	void Parse(XMLTag& tag);

protected:
	bool ParseCommonParams(XMLTag& tag);
};

//-----------------------------------------------------------------------------
// Control Section
class FEBioControlSection3 : public FEBioFileSection
{
public:
	FEBioControlSection3(FEBioImport* pim) : FEBioFileSection(pim) {}
	void Parse(XMLTag& tag);

protected:
	bool ParseCommonParams(XMLTag& tag);
	void ParseIntegrationRules(XMLTag& tag);
};
