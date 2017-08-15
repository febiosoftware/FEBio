#pragma once
#include "FEBioImport.h"

//-----------------------------------------------------------------------------
// Step Section (old format)
class FEBioStepSection : public FEBioFileSection
{
public:
	FEBioStepSection(FEBioImport* pim) : FEBioFileSection(pim){}
	void Parse(XMLTag& tag);
};

//-----------------------------------------------------------------------------
// Step Section (2.0 format)
class FEBioStepSection2 : public FEBioFileSection
{
public:
	FEBioStepSection2(FEBioImport* pim) : FEBioFileSection(pim){}
	void Parse(XMLTag& tag);
};

//-----------------------------------------------------------------------------
// Step Section (2.5 format)
class FEBioStepSection25 : public FEFileSection
{
public:
	FEBioStepSection25(FEFileImport* pim) : FEFileSection(pim){}
	void Parse(XMLTag& tag);
};
