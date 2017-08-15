#pragma once
#include "FEBioImport.h"

//-----------------------------------------------------------------------------
// Initial Section
class FEBioInitialSection : public FEFileSection
{
public:
	FEBioInitialSection(FEFileImport* pim) : FEFileSection(pim){}
	void Parse(XMLTag& tag);
};

//-----------------------------------------------------------------------------
// Initial Section
class FEBioInitialSection25 : public FEFileSection
{
public:
	FEBioInitialSection25(FEFileImport* pim) : FEFileSection(pim){}
	void Parse(XMLTag& tag);
};
