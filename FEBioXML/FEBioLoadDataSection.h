#pragma once
#include "FEBioImport.h"

//-----------------------------------------------------------------------------
// LoadData Section
class FEBioLoadDataSection : public FEFileSection
{
public:
	FEBioLoadDataSection(FEFileImport* pim) : FEFileSection(pim) {}
	void Parse(XMLTag& tag);
};
