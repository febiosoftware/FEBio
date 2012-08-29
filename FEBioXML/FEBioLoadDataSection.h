#pragma once
#include "FEBioImport.h"

//-----------------------------------------------------------------------------
// LoadData Section
class FEBioLoadDataSection : public FEBioFileSection
{
public:
	FEBioLoadDataSection(FEFEBioImport* pim) : FEBioFileSection(pim) {}
	void Parse(XMLTag& tag);
};
