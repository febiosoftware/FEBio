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
};
