#pragma once
#include "FEBioImport.h"

//-----------------------------------------------------------------------------
// Step Section
class FEBioStepSection : public FEBioFileSection
{
public:
	FEBioStepSection(FEBioImport* pim) : FEBioFileSection(pim){}
	void Parse(XMLTag& tag);
};
