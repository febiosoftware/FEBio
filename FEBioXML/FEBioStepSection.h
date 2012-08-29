#pragma once
#include "FEBioImport.h"

//-----------------------------------------------------------------------------
// Step Section
class FEBioStepSection : public FEBioFileSection
{
public:
	FEBioStepSection(FEFEBioImport* pim) : FEBioFileSection(pim){}
	void Parse(XMLTag& tag);
};
