#pragma once
#include "FEBioImport.h"

//-----------------------------------------------------------------------------
// FEBio Module Section
class FEBioModuleSection : public FEBioFileSection
{
public:
	FEBioModuleSection(FEBioImport* pim) : FEBioFileSection(pim) {}
	void Parse(XMLTag& tag);
};
