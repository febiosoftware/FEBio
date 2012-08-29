#pragma once
#include "FEBioImport.h"

//-----------------------------------------------------------------------------
// FEBio Module Section
class FEBioModuleSection : public FEBioFileSection
{
public:
	FEBioModuleSection(FEFEBioImport* pim) : FEBioFileSection(pim) {}
	void Parse(XMLTag& tag);
};
