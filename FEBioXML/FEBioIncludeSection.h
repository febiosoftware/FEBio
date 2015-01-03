#pragma once
#include "FEBioImport.h"

//-----------------------------------------------------------------------------
// Include section (new in version 2.0)
class FEBioIncludeSection : public FEBioFileSection
{
public:
	FEBioIncludeSection(FEFEBioImport* pim) : FEBioFileSection(pim){}
	void Parse(XMLTag& tag);
};
