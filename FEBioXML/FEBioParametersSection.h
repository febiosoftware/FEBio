#pragma once
#include "FEBioImport.h"

//-----------------------------------------------------------------------------
// Parameters section (new in version 2.0)
class FEBioParametersSection : public FEBioFileSection
{
public:
	FEBioParametersSection(FEBioImport* pim) : FEBioFileSection(pim){}
	void Parse(XMLTag& tag);
};
