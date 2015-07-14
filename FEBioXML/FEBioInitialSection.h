#pragma once
#include "FEBioImport.h"

//-----------------------------------------------------------------------------
// Initial Section
class FEBioInitialSection : public FEBioFileSection
{
public:
	FEBioInitialSection(FEBioImport* pim) : FEBioFileSection(pim){}
	void Parse(XMLTag& tag);
};
