#pragma once
#include "FEBioImport.h"

//-----------------------------------------------------------------------------
class FEBioCodeSection : public FEBioFileSection
{
public:
	FEBioCodeSection(FEBioImport* pim) : FEBioFileSection(pim){}
	void Parse(XMLTag& tag);
};
