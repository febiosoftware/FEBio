#pragma once
#include "FEBioImport.h"

//-----------------------------------------------------------------------------
// Parameters section
// Allows users to define parameters in the input file
class FEBioParametersSection : public FEFileSection
{
public:
	FEBioParametersSection(FEFileImport* pim) : FEFileSection(pim){}
	void Parse(XMLTag& tag);
};
