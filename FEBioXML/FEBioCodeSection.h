#pragma once
#include "FEBioImport.h"

//-----------------------------------------------------------------------------
// This is an experimental feature. 
// This section allows users to define callbacks from the input file. 
class FEBioCodeSection : public FEFileSection
{
public:
	FEBioCodeSection(FEFileImport* pim) : FEFileSection(pim){}
	void Parse(XMLTag& tag);
};
