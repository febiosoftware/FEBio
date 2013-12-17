#pragma once
#include "FEBioImport.h"

//-----------------------------------------------------------------------------
class FEBioDiscreteSection : public FEBioFileSection
{
public:
	FEBioDiscreteSection(FEFEBioImport* pim) : FEBioFileSection(pim){}
	void Parse(XMLTag& tag);

protected:
	void ParseSpringSection(XMLTag& tag);
};
