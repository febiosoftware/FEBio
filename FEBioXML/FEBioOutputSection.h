#pragma once
#include "FEBioImport.h"

//-----------------------------------------------------------------------------
// Output Section
class FEBioOutputSection : public FEBioFileSection
{
public:
	FEBioOutputSection(FEFEBioImport* pim) : FEBioFileSection(pim){}
	void Parse(XMLTag& tag);

protected:
	void ParseLogfile (XMLTag& tag);
	void ParsePlotfile(XMLTag& tag);
};
