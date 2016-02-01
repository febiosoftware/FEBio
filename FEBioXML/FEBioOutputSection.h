#pragma once
#include "FEBioImport.h"

class FEFacetSet;

//-----------------------------------------------------------------------------
// Output Section
class FEBioOutputSection : public FEBioFileSection
{
public:
	FEBioOutputSection(FEBioImport* pim) : FEBioFileSection(pim){}
	void Parse(XMLTag& tag);

protected:
	void ParseLogfile (XMLTag& tag);
	void ParsePlotfile(XMLTag& tag);
    
protected:
    bool BuildSurface(FESurface& s, FEFacetSet& fs);
};
