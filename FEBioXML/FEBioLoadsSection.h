#pragma once
#include "FEBioImport.h"

//-----------------------------------------------------------------------------
// Loads Section (new in version 1.2)
class FEBioLoadsSection : public FEBioFileSection
{
public:
	FEBioLoadsSection(FEFEBioImport* pim) : FEBioFileSection(pim){}
	void Parse(XMLTag& tag);

protected:
	void ParseBCForce    (XMLTag& tag);
	void ParseBodyForce  (XMLTag& tag);
	void ParseHeatSource (XMLTag& tag);
	void ParseSurfaceLoad(XMLTag& tag);
};
