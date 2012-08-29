#pragma once
#include "FEBioBoundarySection.h"

//-----------------------------------------------------------------------------
// Loads Section (new in version 1.2)
class FEBioLoadsSection : public FEBioBoundarySection
{
public:
	FEBioLoadsSection(FEFEBioImport* pim) : FEBioBoundarySection(pim){}
	void Parse(XMLTag& tag);

protected:
	void ParseBodyForce(XMLTag& tag);
};
