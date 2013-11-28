#pragma once
#include "FEBioBoundarySection.h"

//-----------------------------------------------------------------------------
// Loads Section (new in version 1.2)
class FEBioLoadsSection : public FEBioFileSection
{
public:
	FEBioLoadsSection(FEFEBioImport* pim) : FEBioFileSection(pim){}
	void Parse(XMLTag& tag);

public:
	void ParseBCForce             (XMLTag& tag);
	void ParseBCPressure          (XMLTag& tag);
	void ParseBCTraction          (XMLTag& tag);
	void ParseBCPoroNormalTraction(XMLTag& tag);
	void ParseBCFluidFlux         (XMLTag& tag);
	void ParseBCSoluteFlux        (XMLTag &tag);
	void ParseBCHeatFlux          (XMLTag& tag);
	void ParseBCConvectiveHeatFlux(XMLTag& tag);

protected:
	void ParseBodyForce (XMLTag& tag);
	void ParseHeatSource(XMLTag& tag);
};
