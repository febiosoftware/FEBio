#pragma once
#include "FEBioImport.h"

//-----------------------------------------------------------------------------
//! Globals Section
class FEBioGlobalsSection : public FEBioFileSection
{
public:
	FEBioGlobalsSection(FEFEBioImport* pim) : FEBioFileSection(pim){}
	void Parse            (XMLTag& tag);

protected:
	void ParseBodyForce   (XMLTag& tag);	// only for versions < 1.2
	void ParseConstants   (XMLTag& tag);
	void ParseGSSoluteData(XMLTag& tag);
	void ParseMGData      (XMLTag &tag);
};
