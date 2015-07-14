#pragma once
#include "FEBioImport.h"

//-----------------------------------------------------------------------------
//! Globals Section
class FEBioGlobalsSection : public FEBioFileSection
{
public:
	FEBioGlobalsSection(FEBioImport* pim) : FEBioFileSection(pim){}
	void Parse            (XMLTag& tag);

protected:
	void ParseConstants   (XMLTag& tag);
	void ParseGlobalData  (XMLTag& tag);
};
