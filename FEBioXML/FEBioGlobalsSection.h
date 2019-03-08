#pragma once
#include "FileImport.h"

//-----------------------------------------------------------------------------
//! Globals Section
class FEBIOXML_API FEBioGlobalsSection : public FEFileSection
{
public:
	FEBioGlobalsSection(FEFileImport* pim) : FEFileSection(pim){}
	void Parse(XMLTag& tag);

protected:
	void ParseConstants   (XMLTag& tag);
	void ParseGlobalData  (XMLTag& tag);
};
