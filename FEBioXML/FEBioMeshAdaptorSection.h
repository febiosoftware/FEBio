#pragma once
#include "FEBioImport.h"

class FEBioMeshAdaptorSection : public FEFileSection
{
public:
	FEBioMeshAdaptorSection(FEFileImport* pim) : FEFileSection(pim) {}

	void Parse(XMLTag& tag);

protected:
	void ParseMeshAdaptor(XMLTag& tag);
};
