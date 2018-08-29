#pragma once
#include "FEBioImport.h"

class FEBioRigidSection : public FEBioFileSection
{
public:
	FEBioRigidSection(FEBioImport* pim) : FEBioFileSection(pim){}
	void Parse(XMLTag& tag);
};
