#pragma once
#include "FEBioImport.h"

//-----------------------------------------------------------------------------
// Constraints Section
class FEBioConstraintsSection : public FEBioFileSection
{
public:
	FEBioConstraintsSection(FEFEBioImport* pim) : FEBioFileSection(pim){}
	void Parse(XMLTag& tag);

protected:
	void ParseRigidConstraint(XMLTag& tag);
	void ParsePointConstraint(XMLTag& tag);
};
