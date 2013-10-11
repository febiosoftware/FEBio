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
	bool ParseSurfaceSection (XMLTag& tag, FESurface& s, int nfmt, bool bnodal);
};
