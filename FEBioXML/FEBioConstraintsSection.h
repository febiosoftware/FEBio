#pragma once
#include "FEBioImport.h"

class FEFacetSet;

//-----------------------------------------------------------------------------
// Constraints Section
class FEBioConstraintsSection : public FEBioFileSection
{
public:
	FEBioConstraintsSection(FEBioImport* pim) : FEBioFileSection(pim){}
	void Parse(XMLTag& tag);

	void Parse20(XMLTag& tag);
	void Parse25(XMLTag& tag);

protected:
	void ParseRigidConstraint(XMLTag& tag);
	void ParseRigidConstraint20(XMLTag& tag);
	bool ParseSurfaceSection (XMLTag& tag, FESurface& s, int nfmt, bool bnodal);
	bool BuildSurface(FESurface& s, FEFacetSet& f, bool bnodal);
};
