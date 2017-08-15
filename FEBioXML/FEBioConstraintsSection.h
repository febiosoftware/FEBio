#pragma once
#include "FEBioImport.h"

class FEFacetSet;

//-----------------------------------------------------------------------------
// Constraints Section
// (Base class. Don't use this directly!)
class FEBioConstraintsSection_ : public FEFileSection
{
public:
	FEBioConstraintsSection_(FEFileImport* pim) : FEFileSection(pim){}

protected:
	void ParseRigidConstraint(XMLTag& tag);
	void ParseRigidConstraint20(XMLTag& tag);
	bool ParseSurfaceSection(XMLTag& tag, FESurface& s, int nfmt, bool bnodal);
	bool BuildSurface(FESurface& s, FEFacetSet& f, bool bnodal);
};

//-----------------------------------------------------------------------------
// Constraints Section (format 1.x)
class FEBioConstraintsSection1x : public FEBioConstraintsSection_
{
public:
	FEBioConstraintsSection1x(FEFileImport* pim) : FEBioConstraintsSection_(pim){}
	void Parse(XMLTag& tag);
};

//-----------------------------------------------------------------------------
// Constraints Section (format 2.0)
class FEBioConstraintsSection2 : public FEBioConstraintsSection_
{
public:
	FEBioConstraintsSection2(FEFileImport* pim) : FEBioConstraintsSection_(pim){}
	void Parse(XMLTag& tag);
};

//-----------------------------------------------------------------------------
// Constraints Section (format 2.5)
class FEBioConstraintsSection25 : public FEBioConstraintsSection_
{
public:
	FEBioConstraintsSection25(FEFileImport* pim) : FEBioConstraintsSection_(pim){}
	void Parse(XMLTag& tag);
};
