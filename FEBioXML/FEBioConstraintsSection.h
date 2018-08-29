#pragma once
#include "FEBioImport.h"

class FEFacetSet;

//-----------------------------------------------------------------------------
// Constraints Section
// (Base class. Don't use this directly!)
class FEBioConstraintsSection : public FEFileSection
{
public:
	FEBioConstraintsSection(FEFileImport* pim) : FEFileSection(pim){}

protected:
	bool ParseSurfaceSection(XMLTag& tag, FESurface& s, int nfmt, bool bnodal);
};

//-----------------------------------------------------------------------------
// Constraints Section (format 1.x)
class FEBioConstraintsSection1x : public FEBioConstraintsSection
{
public:
	FEBioConstraintsSection1x(FEFileImport* pim) : FEBioConstraintsSection(pim){}
	void Parse(XMLTag& tag);

protected:
	void ParseRigidConstraint(XMLTag& tag);
};

//-----------------------------------------------------------------------------
// Constraints Section (format 2.0)
class FEBioConstraintsSection2 : public FEBioConstraintsSection
{
public:
	FEBioConstraintsSection2(FEFileImport* pim) : FEBioConstraintsSection(pim){}
	void Parse(XMLTag& tag);
	
protected:
	void ParseRigidConstraint20(XMLTag& tag);
};

//-----------------------------------------------------------------------------
// Constraints Section (format 2.5)
class FEBioConstraintsSection25 : public FEBioConstraintsSection
{
public:
	FEBioConstraintsSection25(FEFileImport* pim) : FEBioConstraintsSection(pim){}
	void Parse(XMLTag& tag);
};
