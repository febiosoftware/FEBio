#pragma once
#include "FEBioImport.h"

//-----------------------------------------------------------------------------
// LoadData Section
class FEBioLoadDataSection : public FEFileSection
{
public:
	FEBioLoadDataSection(FEFileImport* pim);
	void Parse(XMLTag& tag);

	// Set the redefine curves flag.
	// When this flag is set, curves can be redefined by using an existing ID
	void SetRedefineCurvesFlag(bool b) { m_redefineCurves = b; }

protected:
	bool	m_redefineCurves;	// flag to allow redefining curves
};

//-----------------------------------------------------------------------------
// LoadData Section
class FEBioLoadDataSection3 : public FEFileSection
{
public:
	FEBioLoadDataSection3(FEFileImport* pim);
	void Parse(XMLTag& tag);
};
