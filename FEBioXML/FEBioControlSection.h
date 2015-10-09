#pragma once
#include "FEBioImport.h"

//-----------------------------------------------------------------------------
// Control Section
class FEBioControlSection : public FEBioFileSection
{
public:
	FEBioControlSection(FEBioImport* pim) : FEBioFileSection(pim) {}
	void Parse(XMLTag& tag);

protected:
	FESolver* BuildSolver(const char* sztype, FEModel& fem);
	bool ParseCommonParams(XMLTag& tag);
};
