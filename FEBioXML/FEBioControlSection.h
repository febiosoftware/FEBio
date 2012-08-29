#pragma once
#include "FEBioImport.h"

//-----------------------------------------------------------------------------
// Control Section
class FEBioControlSection : public FEBioFileSection
{
public:
	FEBioControlSection(FEFEBioImport* pim) : FEBioFileSection(pim) {}
	void Parse(XMLTag& tag);

protected:
	FESolver* BuildSolver(int nmod, FEModel& fem);
	bool ParseCommonParams(XMLTag& tag);
};
