#pragma once
#include "FECore/FEAnalysis.h"
using namespace FECore;

//-----------------------------------------------------------------------------
//! Analysis class for biphasic-solute problems
class FEBiphasicSoluteAnalysis : public FEAnalysis
{
public:
	FEBiphasicSoluteAnalysis(FEModel* pfem) : FEAnalysis(pfem, FE_POROSOLUTE) {}

	bool Activate();

protected:
	void InitNodes();
};
