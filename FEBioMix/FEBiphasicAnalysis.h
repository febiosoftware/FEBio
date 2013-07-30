#pragma once
#include "FECore/FEAnalysis.h"
using namespace FECore;

//-----------------------------------------------------------------------------
//! Analysis class for biphasic problems
class FEBiphasicAnalysis : public FEAnalysis
{
public:
	FEBiphasicAnalysis(FEModel& fem) : FEAnalysis(fem, FE_BIPHASIC) {}

	bool Init();
};

//-----------------------------------------------------------------------------
//! Analysis class for biphasic-solute problems
class FEBiphasicSoluteAnalysis : public FEAnalysis
{
public:
	FEBiphasicSoluteAnalysis(FEModel& fem) : FEAnalysis(fem, FE_POROSOLUTE) {}

	bool Init();
};
