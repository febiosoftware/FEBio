#pragma once
#include "FECore/FEAnalysis.h"
using namespace FECore;

//-----------------------------------------------------------------------------
//! Analysis class for biphasic problems
class FEBiphasicAnalysis : public FEAnalysis
{
public:
	FEBiphasicAnalysis(FEModel* pfem) : FEAnalysis(pfem, FE_BIPHASIC) {}

	bool Init();

protected:
	void InitNodes();
};

//-----------------------------------------------------------------------------
//! Analysis class for biphasic-solute problems
class FEBiphasicSoluteAnalysis : public FEAnalysis
{
public:
	FEBiphasicSoluteAnalysis(FEModel* pfem) : FEAnalysis(pfem, FE_POROSOLUTE) {}

	bool Init();

protected:
	void InitNodes();
};
