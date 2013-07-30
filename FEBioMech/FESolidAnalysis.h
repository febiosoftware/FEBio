#pragma once
#include "FECore/FEAnalysis.h"
using namespace FECore;

//-----------------------------------------------------------------------------
//! This class describes a analysis step for structural mechanics problems
class FESolidAnalysis : public FEAnalysis
{
public:
	//! constructor
	FESolidAnalysis(FEModel& fem);

	//! initialization
	bool Init();
};

//-----------------------------------------------------------------------------
//! This class describes a analysis step for structural mechanics problems
class FEExplicitSolidAnalysis : public FEAnalysis
{
public:
	//! constructor
	FEExplicitSolidAnalysis(FEModel& fem) : FEAnalysis(fem, FE_EXPLICIT_SOLID){}

	//! initialization
	bool Init();
};

//-----------------------------------------------------------------------------
//! Analysis class for linear elastic problems
class FELinearSolidAnalysis : public FEAnalysis
{
public:
	//! constructor
	FELinearSolidAnalysis(FEModel& fem) : FEAnalysis(fem, FE_LINEAR_SOLID) {}

	//! initialization
	bool Init();
};
