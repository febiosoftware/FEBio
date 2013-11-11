#pragma once
#include "FECore/FEAnalysis.h"
using namespace FECore;

//-----------------------------------------------------------------------------
//! This class describes a analysis step for structural mechanics problems
class FESolidAnalysis : public FEAnalysis
{
public:
	//! constructor
	FESolidAnalysis(FEModel* pfem);

	//! initialization
	bool Init();
};

//-----------------------------------------------------------------------------
//! This class describes a analysis step for structural mechanics problems
class FEExplicitSolidAnalysis : public FEAnalysis
{
public:
	//! constructor
	FEExplicitSolidAnalysis(FEModel* pfem) : FEAnalysis(pfem, FE_EXPLICIT_SOLID){}

	//! initialization
	bool Init();
};

//-----------------------------------------------------------------------------
//! Analysis class for linear elastic problems
class FELinearSolidAnalysis : public FEAnalysis
{
public:
	//! constructor
	FELinearSolidAnalysis(FEModel* pfem) : FEAnalysis(pfem, FE_LINEAR_SOLID) {}

	//! initialization
	bool Init();
};
