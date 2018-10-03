#pragma once
#include "FEOptimizeMethod.h"

//----------------------------------------------------------------------------
//! Basic method that scans the parameter space for a minimum.
class FEScanOptimizeMethod : public FEOptimizeMethod
{
public:
	// this implements the solution algorithm.
	// returns the optimal parameter values in amin
	// returns the optimal measurement vector in ymin
	// returns the optimal objective function value in minObj
	bool Solve(FEOptimizeData* pOpt, vector<double>& amin, vector<double>& ymin, double* minObj) override;

	DECLARE_FECORE_CLASS();
};
