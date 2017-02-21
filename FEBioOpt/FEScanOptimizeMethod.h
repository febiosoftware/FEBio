#pragma once
#include "FEOptimizeMethod.h"

//----------------------------------------------------------------------------
//! Basic method that scans the parameter space for a minimum.
class FEScanOptimizeMethod : public FEOptimizeMethod
{
public:
	bool Solve(FEOptimizeData* pOpt);

protected:
	bool FESolve(vector<double>& x, vector<double>& a, vector<double>& y);

private:
	FEOptimizeData* m_pOpt;
	double	m_inc[32];

	DECLARE_PARAMETER_LIST();
};

