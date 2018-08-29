#pragma once
#include "FEOptimizeMethod.h"

//----------------------------------------------------------------------------
//! Optimization method using Powell's method
class FEPowellOptimizeMethod : public FEOptimizeMethod
{
public:
	bool Solve(FEOptimizeData* pOpt, vector<double>& amin, vector<double>& ymin, double* minObj);

protected:
	double ObjFun(double* p);

	FEOptimizeData*	m_pOpt;
	
	static FEPowellOptimizeMethod*	m_pThis;
	static double objfun(double* p) { return (m_pThis)->ObjFun(p); }
};
