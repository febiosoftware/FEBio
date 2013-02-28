#pragma once
#include "FEOptimizer.h"

//----------------------------------------------------------------------------
//! Optimization method using Powell's method
class FEPowellOptimizeMethod : public FEOptimizeMethod
{
public:
	bool Solve(FEOptimizeData* pOpt);

protected:
	double ObjFun(double* p);

	FEOptimizeData*	m_pOpt;
	
	static FEPowellOptimizeMethod*	m_pThis;
	static double objfun(double* p) { return (m_pThis)->ObjFun(p); }
};
