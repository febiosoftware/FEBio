#pragma once
#include "FEOptimizer.h"
#include <vector>
using namespace std;

//----------------------------------------------------------------------------
//! Optimization method using Levenberg-Marquardt method
class FELMOptimizeMethod : public FEOptimizeMethod
{
public:
	FELMOptimizeMethod();
	bool Solve(FEOptimizeData* pOpt);

protected:
	FEOptimizeData* m_pOpt;

	void ObjFun(vector<double>& x, vector<double>& a, vector<double>& y, matrix& dyda);

	bool FESolve(vector<double>& x, vector<double>& a, vector<double>& y);

	static FELMOptimizeMethod* m_pThis;
	static void objfun(vector<double>& x, vector<double>& a, vector<double>& y, matrix& dyda) { return m_pThis->ObjFun(x, a, y, dyda); }

public:
	double	m_objtol;	// objective tolerance
	double	m_fdiff;	// forward difference step size

protected:
	vector<double>	m_yopt;	// optimal y-values
	vector<double>	m_y0;	// initial (target) y-values
};
