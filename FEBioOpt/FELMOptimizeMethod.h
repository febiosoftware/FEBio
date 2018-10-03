#pragma once
#include "FEOptimizeMethod.h"
#include <FECore/matrix.h>
#include <vector>
using namespace std;

//----------------------------------------------------------------------------
//! Optimization method using Levenberg-Marquardt method
class FELMOptimizeMethod : public FEOptimizeMethod
{
public:
	FELMOptimizeMethod();
	bool Solve(FEOptimizeData* pOpt, vector<double>& amin, vector<double>& ymin, double* minObj) override;

protected:
	FEOptimizeData* m_pOpt;

	void ObjFun(vector<double>& x, vector<double>& a, vector<double>& y, matrix& dyda);

	static FELMOptimizeMethod* m_pThis;
	static void objfun(vector<double>& x, vector<double>& a, vector<double>& y, matrix& dyda) { return m_pThis->ObjFun(x, a, y, dyda); }

public:
	double			m_objtol;	// objective tolerance
	double			m_fdiff;	// forward difference step size
	int				m_nmax;		// maximum number of iterations
	bool			m_bcov;		// flag to print covariant matrix

protected:
	vector<double>	m_yopt;	// optimal y-values

	DECLARE_FECORE_CLASS();
};
