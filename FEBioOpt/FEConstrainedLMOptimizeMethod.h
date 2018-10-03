#pragma once
#include "FEOptimizeMethod.h"
#include <FECore/matrix.h>

//----------------------------------------------------------------------------
//! Optimization method using contrained Levenberg-Marquardt method
#ifdef HAVE_LEVMAR
class FEConstrainedLMOptimizeMethod : public FEOptimizeMethod
{
public:
	FEConstrainedLMOptimizeMethod();
	bool Solve(FEOptimizeData* pOpt, vector<double>& amin, vector<double>& ymin, double* minObj) override;

	FEOptimizeData* GetOptimizeData() { return m_pOpt; }

protected:
	FEOptimizeData* m_pOpt;

	void ObjFun(vector<double>& x, vector<double>& a, vector<double>& y, matrix& dyda);

	static FEConstrainedLMOptimizeMethod* m_pThis;
	static void objfun(vector<double>& x, vector<double>& a, vector<double>& y, matrix& dyda) { return m_pThis->ObjFun(x, a, y, dyda); }

public:
	double	m_tau;		// scale factor for mu
	double	m_objtol;	// objective tolerance
	double	m_fdiff;	// forward difference step size
	int		m_nmax;		// maximum number of iterations
    int     m_loglevel; // log file output level

public:
	vector<double>	m_yopt;	// optimal y-values

	DECLARE_FECORE_CLASS();
};
#endif
