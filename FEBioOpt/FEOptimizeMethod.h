#pragma once
#include <FECore/FEParameterList.h>
#include <FECore/log.h>

//-----------------------------------------------------------------------------
class FEOptimizeData;

//-----------------------------------------------------------------------------
enum {
	PRINT_ITERATIONS,
	PRINT_VERBOSE
};


//-----------------------------------------------------------------------------
//! Base class for optimization algorithms.
//! Derived class implement specific optimization algorithms.
class FEOptimizeMethod : public FEParamContainer
{
public:
	// Implement this function for solve an optimization problem
	// should return the optimal values for the input parameters in a, the optimal
	// values of the measurement vector in ymin and
	// the corresponding objective value in obj.
	// If this function returns false, something went wrong
	virtual bool Solve(FEOptimizeData* pOpt, vector<double>& amin, vector<double>& ymin, double* minObj) = 0;

public:
	Logfile::MODE	m_loglevel;		//!< log file output level
	int				m_print_level;	//!< level of detailed output
};
