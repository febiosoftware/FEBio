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
	virtual bool Solve(FEOptimizeData* pOpt) = 0;

public:
	Logfile::MODE	m_loglevel;		//!< log file output level
	int				m_print_level;	//!< level of detailed output
};
