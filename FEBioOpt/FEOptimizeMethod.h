/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/



#pragma once
#include <FECore/FECoreClass.h>

//-----------------------------------------------------------------------------
class FEOptimizeData;

//-----------------------------------------------------------------------------
enum {
	PRINT_ITERATIONS,
	PRINT_VERBOSE
};

//-----------------------------------------------------------------------------
enum LogLevel {
	LOG_NEVER,
	LOG_FILE,
	LOG_SCREEN,
	LOG_FILE_AND_SCREEN
};

//-----------------------------------------------------------------------------
//! Base class for optimization algorithms.
//! Derived class implement specific optimization algorithms.
class FEOptimizeMethod : public FECoreClass
{
	FECORE_BASE_CLASS(FEOptimizeMethod)

public:
	FEOptimizeMethod(FEModel* fem);

	// Implement this function for solve an optimization problem
	// should return the optimal values for the input parameters in a, the optimal
	// values of the measurement vector in ymin and
	// the corresponding objective value in obj.
	// If this function returns false, something went wrong
	virtual bool Solve(FEOptimizeData* pOpt, vector<double>& amin, vector<double>& ymin, double* minObj) = 0;

public:
	int		m_loglevel;		//!< log file output level
	int		m_print_level;	//!< level of detailed output
};
