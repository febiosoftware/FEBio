// LinSolver.cpp: implementation of the LinSolver class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "LinearSolver.h"

//-----------------------------------------------------------------------------
// The variable m_numthreads determines the number of threads to request for
// solvers that can take advantage of multiple processors.

int FECore::LinearSolver::m_numthreads = 1;

