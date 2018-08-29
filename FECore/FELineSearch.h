#pragma once
#include "DumpStream.h"

class FENewtonSolver;

class FELineSearch
{
public:
	FELineSearch(FENewtonSolver* pns);

	// Do a line search. ls is the initial search step
	double DoLineSearch(double ls = 1.0);

	// serialization
	void Serialize(DumpStream& ar);

public:
	double	m_LSmin;		//!< minimum line search step
	double	m_LStol;		//!< line search tolerance
	int		m_LSiter;		//!< max nr of line search iterations

private:
	FENewtonSolver*	m_pns;
};
