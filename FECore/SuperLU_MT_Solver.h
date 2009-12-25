#pragma once
#include "LinearSolver.h"

//-----------------------------------------------------------------------------
//! This class implements a wrapper class for the SuperLU_MT solver

#ifdef SUPERLU_MT
	#include "pdsp_defs.h"
#endif

class SuperLU_MT_Solver : public LinearSolver
{
public:
	bool PreProcess();
	bool Factor();
	bool Solve(vector<double>& x, vector<double>& b);
	void Destroy();

	SparseMatrix* CreateSparseMatrix(int ntype) { return (m_pA = new CompactUnSymmMatrix()); }

	SuperLU_MT_Solver();

#ifdef SUPERLU_MT

protected:

	bool m_balloc;
	bool m_bfact;

	SuperMatrix m_A, m_L, m_U, m_B, m_X;
	vector<int>	m_perm_c;
	vector<int>	m_perm_r;
	vector<int>	etree;

    superlumt_options_t		m_ops;
	superlu_memusage_t		m_mem;

	double	rpg, rcond;
	double	ferr, berr;
	int		info;
	equed_t	equed;

#endif // SUPERLU_MT
};
