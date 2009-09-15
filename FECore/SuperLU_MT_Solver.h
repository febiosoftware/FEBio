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
	bool PreProcess(SparseMatrix& K);
	bool Factor(SparseMatrix& K);
	bool Solve(SparseMatrix& K, vector<double>& x, vector<double>& b);
	bool Solve(SparseMatrix& K, matrix& x, matrix& b) { return false; }
	void Destroy(SparseMatrix& K);

	SparseMatrix* GetMatrix(int ntype) { return new CompactUnSymmMatrix(); }

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
