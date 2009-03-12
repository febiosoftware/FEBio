#pragma once

#include "SparseMatrix.h"
#include "vector.h"
#include "matrix.h"

//-----------------------------------------------------------------------------
//! base class for the linear solver classes

//! This class defines several virtual functions that need to be overriden
//! in the derived class

class LinearSolver
{
public:
	LinearSolver() { m_bvalid = false; }
	virtual ~LinearSolver() { Destroy(); }

	virtual bool PreProcess(SparseMatrix& K) { m_bvalid = true; return true; }
	virtual bool Factor(SparseMatrix& K) = 0;
	virtual bool Solve(SparseMatrix& K, vector<double>& x, vector<double>& b) = 0;
	virtual bool Solve(SparseMatrix& K, matrix& x, matrix& b) = 0;
	virtual void Destroy() { m_bvalid = false; };

	virtual SparseMatrix* GetMatrix() = 0;

protected:
	bool	m_bvalid;	// flag indication wether a valid matrix structure is ready
};

//-----------------------------------------------------------------------------
//! LU decomposition solver

//! This solver performs an LU decomposition and uses a backsolving algorithm
//! to solve the equations.
//! This solver uses the FullMatrix class and therefore is not the preferred
//! solver. It should only be used for small problems and only when the other
//! solvers are not adequate.

class LUSolver : public LinearSolver
{
	virtual bool PreProcess(SparseMatrix& K);
	virtual bool Factor(SparseMatrix& K);
	virtual bool Solve(SparseMatrix& K, vector<double>& x, vector<double>& b);
	virtual bool Solve(SparseMatrix& K, matrix& x, matrix& b);
	virtual void Destroy();

	virtual SparseMatrix* GetMatrix() { return new FullMatrix(); }

protected:
	vector<int>	indx;
};


//-----------------------------------------------------------------------------
//! Implements a linear solver that uses a skyline format

class SkylineSolver : public LinearSolver
{
public:
	virtual bool PreProcess(SparseMatrix& K);
	virtual bool Factor(SparseMatrix& K);
	virtual bool Solve(SparseMatrix& K, vector<double>& x, vector<double>& b);
	virtual bool Solve(SparseMatrix& K, matrix& x, matrix& b);
	virtual void Destroy();

	virtual SparseMatrix* GetMatrix() { return new SkylineMatrix(); }

};

//-----------------------------------------------------------------------------
//! Implements a linear solver that uses a compact column storage format.

//! This solver can only be used on systems where it is available,
//! such as SGI IRIX, SGI ALTIX, ...

#ifdef PSLDLT
extern "C" {
	void PSLDLT_Ordering(int token, int method);
	void PSLDLT_Preprocess(int, int, int*, int*, int*, double*);
	void PSLDLT_Factor(int, int, int*, int*, double*);
	void PSLDLT_Solve(int, double*, double*);
	void PSLDLT_Destroy(int token);
}
#endif // PSLDLT

class PSLDLTSolver : public LinearSolver
{
public:
	virtual bool PreProcess(SparseMatrix& K);
	virtual bool Factor(SparseMatrix& K);
	virtual bool Solve(SparseMatrix& K, vector<double>& x, vector<double>& b);
	virtual bool Solve(SparseMatrix& K, matrix& x, matrix& b);
	virtual void Destroy();

	virtual SparseMatrix* GetMatrix() { return new CompactMatrix(); }
};

//-----------------------------------------------------------------------------
//! Implements a wrapper class for the SuperLU library

//! This solver can only be used on systems where it is available.
//! This solver also uses some of the BLAS routines so this package also needs
//! to be available on the system. Although SuperLU comes with a stripped down
//! version of BLAS.

#ifdef SUPERLU
		#include "slu_ddefs.h"
#endif

class SuperLUSolver : public LinearSolver
{
public:
	virtual bool PreProcess(SparseMatrix& K);
	virtual bool Factor(SparseMatrix& K);
	virtual bool Solve(SparseMatrix& K, vector<double>& x, vector<double>& b);
	virtual bool Solve(SparseMatrix& K, matrix& x, matrix& b);
	virtual void Destroy();

	virtual SparseMatrix* GetMatrix() { return new CompactUnSymmMatrix(); }

	SuperLUSolver() { m_balloc = false; m_bfact = false; m_bcond = false; }

#ifdef SUPERLU
protected:
	double norm(SparseMatrix& K); // calculates the 1-norm of the matrix A
#endif

protected:

	bool m_balloc;
	bool m_bfact;
	bool m_bcond;	// calculate condition numbers


#ifdef SUPERLU

	SuperMatrix A, L, U, B, X;
	vector<int>	perm_c;
	vector<int>	perm_r;
	vector<int>	etree;

	superlu_options_t	options;
	SuperLUStat_t	stat;
	mem_usage_t	mem_usage;

	double	rpg, rcond;
	double	ferr, berr;
	int		info;
	char	equed[1];

#endif // SUPERLU
};

//-----------------------------------------------------------------------------
//! This class implements a wrapper class for the SuperLU_MT solver

#ifdef SUPERLU_MT
	#include "pdsp_defs.h"
#endif

class SuperLU_MT_Solver : public LinearSolver
{
public:
	virtual bool PreProcess(SparseMatrix& K);
	virtual bool Factor(SparseMatrix& K);
	virtual bool Solve(SparseMatrix& K, vector<double>& x, vector<double>& b);
	virtual bool Solve(SparseMatrix& K, matrix& x, matrix& b) { return false; }
	virtual void Destroy();

	virtual SparseMatrix* GetMatrix() { return new CompactUnSymmMatrix(); }

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

//-----------------------------------------------------------------------------
//! this class implements an iterative solver that implements the
//! conjugate gradient method

class ConjGradIterSolver : public LinearSolver
{
public:
	ConjGradIterSolver();

	virtual bool PreProcess(SparseMatrix& K);
	virtual bool Factor(SparseMatrix& K);
	virtual bool Solve(SparseMatrix& K, vector<double>& x, vector<double>& b);
	virtual bool Solve(SparseMatrix& K, matrix& x, matrix& b);
	virtual void Destroy();

	virtual SparseMatrix* GetMatrix() { return new CompactMatrix(); }

public:
	double	m_tol;		// convergence tolerance
	int		m_kmax;		// max iterations
	int		m_nprint;	// printing level

	vector<double>	m_P;	// preconditioning vector
};
