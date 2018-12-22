#pragma once
#include <FECore/LinearSolver.h>
#include <FECore/SparseMatrix.h>
#include <FECore/CompactUnSymmMatrix.h>

//-----------------------------------------------------------------------------
// This class implements the bi-partitioned iterative solver, by:
// Esmaily-Moghadam, Bazilevs, Marsden, Comput. Methods Appl. Mech. Engrg. 286(2015) 40-62
//
class BIPNSolver : public LinearSolver
{
public:
	// constructor
	BIPNSolver(FEModel* fem);

	// set the output level
	void SetPrintLevel(int n) override;

	// set the max nr of BIPN iterations
	void SetMaxIterations(int n);

	// Set the BIPN convergence tolerance
	void SetTolerance(double eps);

	// Set split row/column
	void SetPartition(int n) override;

	// Use CG for step 2 or not
	void UseConjugateGradient(bool b);

	// set the CG convergence parameters
	void SetCGParameters(int maxiter, double tolerance, bool doResidualStoppingTest);

	// set the GMRES convergence parameters
	void SetGMRESParameters(int maxiter, double tolerance, bool doResidualStoppingTest, bool precondition);

public:
	// allocate storage
	bool PreProcess() override;

	//! Factor the matrix (for iterative solvers, this can be used for creating pre-conditioner)
	bool Factor() override;

	//! Calculate the solution of RHS b and store solution in x
	bool BackSolve(double* x, double* y) override;

	//! Return a sparse matrix compatible with this solver
	SparseMatrix* CreateSparseMatrix(Matrix_Type ntype) override;

private:
	bool step2_cgsolve(vector<double>& x, vector<double>& b);
	bool step2_gmressolve(vector<double>& x, vector<double>& b);
	bool gmressolve(vector<double>& x, vector<double>& b);

private:
	CRSSparseMatrix*		m_A;
	std::vector<double>		m_W;

	// blocks
	CSRMatrix	K, G, D, L;

	std::vector<double>		Wm, Wc;
	std::vector<double>		Rm, Rc;
	std::vector<double>		Rm_n, Rc_n;
	std::vector<double>		yu, yp;
	std::vector<double>		yu_n, yp_n;

	std::vector< std::vector<double> >	Yu, Yp;
	std::vector<double>	au, ap;
	std::vector<double> du, dp;

	vector< vector<double> > RM;
	vector< vector<double> > RC;

	vector< vector<double> > Rmu;
	vector< vector<double> > Rmp;
	vector< vector<double> > Rcu;
	vector< vector<double> > Rcp;

	int		m_print_level;	//!< level of output (0 is no output)
	int		m_split;		//!< set the split row index
	int		m_maxiter;		//!< max nr of BIPN iterations
	double	m_tol;			//!< BPIN convergence tolerance

	bool	m_use_cg;		//!< use CG for step 2, otherwise GMRES is used

	// CG data
	int		m_cg_maxiter;			//!< max CG iterations
	double	m_cg_tol;				//!< CG tolerance
	bool	m_cg_doResidualTest;	//!< do the residual stopping test

	vector<double>	cg_tmp;

	// GMRES data
	int		m_gmres_maxiter;		//!< max GMRES iterations
	double	m_gmres_tol;			//!< GMRES tolerance
	bool	m_gmres_doResidualTest;	//!< do the residual stopping test
	bool	m_gmres_ilu0;			//!< Use ILU0 preconditioner?

	vector<double> gmres_tmp;
};
