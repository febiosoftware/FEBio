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
#include <FECore/LinearSolver.h>
#include "BlockMatrix.h"

class FGMRESSolver;

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

	// Use CG for step 2 or not
	void UseConjugateGradient(bool b);

	// set the CG convergence parameters
	void SetCGParameters(int maxiter, double tolerance, bool doResidualStoppingTest);

	// set the GMRES convergence parameters
	void SetGMRESParameters(int maxiter, double tolerance, bool doResidualStoppingTest, int precondition);

	// Do Jacobi preconditioner
	void DoJacobiPreconditioner(bool b);

	// set the schur preconditioner option
	void SetSchurPreconditioner(int n);

public:
	// allocate storage
	bool PreProcess() override;

	//! Factor the matrix (for iterative solvers, this can be used for creating pre-conditioner)
	bool Factor() override;

	//! Calculate the solution of RHS b and store solution in x
	bool BackSolve(double* x, double* y) override;

	//! Return a sparse matrix compatible with this solver
	SparseMatrix* CreateSparseMatrix(Matrix_Type ntype) override;

	// set the sparse matrix
	bool SetSparseMatrix(SparseMatrix* A) override;

private:
	int cgsolve(SparseMatrix* K, LinearSolver* PC, vector<double>& x, vector<double>& b);
	int gmressolve(SparseMatrix* K, LinearSolver* PC, vector<double>& x, vector<double>& b);

private:
	BlockMatrix*	m_A;		//!< the block matrix
	FGMRESSolver*	m_Asolver;	//!< the solver for the A - block
	LinearSolver*	m_PS;		//!< Schur complement preconditioner

	std::vector<double>		Kd;
	std::vector<double>		Wm, Wc;
	std::vector<double>		Rm0, Rc0;
	std::vector<double>		Rm_n, Rc_n;
	std::vector<double>		yu, yp;
	std::vector<double>		yu_n, yp_n;

	std::vector< std::vector<double> >	Yu, Yp;
	std::vector<double>	au, ap;
	std::vector<double> du, dp;

	vector<double> RM;
	vector<double> RC;

	vector< vector<double> > Rmu;
	vector< vector<double> > Rmp;
	vector< vector<double> > Rcu;
	vector< vector<double> > Rcp;

	int		m_print_level;	//!< level of output (0 is no output)
	int		m_maxiter;		//!< max nr of BIPN iterations
	double	m_tol;			//!< BPIN convergence tolerance

	bool	m_use_cg;		//!< use CG for step 2, otherwise GMRES is used

	bool	m_do_jacobi;	//!< do Jacobi precondition

	int		m_precondition_schur;	//!< preconditioner the Schur solver (0 = none, 1 = diag(M), 2 = ICHOL(M))

	// CG data
	int		m_cg_maxiter;			//!< max CG iterations
	double	m_cg_tol;				//!< CG tolerance
	bool	m_cg_doResidualTest;	//!< do the residual stopping test
	int		m_cg_iters;				//!< iterations of CG solve

	vector<double>	cg_tmp;

	// GMRES data
	int		m_gmres_maxiter;		//!< max GMRES iterations
	double	m_gmres_tol;			//!< GMRES tolerance
	bool	m_gmres_doResidualTest;	//!< do the residual stopping test
	int		m_gmres_pc;				//!< preconditioner?
	int		m_gmres1_iters, m_gmres2_iters;

	vector<double> gmres_tmp;

	DECLARE_FECORE_CLASS();
};
