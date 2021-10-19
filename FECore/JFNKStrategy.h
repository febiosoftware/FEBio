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
#include "FENewtonStrategy.h"
#include "SparseMatrix.h"

class JFNKMatrix;

//-----------------------------------------------------------------------------
// Implements a Jacobian-Free Newton-Krylov strategy
class JFNKStrategy : public FENewtonStrategy
{
public:
	JFNKStrategy(FEModel* fem);

	//! New initialization method
	bool Init() override;

	//! initialize the linear system
	SparseMatrix* CreateSparseMatrix(Matrix_Type mtype) override;

	//! perform a BFGS udpate
	bool Update(double s, vector<double>& ui, vector<double>& R0, vector<double>& R1) override;

	//! solve the equations
	void SolveEquations(vector<double>& x, vector<double>& b) override;

	//! Overide reform stiffness because we don't want to do any reformations
	bool ReformStiffness() override;

	//! override so we can store a copy of the residual before we add Fd
	bool Residual(std::vector<double>& R, bool binit) override;

private:
	double				m_jfnk_eps;			//!< JFNK epsilon

public:
	// keep a pointer to the linear solver
	LinearSolver*	m_plinsolve;		//!< pointer to linear solver
	int				m_neq;				//!< number of equations
	bool			m_bprecondition;	//!< the solver requires preconditioning

	JFNKMatrix*		m_A;

	DECLARE_FECORE_CLASS();
};
