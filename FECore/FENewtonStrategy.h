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
#include <vector>
#include "fecore_api.h"
#include "fecore_enum.h"
#include "FECoreBase.h"

//-----------------------------------------------------------------------------
class FENewtonSolver;
class SparseMatrix;
class LinearSolver;

//-----------------------------------------------------------------------------
//! A Base class for newton-type solution strategies
class FECORE_API FENewtonStrategy : public FECoreBase
{
	FECORE_SUPER_CLASS(FENEWTONSTRATEGY_ID)
	FECORE_BASE_CLASS(FENewtonStrategy)

public:
	FENewtonStrategy(FEModel* fem);
	virtual ~FENewtonStrategy();

	void SetNewtonSolver(FENewtonSolver* solver);

	void Serialize(DumpStream& ar) override;

	//! reset data for new run
	virtual void Reset();

public:
	//! initialize the linear system
	virtual SparseMatrix* CreateSparseMatrix(Matrix_Type mtype);

	//! Presolve update
	virtual void PreSolveUpdate() {}

	//! perform a Newton udpate (returning false will force a matrix reformations)
	virtual bool Update(double s, std::vector<double>& ui, std::vector<double>& R0, std::vector<double>& R1) = 0;

	//! solve the equations
	virtual void SolveEquations(std::vector<double>& x, std::vector<double>& b) = 0;

	//! reform the stiffness matrix
	virtual bool ReformStiffness();

	//! calculate the residual
	virtual bool Residual(std::vector<double>& R, bool binit);

public:
	int		m_maxups;		//!< max nr of QN iters permitted between stiffness reformations
	int		m_max_buf_size;	//!< max buffer size for update vector storage
	bool	m_cycle_buffer;	//!< recycle the buffer when updates is larger than buffer size
	double	m_cmax;			//!< maximum value for the condition number
	int		m_nups;			//!< nr of stiffness updates

protected:
	FENewtonSolver*	m_pns;
};
