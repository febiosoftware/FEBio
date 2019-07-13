/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, The Trustees of Columbia University in 
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
#include "SparseMatrix.h"
#include "FECoreBase.h"

//-----------------------------------------------------------------------------
class CRSSparseMatrix;
class CompactMatrix;

//-----------------------------------------------------------------------------
// Base class for preconditioners for iterative linear solvers
class FECORE_API Preconditioner : public FECoreBase
{
	FECORE_SUPER_CLASS

public:
	Preconditioner(FEModel* fem);
	virtual ~Preconditioner();

	// return the sparse matrix
	SparseMatrix* GetSparseMatrix();

	// set the sparse matrix
	void SetSparseMatrix(SparseMatrix* A);

	// create a preconditioner for a sparse matrix
	virtual bool Create() = 0;
	
	// apply to vector P x = y
	virtual bool mult_vector(double* x, double* y) = 0;

private:
	SparseMatrix*	m_K;
};

//-----------------------------------------------------------------------------
class FECORE_API DiagonalPreconditioner : public Preconditioner
{
public:
	DiagonalPreconditioner(FEModel* fem);

	// take square root of diagonal entries
	void CalculateSquareRoot(bool b);

	// create a preconditioner for a sparse matrix
	bool Create() override;

	bool Create(double d);

	// apply to vector P x = y
	bool mult_vector(double* x, double* y) override;

private:
	vector<double>	m_D;

	bool	m_bsqr;		// Take square root
};
