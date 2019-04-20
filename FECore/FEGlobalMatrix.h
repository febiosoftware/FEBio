/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, Columbia University, and others.

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
#include "FESolver.h"
#include <vector>

//-----------------------------------------------------------------------------
class FEModel;
class FEMesh;
class FESurface;

//-----------------------------------------------------------------------------
class FEElementMatrix : public matrix
{
	FEElementMatrix(vector<int>& lmi, vector<int>& lmj) : matrix(lmi.size(), lmj.size())
	{
		m_lmi = lmi;
		m_lmj = lmj;
	};

public:
	std::vector<int>	m_lmi;
	std::vector<int>	m_lmj;
};

//-----------------------------------------------------------------------------
//! This class implements a global system matrix.

//! The global system matrix is usually created by the discretization of the FE 
//! equations into a linear system of equations. The structure of it depends greatly
//! on the element connectivity and usually results in a sparse matrix structure. 
//! Several sparse matrix structures are supported (Compact, Skyline, etc.) and to 
//! simplify the creation of the specific matrix structure, the FEGlobalMatrix offers
//! functionality to create the global matrix structure without the need to know 
//! what particular sparse matrix format is used by the linear solver.

//! \todo I think the SparseMatrixProfile can handle all of the build functions.

class FECORE_API FEGlobalMatrix
{
protected:
	enum { MAX_LM_SIZE = 64000 };

public:
	//! constructor
	FEGlobalMatrix(SparseMatrix* pK, bool del = true);

	//! destructor
	virtual ~FEGlobalMatrix();

	//! construct the stiffness matrix from a FEM object
	bool Create(FEModel* pfem, int neq, bool breset);

	//! construct the stiffness matrix from a mesh
	bool Create(FEMesh& mesh, int neq);

	//! construct the stiffness matrix from a mesh
	bool Create(FEMesh& mesh, int nstart, int nend);

	//! construct a stiffness matrix from a surface
	//! The boundary array is a list of equation numbers.
	bool Create(FESurface& surf, const std::vector<int>& equationIDs);

	//! clears the sparse matrix that stores the stiffness matrix
	void Clear();

	//! assemble an element stiffness matrix into the global stiffness matrix
	void Assemble(matrix& ke, vector<int>& lm) { m_pA->Assemble(ke, lm); }

	//! more general assembly routine
	void Assemble(matrix& ke, vector<int>& lmi, vector<int>& lmj) { m_pA->Assemble(ke, lmi, lmj); }

	//! Even more general assembly routine
	void Assemble(FEElementMatrix& ke) { m_pA->Assemble(ke, ke.m_lmi, ke.m_lmj); }

	//! return the nonzeroes in the sparse matrix
	int NonZeroes() { return m_pA->NonZeroes(); }

	//! return the number of rows
	int Rows() { return m_pA->Rows(); }

	//! converts a FEGlobalMatrix to a SparseMatrix
	operator SparseMatrix* () { return m_pA; }

	//! converts a FEGlobalMatrix to a SparseMatrix
	operator SparseMatrix& () { return *m_pA;}

	//! return a pointer to the sparse matrix
	SparseMatrix* GetSparseMatrixPtr() { return m_pA; }

	//! zero the sparse matrix
	void Zero() { m_pA->Zero(); }

public:
	void build_begin(int neq);
	void build_add(std::vector<int>& lm);
	void build_end();
	void build_flush();

protected:
	SparseMatrix*	m_pA;	//!< the actual global stiffness matrix
	bool			m_delA;	//!< delete A in destructor

	// The following data structures are used to incrementally
	// build the profile of the sparse matrix

	SparseMatrixProfile*	m_pMP;		//!< profile of sparse matrix
	SparseMatrixProfile		m_MPs;		//!< the "static" part of the matrix profile
	vector< vector<int> >	m_LM;		//!< used for building the stiffness matrix
	int	m_nlm;				//!< nr of elements in m_LM array
};
