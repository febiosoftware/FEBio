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

#include "SparseMatrix.h"
#include "FESolver.h"
#include <vector>

//-----------------------------------------------------------------------------
class FEModel;
class FEMesh;
class FESurface;
class FEElement;

//-----------------------------------------------------------------------------
//! This class represents an element matrix, i.e. a matrix of values and the row and
//! column indices of the corresponding matrix elements in the global matrix. 
class FECORE_API FEElementMatrix : public matrix
{
public:
	// default constructor
	FEElementMatrix(){}
	FEElementMatrix(int nr, int nc) : matrix(nr, nc) {}
	FEElementMatrix(const FEElement& el);

	// constructor for symmetric matrices
	FEElementMatrix(const FEElement& el, const vector<int>& lmi);

	// constructor
	FEElementMatrix(const FEElement& el, vector<int>& lmi, vector<int>& lmj);

	// copy constructor
	FEElementMatrix(const FEElementMatrix& ke);
	FEElementMatrix(const FEElementMatrix& ke, double scale);

	// assignment operator
	void operator = (const matrix& ke);

	// row indices
	std::vector<int>& RowIndices() { return m_lmi; }
	const std::vector<int>& RowIndices() const { return m_lmi; }

	// column indices
	std::vector<int>& ColumnsIndices() { return m_lmj; }
	const std::vector<int>& ColumnsIndices() const { return m_lmj; }

	// set the row and columnd indices (assuming they are the same)
	void SetIndices(const std::vector<int>& lm) { m_lmi = m_lmj = lm; }

	// set the row and columnd indices
	void SetIndices(const std::vector<int>& lmr, const std::vector<int>& lmc) { m_lmi = lmr; m_lmj = lmc; }

	// Set the node indices
	void SetNodes(const std::vector<int>& en) { m_node = en; }

	// get the nodes
	const std::vector<int>& Nodes() const { return m_node; }

private:
	std::vector<int>	m_node;	//!< node indices
	std::vector<int>	m_lmi;	//!< row indices
	std::vector<int>	m_lmj;	//!< column indices
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
	bool Create(const FESurface& surf, const std::vector<int>& equationIDs);

	//! clears the sparse matrix that stores the stiffness matrix
	void Clear();

	//! Assembly routine
	virtual void Assemble(const FEElementMatrix& ke);

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

	//! get the sparse matrix profile
	SparseMatrixProfile* GetSparseMatrixProfile() { return m_pMP; }

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
