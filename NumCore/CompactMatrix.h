#pragma once
#include "FECore/SparseMatrix.h"
#include "CSRMatrix.h"

//=============================================================================
//! This class stores a sparse matrix in Harwell-Boeing format.

//! This is the base class for the symmetric and unsymmetric classes

class CompactMatrix : public SparseMatrix
{
public:
	//! constructor
	CompactMatrix( int offset );

	//! destructor
	virtual ~CompactMatrix();

	//! zero matrix elements
	void Zero() override;

	//! Clear 
	void Clear() override;

public:
	//! Pointer to matrix values
	double* Values  () override { return m_pd;   }

	//! Pointer to matrix indices
	int*    Indices() override { return m_pindices; }

	//! pointer to matrix row pointers
	int*    Pointers() override { return m_ppointers; }

	//! return the index offset (is 0 or 1)
	int     Offset() const override { return m_offset; }

public:
	//! Create the matrix
	void alloc(int nr, int nc, int nz, double* pv, int *pi, int* pp, bool bdel = true);

	//! is the matrix symmetric or not
	virtual bool isSymmetric() = 0;

	//! is this a row-based format or not
	virtual bool isRowBased() = 0;

protected:
	double*	m_pd;			//!< matrix values
	int*	m_pindices;		//!< indices
	int*	m_ppointers;	//!< pointers
	int		m_offset;		//!< adjust array indices for fortran arrays
	bool	m_bdel;			//!< delete data arrays in destructor
};
