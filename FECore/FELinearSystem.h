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
#include "FEGlobalMatrix.h"
#include "matrix.h"
#include <vector>

class FESolver;

//-----------------------------------------------------------------------------
// Experimental class to see if all the assembly operations can be moved to a class
// and out of the solver class.
class FECORE_API FELinearSystem
{
public:
	// Constructor
	// Takes a FEGlobalMatrix class K that will store the actual stiffness matrix
	// and a vector F which contains the assembled contribution of the prescribed 
	// degrees of freedom. The F vector must be added to the "force" vector. The u 
	// vector contains the nodal values of the prescribed degrees of freedom.
	FELinearSystem(FESolver* solver, FEGlobalMatrix& K, std::vector<double>& F, std::vector<double>& u, bool bsymm);

	// virtual destructor
	virtual ~FELinearSystem();

	// get symmetry flag
	bool IsSymmetric() const;

	// Get the solver that is using this linear system
	FESolver* GetSolver();

public:
	// Assembly routine
	// This assembles the element stiffness matrix ke into the global matrix.
	// The contributions of prescribed degrees of freedom will be stored in m_F
	virtual void Assemble(const FEElementMatrix& ke);

	// This assembles a matrix to the RHS by pre-multiplying the matrix with the 
	// prescribed value array U and then adding it to F
	void AssembleRHS(std::vector<int>& lm, matrix& ke, std::vector<double>& U);

	// This assembles a vetor to the RHS
	void AssembleRHS(std::vector<int>& lm, std::vector<double>& fe);

protected:
	bool					m_bsymm;	//!< symmetry flag
	FESolver*				m_solver;
	FEGlobalMatrix&			m_K;	//!< The global stiffness matrix
	std::vector<double>&	m_F;	//!< Contributions from prescribed degrees of freedom
	std::vector<double>&	m_u;	//!< the array with prescribed values
};
