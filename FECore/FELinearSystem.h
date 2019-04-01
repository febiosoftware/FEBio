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
#include "FEGlobalMatrix.h"
#include "matrix.h"
#include <vector>
using namespace std;

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
	FELinearSystem(FEModel* fem, FEGlobalMatrix& K, vector<double>& F, vector<double>& u);

public:
	// Assembly routine
	// This assembles the element stiffness matrix ke into the global matrix.
	// The contributions of prescribed degrees of freedom will be stored in m_F
	void AssembleLHS(vector<int>& en, vector<int>& lm, matrix& ke);

	// This assembles a matrix to the RHS by pre-multiplying the matrix with the 
	// prescribed value array U and then adding it to F
	void AssembleRHS(vector<int>& lm, matrix& ke, vector<double>& U);

	// This assembles a vetor to the RHS
	void AssembleRHS(vector<int>& lm, vector<double>& fe);

private:
	FEModel*		m_fem;
	FEGlobalMatrix& m_K;	//!< The global stiffness matrix
	vector<double>&	m_F;	//!< Contributions from prescribed degrees of freedom
	vector<double>&	m_u;	//!< the array with prescribed values
};
