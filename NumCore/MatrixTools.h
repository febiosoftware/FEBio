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
#include <ostream>
#include "CompactUnSymmMatrix.h"

namespace NumCore
{
	// write a matrix to file
	bool write_hb(CompactMatrix& K, const char* szfile);

	// write a vector to file
	bool write_vector(const vector<double>& a, const char* szfile);

	// calculate inf-norm of inverse matrix (only works with CRSSparsMatrix(1))
	double inverse_infnorm(CompactMatrix* A);

	// calculate condition number of a CRSSparseMatrix(1) (Very expensive!)
	double conditionNumber(CRSSparseMatrix* A);

	// estimate condition number
	double estimateConditionNumber(SparseMatrix* A);

	// create a random vector
	void randomVector(vector<double>& R, double vmin = 0.0, double vmax = 1.0);

	// inf-norm of a vector
	double infNorm(const std::vector<double>& x);

	// print matrix sparsity pattern to svn file
	void print_svg(CompactMatrix* m, std::ostream &out, int i0 = 0, int j0 = 0, int i1 = -1, int j1 = -1);

} // namespace NumCore
