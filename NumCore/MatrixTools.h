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
#include <ostream>
#include <FECore/CompactUnSymmMatrix.h>
#include "numcore_api.h"

namespace NumCore
{
	// write a matrix to file
	// mode:
	// - 0 = binary mode
	// - 1 = text mode
	NUMCORE_API bool write_hb(CompactMatrix& K, const char* szfile, int mode = 0);

	// write a vector to file
	NUMCORE_API bool write_vector(const std::vector<double>& a, const char* szfile, int mode = 0);

	// read compact matrix (binary mode only)
	NUMCORE_API CompactMatrix* read_hb(const char* szfile);

	// read vector<double> from file
	NUMCORE_API bool read_vector(std::vector<double>& a, const char* szfile);

	// calculate inf-norm of inverse matrix (only works with CRSSparsMatrix(1))
	NUMCORE_API double inverse_infnorm(CompactMatrix* A);

	// calculate condition number of a CRSSparseMatrix(1) (Very expensive!)
	NUMCORE_API double conditionNumber(CRSSparseMatrix* A);

	// estimate condition number
	NUMCORE_API double estimateConditionNumber(SparseMatrix* A);

	// create a random vector
	NUMCORE_API void randomVector(std::vector<double>& R, double vmin = 0.0, double vmax = 1.0);

	// inf-norm of a vector
	NUMCORE_API double infNorm(const std::vector<double>& x);

	// 1-norm of a vector
	NUMCORE_API double oneNorm(const std::vector<double>& x);

	// print matrix sparsity pattern to svn file
	NUMCORE_API void print_svg(CompactMatrix* m, std::ostream &out, int i0 = 0, int j0 = 0, int i1 = -1, int j1 = -1);

} // namespace NumCore
