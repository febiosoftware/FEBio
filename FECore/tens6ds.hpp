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
// NOTE: This file is automatically included from tens6d.h
// Users should not include this file manually!

// access operator
inline double tens6ds::operator() (int i, int j, int k, int l, int m, int n)
{
	// lookup table convering triplets to a row/column index
	const int LUT[3][3][3] = {
		{{0,1,2},{1,3,4},{2,4,5}},
		{{1,3,4},{3,6,7},{4,7,8}},
		{{2,4,5},{4,7,8},{5,8,9}}};

	// index to start of columns
	const int M[10] = {0,1,3,6,10,15,21,28,37,46};

	int I = LUT[i][j][k];
	int J = LUT[l][m][n];
	return (I <= J ? d[M[J]+I] : d[M[I]+J]);
}
