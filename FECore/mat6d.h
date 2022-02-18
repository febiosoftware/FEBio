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

//-----------------------------------------------------------------------------
// Class describing a symmetric 6x6 matrix.
// Only the upper triangular matrix is stored in column major order
//     / 0   1   3   6   10   15  \
//     |     2   4   7   11   16  |
//     |         5   8   12   17  |
// A = |             9   13   18  |
//     |                 14   19  |
//     \                      20  /

class mat6ds
{
public:
	enum { NNZ = 21 };

public:
	mat6ds() {}

	//! operators
	mat6ds& operator *= (double g);

public:
	// access operators
	double& operator () (int i, int j);
	const double& operator () (int i, int j) const;

public:
	//! initialize to zero
	void zero();

private:
	double	d[NNZ];
};

//-----------------------------------------------------------------------------
// Class describing a 6x6 matrix
class mat6d
{
public:
	//! default constructor
	mat6d() {}

	//! operators
	mat6d& operator *= (double g);

public:
	// access operator
	double* operator [] (int i) { return d[i]; }

public:
	//! initialize to zero
	void zero();

private:
	double	d[6][6];	//!< matrix data
};

//-----------------------------------------------------------------------------
inline double& mat6ds::operator()(int i, int j)
{
	const int n[] = {0, 1, 3, 6, 10, 15};
	return (i<=j? d[n[j]+i] : d[n[i]+j]);
}

//-----------------------------------------------------------------------------
inline const double& mat6ds::operator()(int i, int j) const
{
	const int n[] = {0, 1, 3, 6, 10, 15};
	return (i<=j? d[n[j]+i] : d[n[i]+j]);
}

//-----------------------------------------------------------------------------
inline mat6ds& mat6ds::operator *= (double g)
{
	d[ 0] *= g; d[ 1] *= g; d[ 2] *= g; d[ 3] *= g; d[ 4] *= g; d[ 5] *= g;
	d[ 6] *= g; d[ 7] *= g; d[ 8] *= g; d[ 9] *= g; d[10] *= g;
	d[11] *= g; d[12] *= g; d[13] *= g; d[14] *= g;
	d[15] *= g; d[16] *= g; d[17] *= g;
	d[18] *= g; d[19] *= g;
	d[20] *= g;
	return *this;
}

//-----------------------------------------------------------------------------
inline void mat6ds::zero()
{
	d[ 0] = d[ 1] = d[ 2] = d[ 3] = d[ 4] = d[ 5] = 0.0;
	d[ 6] = d[ 7] = d[ 8] = d[ 9] = d[10] = 0.0;
	d[11] = d[12] = d[13] = d[14] = 0.0;
	d[15] = d[16] = d[17] = 0.0;
	d[18] = d[19] = 0.0;
	d[20] = 0.0;
}

//-----------------------------------------------------------------------------
inline void mat6d::zero()
{
	d[0][0] = d[0][1] = d[0][2] = d[0][3] = d[0][4] = d[0][5] = 0.0;
	d[1][0] = d[1][1] = d[1][2] = d[1][3] = d[1][4] = d[1][5] = 0.0;
	d[2][0] = d[2][1] = d[2][2] = d[2][3] = d[2][4] = d[2][5] = 0.0;
	d[3][0] = d[3][1] = d[3][2] = d[3][3] = d[3][4] = d[3][5] = 0.0;
	d[4][0] = d[4][1] = d[4][2] = d[4][3] = d[4][4] = d[4][5] = 0.0;
	d[5][0] = d[5][1] = d[5][2] = d[5][3] = d[5][4] = d[5][5] = 0.0;
}

//-----------------------------------------------------------------------------
mat6d& mat6d::operator *= (double g)
{
	d[0][0] *= g; d[0][1] *= g; d[0][2] *= g; d[0][3] *= g; d[0][4] *= g; d[0][5] *= g;
	d[1][0] *= g; d[1][1] *= g; d[1][2] *= g; d[1][3] *= g; d[1][4] *= g; d[1][5] *= g;
	d[2][0] *= g; d[2][1] *= g; d[2][2] *= g; d[2][3] *= g; d[2][4] *= g; d[2][5] *= g;
	d[3][0] *= g; d[3][1] *= g; d[3][2] *= g; d[3][3] *= g; d[3][4] *= g; d[3][5] *= g;
	d[4][0] *= g; d[4][1] *= g; d[4][2] *= g; d[4][3] *= g; d[4][4] *= g; d[4][5] *= g;
	d[5][0] *= g; d[5][1] *= g; d[5][2] *= g; d[5][3] *= g; d[5][4] *= g; d[5][5] *= g;
	return *this;
}
