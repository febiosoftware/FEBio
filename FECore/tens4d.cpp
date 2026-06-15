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



#include "stdafx.h"
#include "tens4d.h"
#include <math.h>


//-----------------------------------------------------------------------------
//! This function checks the positive definiteness of a 4th order tensor
//! having both major and minor symmetries. The function does not do an
//! exhaustive test, in the sense it can only detect failure. If a tensor passes
//! it is not guaranteed that the tensor is indeed positive-definite.
bool IsPositiveDefinite(const tens4ds& t)
{
	// test 1. all diagonal entries have to be positive
	if (t(0,0) <= 0) return false;
	if (t(1,1) <= 0) return false;
	if (t(2,2) <= 0) return false;
	if (t(3,3) <= 0) return false;
	if (t(4,4) <= 0) return false;
	if (t(5,5) <= 0) return false;

	// test 2. t(i,i)+t(j,j) > 2t(i,j)
	int i, j;
	for (i=0; i<6; ++i)
	{
		for (j=i+1; j<6; ++j)
		{
			if (t(i,i)+t(j,j) <= 2*t(i,j))
			{
				return false;
			}
		}
	}

	// test 3. the element with largest modulus lies on the main diagonal
	double l = -1, v;
	bool d = false;
	for (i=0; i<6; ++i)
	{
		for (j=i; j<6; ++j)
		{
			v = fabs(t(i,j));
			if (v > l)
			{
				l = v;
				d = (i==j);
			}
		}
	}

	if (d == false)
	{
		return false;
	}

	// if all tests pass, it is not guaranteed that the tensor is indeed positive-definite
	// but we'd have some good reasons to believe so.
	return true;
}

//-----------------------------------------------------------------------------
// (a ddot b)_ijkl = a_ijmn b_mnkl
tens4d ddot(const tens4d& a, const tens4d& b)
{
    double A[9][9] = {0};
    double B[9][9] = {0};
    double C[9][9] = {0};
    tens4d aa = a, bb = b;
    aa.extract(A); bb.extract(B);
    
    // perform matrix multiplication
    for (int  i=0; i<9; ++i) {
        for (int j=0; j<9; ++j) {
            C[i][j] = 0;
            for (int k=0; k<9; ++k) {
                C[i][j] += A[i][k]*B[k][j];
            }
        }
    }
    tens4d c(C);
    return c;
    
}

//-----------------------------------------------------------------------------
// inverse
tens4dmm tens4dmm::inverse() const
{
    // convert this tens4dmm into a tens4d object
    const tens4d t(*this);
    // invert the tens4d object
    const tens4d ti = t.inverse();
    
    // debug check
//    tens4d tIoI = ddot(t,ti);
    
    // create a tens4dmm object from this ten4d inverse
    tens4dmm S = ti.mmsymm();
    return S;
}
