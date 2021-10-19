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
#include "tensor_base.h"

//-----------------------------------------------------------------------------
// Tensor classes defined in this file
class tens6ds;
class tens6d;

//-----------------------------------------------------------------------------
// traits for these classes defining the number of components
template <> class tensor_traits<tens6ds>{public: enum { NNZ =  55}; };
template <> class tensor_traits<tens6d >{public: enum { NNZ = 729}; };

//-----------------------------------------------------------------------------
//! Class for 6th order tensors. This class assumes the following symmetries:
//  - full (minor) symmetry in the first three legs 
//  - full (minor) symmetry in the last three legs 
//  - major symmetry (A[i,j,k;l,m,n] = A[l,m,n;i,j,k])
//
// This tensor is effectively stored as an upper-triangular 10x10 matrix
// using column major ordering.
//    
//       / 0  1  3  6 10 15 21 28 37 46 \     / A00 A01 A02  ...   A09 \
//      |     2  4  7 11 16 22 29 38 47  |   |      A11 A12  ...   A19  |
//      |        5  8 12 17 23 30 39 48  |   |          A22  ...   A18  |
//      |           9 13 18 24 31 40 49  |   |            .             |
//      |             14 19 25 32 41 50  |   |              .           |
//  A = |                20 26 33 42 51  | = |                .         |
//      |                   27 34 43 52  |   |                          |
//      |                      35 44 53  |   |                          |
//      |                         45 54  |   |                          |
//      \                            55 /     \                    A99 / 
//
//  where A[I,J] = A[i,j,k;l,m,n] using the following convention
//
//    I/J  |  i/l   j/m    k/n
// --------+-------------------
//     0   |   0      0      0
//     1   |   0      0      1
//     2   |   0      0      2
//     3   |   0      1      1
//     4   |   0      1      2
//     5   |   0      2      2
//     6   |   1      1      1
//     7   |   1      1      2
//     8   |   1      2      2
//     9   |   2      2      2
//
class tens6ds : public tensor_base<tens6ds>
{
public:
	// constructors
	tens6ds(){}

public:
	// access operator
	double operator () (int i, int j, int k, int l, int m, int n);
};

void calculate_e2O(tens6ds& e, double K[3][3], double Ri[3], double Rj[3] );


//-----------------------------------------------------------------------------
// class for general 6-th order tensors. No assumed symmetries
class tens6d : public tensor_base<tens6d>
{
public:
	// default constructor
	tens6d() {}

	// access operators
	double  operator () (int i, int j, int k, int l, int m, int n) const;
	double& operator () (int i, int j, int k, int l, int m, int n);
};

// The following file contains the actual definition of the class functions
#include "tens6d.hpp"
#include "tens6ds.hpp"
