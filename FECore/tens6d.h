#pragma once
#include "tensor_base.h"

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
class tens6ds : public tensor_base<tens6ds, 55>
{
public:
	enum { NNZ = 55 };

public:
	// constructors
	tens6ds(){}

public:
	// access operator
	double operator () (int i, int j, int k, int l, int m, int n);
};

// The following file contains the actual definition of the class functions
#include "tens6ds.hpp"
