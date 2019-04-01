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
#include "tensor_base.h"
#include "mat3d.h"

//-----------------------------------------------------------------------------
class tens4ds;
class tens4dms;
class tens4d;

//-----------------------------------------------------------------------------
//! Class for 4th order tensors with major and minor symmetries (i.e., super-symmetry)

// Due to the major symmetry we can store this tensor as a 6x6 matrix.
// The tensor is stored in column major order:
//
//     / 0   1   3   6   10   15  \   / C0000  C0011  C0022  C0001  C0012  C0002 \
//     |     2   4   7   11   16  |   |        C1111  C1122  C1101  C1112  C1102 |
//     |         5   8   12   17  |   |               C2222  C2201  C2212  C2202 |
// A = |             9   13   18  | = |                      C0101  C0112  C0102 |
//     |                 14   19  |   |                             C1212  C1202 |
//     \                      20  /   \                                    C0202 /
//
// Note that due to the major symmetry we only store the upper triangular matrix of this tensor
//

class tens4ds
{
public:
	enum { NNZ = 21 };

	// default constructor
	tens4ds() {}
	explicit tens4ds(const double g);
	tens4ds(double m[6][6]);

	double& operator () (int i, int j, int k, int l);
	double operator () (int i, int j, int k, int l) const;

	// TODO: remove this?
	double& operator () (int i, int j);
	double operator () (int i, int j) const;

	// arithmetic operators
	tens4ds operator + (const tens4ds& t) const;
	tens4ds operator - (const tens4ds& t) const;
	tens4ds operator * (double g) const;
	tens4ds operator / (double g) const;

	// arithmetic assignment operators
	tens4ds& operator += (const tens4ds& t);
	tens4ds& operator -= (const tens4ds& t);
	tens4ds& operator *= (double g);
	tens4ds& operator /= (double g);

	// unary operators
	tens4ds operator - () const;
	
	// double dot product with 2nd order tensor
	mat3ds dot(const mat3ds& m) const;
    mat3ds dot(const mat3dd& m) const { return dot(mat3ds(m)); }
    mat3ds dot(const mat3d& m) const;
    mat3ds dot2(const mat3d& m) const;

	// trace
	double tr() const;
	
	// initialize to zero
	void zero();

	// extract 6x6 matrix
	void extract(double d[6][6]);

	// calculates the inverse
	tens4ds inverse() const;
    
    // evaluate push/pull operation
    tens4ds pp(const mat3d& F);

public:
	double d[NNZ];	// stored in column major order
};

//! Check positive definiteness of a 4th-order symmetric tensor
bool IsPositiveDefinite(const tens4ds& t);

// outer (dyadic) products for symmetric matrices
tens4ds dyad1s(const mat3dd& a);
tens4ds dyad1s(const mat3ds& a);
tens4ds dyad1s(const mat3dd& a, const mat3dd& b); 
tens4ds dyad1s(const mat3ds& a, const mat3dd& b);
tens4ds dyad1s(const mat3ds& a, const mat3ds& b);
inline tens4ds dyad1s(const mat3dd& a, const mat3ds& b) { return dyad1s(b, a); }
tens4ds dyad2s(const mat3ds& a);
tens4ds dyad2s(const mat3ds& a, const mat3ds& b);
tens4ds dyad4s(const mat3dd& a);
tens4ds dyad4s(const mat3ds& a);
tens4ds dyad4s(const mat3ds& a, const mat3ds& b);
tens4ds dyad4s(const mat3ds& a, const mat3dd& b);
tens4ds dyad4s(const vec3d& a, const mat3d& K, const vec3d& b);
tens4ds ddots(const tens4ds& a, const tens4ds& b);
mat3d vdotTdotv(const vec3d& a, const tens4ds& T, const vec3d& b);

inline tens4ds operator * (const double g, const tens4ds& a) { return a*g; }

// The following file contains the actual definition of the class functions
#include "tens4ds.hpp"

//-----------------------------------------------------------------------------
//! Class for 4th order tensors with major symmetry only

// Due to the lack of minor symmetry, we have to store additional components of this tensor as a 9x9 matrix.
// Major symmetry ensures that this storage matrix is symmetric about its main diagonal.
// The tensor is stored in column major order:
//
//     / 0   1   3   6   10   15  21   28   36  \
//     |     2   4   7   11   16  22   29   37  |
//     |         5   8   12   17  23   30   38  |
// A = |             9   13   18  24   31   39  |
//     |                 14   19  25   32   40  |
//     |                      20  26   33   41  |
//     |                          27   34   42  |
//     |                               35   43  |
//     \                                    44  /
//


class tens4dms
{
public:
	enum { NNZ = 45 };

	// Default constructor
	tens4dms() {}
	tens4dms(const double g);
	tens4dms(double m[9][9]);

	// access operators
	double& operator () (int i, int j, int k, int l);
	double operator () (int i, int j, int k, int l) const;
	double& operator () (int i, int j);
	double operator () (int i, int j) const;

	// arithmetic operators
	tens4dms operator + (const tens4dms& t) const;
	tens4dms operator - (const tens4dms& t) const;
	tens4dms operator * (double g) const;
	tens4dms operator / (double g) const;

	// arithmetic assignment operators
	tens4dms& operator += (const tens4dms& t);
	tens4dms& operator -= (const tens4dms& t);
	tens4dms& operator *= (double g);
	tens4dms& operator /= (double g);

	// unary operators
	tens4dms operator - () const;
	
	// trace
	double tr() const;
	
	// initialize to zero
	void zero();

	// extract 9x9 matrix
	void extract(double d[9][9]);

	// compute the super-symmetric (major and minor symmetric) component of the tensor
	tens4ds supersymm() const;

	//// calculates the inverse
	//tens4dms inverse() const;

public:
	double d[NNZ];	// stored in column major order
};

// outer (dyadic) products for symmetric and non-symmetric matrices
tens4dms dyad1(const mat3d& a);
tens4dms dyad1(const mat3ds& a, const mat3ds& b);

// The following file contains the actual definition of the class functions
#include "tens4dms.hpp"

//-----------------------------------------------------------------------------
//! Class for 4th order tensors without symmetry

// We store components of this tensor as a 9x9 matrix.
// The tensor is stored in column major order:
//
//       00  11  22   01   12   02   10   21   20    |
//       --------------------------------------------+---
//     / 0   9   18   27   36   45   54   63   72  \ | 00
//     | 1   10  19   28   37   46   55   64   73  | | 11
//     | 2   11  20   29   38   47   56   65   74  | | 22
// A = | 3   12  21   30   39   48   57   66   75  | | 01
//     | 4   13  22   31   40   49   58   67   76  | | 12
//     | 5   14  23   32   41   50   59   68   77  | | 02
//     | 6   15  24   33   42   51   60   69   78  | | 10
//     | 7   16  25   34   43   52   61   70   79  | | 21
//     \ 8   17  26   35   44   53   62   71   80  / | 20
//

template <> class tensor_traits<tens4d> {public: enum { NNZ = 81}; };

class tens4d : public tensor_base<tens4d>
{
public:
	// constructors
	tens4d() {}

public:
	// access operators
	double& operator () (int i, int j, int k, int l);
	double operator () (int i, int j, int k, int l) const;

private:
	double& operator () (int i, int j);
	double operator () (int i, int j) const;
};

// The following file contains the actual definition of the class functions
#include "tens4d.hpp"
