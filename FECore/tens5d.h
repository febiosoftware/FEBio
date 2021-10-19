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
// Classes defined in this file
class tens5ds;
class tens5d;

//-----------------------------------------------------------------------------
// traits for these classes defining the number of components
template <> class tensor_traits<tens5ds> {public: enum { NNZ =  21}; };
template <> class tensor_traits<tens5d > {public: enum { NNZ = 243}; };

//-----------------------------------------------------------------------------
//! Class for 5th order tensors with full symmetry (any pair of indices can be swapped)
//
class tens5ds : public tensor_base<tens5ds>
{
public:
	// default constructor
	tens5ds() {}

	// access operator
	// TODO: implement this
//	double operator () (int i, int j, int k, int l, int m) const;
};

//-----------------------------------------------------------------------------
//! Class for 5th order tensors (no symmetries)
//! This is stored as a 27x9 matrix
//
//      /  0  27  54  81  108  135 162 189 216  \
//  A = |  .                                 .  |
//      |  .                                 .  |
//      |  .                                 .  |
//      \  26                               242 /
//
class tens5d : public tensor_base<tens5d>
{
public:
	// default constructor
	tens5d() {}

	// access operators
	double  operator () (int i, int j, int k, int l, int m) const;
	double& operator () (int i, int j, int k, int l, int m);
};

// The following file contains the actual definition of the class functions
#include "tens5d.hpp"
#include "tens5ds.hpp"
