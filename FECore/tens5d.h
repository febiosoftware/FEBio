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
