#pragma once
#include "tensor_base.h"

//-----------------------------------------------------------------------------
//! Class for 5th order tensors with full symmetry (any pair of indices can be swapped)
//
class tens5ds : public tensor_base<tens5ds, 21>
{
public:
	enum { NNZ = 21 };

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
class tens5d : public tensor_base<tens5d, 243>
{
public:
	enum { NNZ = 243 };

public:
	// default constructor
	tens5d() {}

	// access operators
	double operator () (int i, int j, int k, int l, int m) const;
	double& operator () (int i, int j, int k, int l, int m);
};

void calculate_d2O(tens5ds& d, double K[3][3], double Ri[3], double Rj[3] );

// The following file contains the actual definition of the class functions
#include "tens5d.hpp"
#include "tens5ds.hpp"
