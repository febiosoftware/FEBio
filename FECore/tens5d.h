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
	double operator () (int i, int j, int k, int l, int m) const;
};

//-----------------------------------------------------------------------------
//! Class for 5th order tensors (no symmetries)
//
class tens5d  : public tensor_base<tens5ds, 243>
{
public:
	enum { NNZ = 243 };

public:
	// default constructor
	tens5d() {}

	// access operator
	double operator () (int i, int j, int k, int l, int m) const;
};


// The following file contains the actual definition of the class functions
#include "tens5d.hpp"
#include "tens5ds.hpp"
