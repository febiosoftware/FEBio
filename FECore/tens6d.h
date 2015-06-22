#pragma once

#include "mat3d.h"
#include "tens3d.h"

//-----------------------------------------------------------------------------
//! Class for 6th order tensor based on symmetry requirements for second-order hyperelasticity:
//  - E = d2Phi/dGijkdGlmn  
//  - full symmetry in the first three legs 
//  - full symmetry in the last three legs 
//  - major-like symmetry when switching the order of the partial derivatives (Eijklmn = Elmnijk)
//  - 46 out of 729 compenents are unique
//
// E111111	d[  0]	E113123	d[ 18]	E133233	d[ 36]			
// E111112	d[  1]	E113133	d[ 19]	E133333	d[ 37]	
// E111113	d[  2]	E113222	d[ 20]	E222222	d[ 38]	
// E111122	d[  3]	E113223	d[ 21]	E222223	d[ 39]	
// E111123	d[  4]	E113233	d[ 22]	E222233	d[ 40]	
// E111133	d[  5]	E113333	d[ 23]	E222333	d[ 41]	
// E111222	d[  6]	E122133	d[ 24]	E223233	d[ 42]	
// E111223	d[  7]	E122222	d[ 25]	E223333	d[ 43]	
// E111233	d[  8]	E122223	d[ 26]	E233333	d[ 44]	
// E111333	d[  9]	E122233	d[ 27]	E333333	d[ 45]	
// E112122	d[ 10]	E122333	d[ 28]	
// E112123	d[ 11]	E123133	d[ 29]	
// E112133	d[ 12]	E123222	d[ 30]	
// E112222	d[ 13]	E123223	d[ 31]	
// E112223	d[ 14]	E123233	d[ 32]	
// E112233	d[ 15]	E123333	d[ 33]	
// E112333	d[ 16]	E133222	d[ 34]	
// E113122	d[ 17]	E133223	d[ 35]	

class tens6ds
{
public:
	enum { NNZ = 46 };

	// default constructor
	tens6ds(){zero();}
	
	tens6ds(const double g)
	{
		for (int i = 0; i < NNZ; i++)
			d[i] = g;
	}

	tens6ds(double m[46])
	{
		for (int i = 0; i < NNZ; i++)
			d[i] = m[i];
	}

	// arithmetic operators
	tens6ds operator + (const tens6ds& t) const;
	tens6ds operator - (const tens6ds& t) const;
	tens6ds operator * (double g) const;
	tens6ds operator / (double g) const;

	// arithmetic assignment operators
	tens6ds& operator += (const tens6ds& t);
	tens6ds& operator -= (const tens6ds& t);
	tens6ds& operator *= (double g);
	tens6ds& operator /= (double g);

	// unary operators
	tens6ds operator - () const;
	
	// initialize to zero
	void zero();

	//void unit();

	//tens3drs contract3s(tens3drs H);


public:
	double d[NNZ];	// stored in column major order
};

// The following file contains the actual definition of the class functions
#include "tens6ds.hpp"