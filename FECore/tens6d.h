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
// E111111	d[  0]	E112111	d[ 18]	E113111	d[ 36]			
// E111112	d[  1]	E112112	d[ 19]	E113112	d[ 37]	
// E111113	d[  2]	E112113	d[ 20]	E113113	d[ 38]	
// E111122	d[  3]	E112122	d[ 21]	E113122	d[ 39]	
// E111123	d[  4]	E112123	d[ 22]	E113123	d[ 40]	
// E111133	d[  5]	E112133	d[ 23]	E113133	d[ 41]	
// E111211	d[  6]	E112211	d[ 24]	E113211	d[ 42]	
// E111212	d[  7]	E112212	d[ 25]	E113212	d[ 43]	
// E111213	d[  8]	E112213	d[ 26]	E113213	d[ 44]	
// E111222	d[  9]	E112222	d[ 27]	E113222	d[ 45]	
// E111223	d[ 10]	E112223	d[ 28]	
// E111233	d[ 11]	E112233	d[ 29]	
// E111311	d[ 12]	E112311	d[ 30]	
// E111312	d[ 13]	E112312	d[ 31]	
// E111313	d[ 14]	E112313	d[ 32]	
// E111322	d[ 15]	E112322	d[ 33]	
// E111323	d[ 16]	E112323	d[ 34]	
// E111333	d[ 17]	E112333	d[ 35]	

class tens6ds
{
public:
	enum { NNZ = 46 };

	// default constructor
	tens6ds(){}
	
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