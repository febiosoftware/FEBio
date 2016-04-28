#pragma once

//-----------------------------------------------------------------------------
//! Class for 5th order tensor based on symmetry requirements for second-order hyperelasticity and Cauchy stress:
//  - D = d2Phi/ddDijdHklm  
//  - symmetry in the first two legs (Dij = Dji) 
//  - fullsymmetry in the last two legs (Tijk = Tjik = Tkji = Tikj = Tkij = Tjki), 
//  - major-like symmetry when switching the order of the partial derivatives (Dijklm = Dklmij)
//  - 29 out of 243 compenents are unique
//
// D11111	d[  0]	D13222	d[ 18]				
// D11112	d[  1]	D13223	d[ 19]		
// D11113	d[  2]	D13233	d[ 20]		
// D11122	d[  3]	D13333	d[ 21]	
// D11123	d[  4]	D22222	d[ 22]	
// D11133	d[  5]	D22223	d[ 23]	
// D11222	d[  6]	D22233	d[ 24]	
// D11223	d[  7]	D22333	d[ 25]		
// D11233	d[  8]	D23233	d[ 26]	
// D11333	d[  9]	D23333	d[ 27]		
// D12122	d[ 10]	D33333	d[ 28]		
// D12123	d[ 11]	
// D12133	d[ 12]			
// D12222	d[ 13]		
// D12223	d[ 14]		
// D12233	d[ 15]		
// D12333	d[ 16]			
// D13133	d[ 17]			


class tens5ds
{
public:
	enum { NNZ = 29 };

	// default constructor
	tens5ds() {}

	// access operator
	double operator () (int i, int j, int k, int l, int m) const;

	// arithmetic operators
	tens5ds operator + (const tens5ds& t) const;
	tens5ds operator - (const tens5ds& t) const;
	tens5ds operator * (double g) const;
	tens5ds operator / (double g) const;

	// arithmetic assignment operators
	tens5ds& operator += (const tens5ds& t);
	tens5ds& operator -= (const tens5ds& t);
	tens5ds& operator *= (double g);
	tens5ds& operator /= (double g);

	// unary operators
	tens5ds operator - () const;
	
	// initialize to zero
	void zero();

public:
	double d[NNZ];	// stored in column major order
};

// The following file contains the actual definition of the class functions
#include "tens5ds.hpp"
