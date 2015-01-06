#pragma once

#include "mat3d.h"
#include "tens3d.h"

//-----------------------------------------------------------------------------
//! Class for 5th order tensor based on symmetry requirements for second-order hyperelasticity:
//  - D = d2Phi/dFijdGklm  
//  - symmetry in the first two legs (Dijklm = Djiklm) 
//  - symmetry in the last two legs (Dijklm = Dijkml), 
//  - major-like symmetry when switching the order of the partial derivatives (Dijklm = Dklmij)
//  - 108 out of 243 compenents are unique
//
// D11111	d[  0]	D12111	d[ 18]	D13111	d[ 36]	D22111	d[ 54]	D23111	d[ 72]	D33111	d[ 90]		
// D11112	d[  1]	D12112	d[ 19]	D13112	d[ 37]	D22112	d[ 55]	D23112	d[ 73]	D33112	d[ 91]
// D11113	d[  2]	D12113	d[ 20]	D13113	d[ 38]	D22113	d[ 56]	D23113	d[ 74]	D33113	d[ 92]
// D11122	d[  3]	D12122	d[ 21]	D13122	d[ 39]	D22122	d[ 57]	D23122	d[ 75]	D33122	d[ 93]
// D11123	d[  4]	D12123	d[ 22]	D13123	d[ 40]	D22123	d[ 58]	D23123	d[ 76]	D33123	d[ 94]
// D11133	d[  5]	D12133	d[ 23]	D13133	d[ 41]	D22133	d[ 59]	D23133	d[ 77]	D33133	d[ 95]
// D11211	d[  6]	D12211	d[ 24]	D13211	d[ 42]	D22211	d[ 60]	D23211	d[ 78]	D33211	d[ 96]
// D11212	d[  7]	D12212	d[ 25]	D13212	d[ 43]	D22212	d[ 61]	D23212	d[ 79]	D33212	d[ 97]
// D11213	d[  8]	D12213	d[ 26]	D13213	d[ 44]	D22213	d[ 62]	D23213	d[ 80]	D33213	d[ 98]
// D11222	d[  9]	D12222	d[ 27]	D13222	d[ 45]	D22222	d[ 63]	D23222	d[ 81]	D33222	d[ 99]
// D11223	d[ 10]	D12223	d[ 28]	D13223	d[ 46]	D22223	d[ 64]	D23223	d[ 82]	D33223	d[100]
// D11233	d[ 11]	D12233	d[ 29]	D13233	d[ 47]	D22233	d[ 65]	D23233	d[ 83]	D33233	d[101]
// D11311	d[ 12]	D12311	d[ 30]	D13311	d[ 48]	D22311	d[ 66]	D23311	d[ 84]	D33311	d[102]
// D11312	d[ 13]	D12312	d[ 31]	D13312	d[ 49]	D22312	d[ 67]	D23312	d[ 85]	D33312	d[103]
// D11313	d[ 14]	D12313	d[ 32]	D13313	d[ 50]	D22313	d[ 68]	D23313	d[ 86]	D33313	d[104]
// D11322	d[ 15]	D12322	d[ 33]	D13322	d[ 51]	D22322	d[ 69]	D23322	d[ 87]	D33322	d[105]
// D11323	d[ 16]	D12323	d[ 34]	D13323	d[ 52]	D22323	d[ 70]	D23323	d[ 88]	D33323	d[106]
// D11333	d[ 17]	D12333	d[ 35]	D13333	d[ 53]	D22333	d[ 71]	D23333	d[ 89]	D33333	d[107]

class tens5d
{
public:
	enum { NNZ = 108 };

	// default constructor
	tens5d(){}
	
	tens5d(const double g)
	{
		for (int i = 0; i < NNZ; i++)
			d[i] = g;
	}

	tens5d(double m[108])
	{
		for (int i = 0; i < NNZ; i++)
			d[i] = m[i];
	}

	// arithmetic operators
	tens5d operator + (const tens5d& t) const;
	tens5d operator - (const tens5d& t) const;
	tens5d operator * (double g) const;
	tens5d operator / (double g) const;

	// arithmetic assignment operators
	tens5d& operator += (const tens5d& t);
	tens5d& operator -= (const tens5d& t);
	tens5d& operator *= (double g);
	tens5d& operator /= (double g);

	// unary operators
	tens5d operator - () const;
	
	// initialize to zero
	void zero();

	void unit();

	mat3ds contract3rs(tens3drs H);
	tens3drs contract2s(mat3ds s);

public:
	double d[NNZ];	// stored in column major order
};

// The following file contains the actual definition of the class functions
#include "tens5d.hpp"