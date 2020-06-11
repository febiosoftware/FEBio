/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2020 University of Utah, The Trustees of Columbia University in 
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
// NOTE: This file is automatically included from tens4d.h
// Users should not include this file manually!

#include "matrix.h"

inline tens4dms::tens4dms(const double g)
{
	for (int i = 0; i < NNZ; i++)
		d[i] = g;
}

inline tens4dms::tens4dms(double m[9][9])
{
	d[ 0] = m[0][0];
	d[ 1] = m[0][1]; d[ 2] = m[1][1];
	d[ 3] = m[0][2]; d[ 4] = m[1][2]; d[ 5] = m[2][2];
	d[ 6] = m[0][3]; d[ 7] = m[1][3]; d[ 8] = m[2][3]; d[ 9] = m[3][3];
	d[10] = m[0][4]; d[11] = m[1][4]; d[12] = m[2][4]; d[13] = m[3][4]; d[14] = m[4][4];
	d[15] = m[0][5]; d[16] = m[1][5]; d[17] = m[2][5]; d[18] = m[3][5]; d[19] = m[4][5]; d[20] = m[5][5];
	d[21] = m[0][6]; d[22] = m[1][6]; d[23] = m[2][6]; d[24] = m[3][6]; d[25] = m[4][6]; d[26] = m[5][6]; d[27] = m[6][6];
	d[28] = m[0][7]; d[29] = m[1][7]; d[30] = m[2][7]; d[31] = m[3][7]; d[32] = m[4][7]; d[33] = m[5][7]; d[34] = m[6][7]; d[35] = m[7][7];
	d[36] = m[0][8]; d[37] = m[1][8]; d[38] = m[2][8]; d[39] = m[3][8]; d[40] = m[4][8]; d[41] = m[5][8]; d[42] = m[6][8]; d[43] = m[7][8]; d[44] = m[8][8];
}

inline double& tens4dms::operator () (int i, int j, int k, int l)
{
	const int m[3][3] = {{0,3,5},{6,1,4},{8,7,2}};
	tens4dms& T = (*this);
	return T(m[i][j], m[k][l]);
}

inline double tens4dms::operator () (int i, int j, int k, int l) const
{
	const int m[3][3] = {{0,3,5},{6,1,4},{8,7,2}};
	const tens4dms& T = (*this);
	return T(m[i][j], m[k][l]);
}

inline double& tens4dms::operator () (int i, int j)
{
	const int m[9] = {0, 1, 3, 6, 10, 15, 21, 28, 36};
	if (i<=j) return d[m[j]+i]; else return d[m[i]+j];
}

inline double tens4dms::operator () (int i, int j) const
{
	const int m[9] = {0, 1, 3, 6, 10, 15, 21, 28, 36};
	if (i<=j) return d[m[j]+i]; else return d[m[i]+j];
}

// operator +
inline tens4dms tens4dms::operator + (const tens4dms& t) const
{
	tens4dms s;
	for (int i=0; i<NNZ; i++)
		s.d[i] = d[i] + t.d[i];
	
	return s;
}

// operator -
inline tens4dms tens4dms::operator - (const tens4dms& t) const
{
	tens4dms s;
	for (int i=0; i<NNZ; i++)
		s.d[i] = d[i] - t.d[i];

	return s;
}

// operator *
inline tens4dms tens4dms::operator * (double g) const
{
	tens4dms s;
	for (int i=0; i<NNZ; i++)
		s.d[i] = g*d[i];
	
	return s;
}

// operator /
inline tens4dms tens4dms::operator / (double g) const
{
	tens4dms s;
	for (int i=0; i<NNZ; i++)
		s.d[i] = d[i]/g;
	
	return s;
}

// assignment operator +=
inline tens4dms& tens4dms::operator += (const tens4dms& t)
{
	for (int i=0; i<NNZ; i++)
		d[i] += t.d[i];
	
	return (*this);
}

// assignment operator -=
inline tens4dms& tens4dms::operator -= (const tens4dms& t)
{
	for (int i=0; i<NNZ; i++)
		d[i] -= t.d[i];
	
	return (*this);
}

// assignment operator *=
inline tens4dms& tens4dms::operator *= (double g)
{
	for (int i=0; i<NNZ; i++)
		d[i] *= g;
	
	return (*this);
}

// assignment operator /=
inline tens4dms& tens4dms::operator /= (double g)
{
	for (int i=0; i<NNZ; i++)
		d[i] /= g;
	
	return (*this);
}

// unary operator -
inline tens4dms tens4dms::operator - () const
{
	tens4dms s;
	for (int i = 0; i < NNZ; i++)
		s.d[i] = -d[i];

	return s;
}

// trace
// C.tr() = I:C:I
// No need to change the trace operation because major symmetry still applies
inline double tens4dms::tr() const
{
	return (d[0]+d[2]+d[5]+2*(d[1]+d[3]+d[4]));
}

// intialize to zero
inline void tens4dms::zero()
{
	for (int i = 0; i < NNZ; i++)
		d[i] = 0;
}

// extract 9x9 matrix
inline void tens4dms::extract(double D[9][9])
{
	D[0][0] = d[0];  D[0][1] = d[1];  D[0][2] = d[3];  D[0][3] = d[6];  D[0][4] = d[10]; D[0][5] = d[15]; D[0][6] = d[21]; D[0][7] = d[28]; D[0][8] = d[36];
	D[1][0] = d[1];  D[1][1] = d[2];  D[1][2] = d[4];  D[1][3] = d[7];  D[1][4] = d[11]; D[1][5] = d[16]; D[1][6] = d[22]; D[1][7] = d[29]; D[1][8] = d[37];
	D[2][0] = d[3];  D[2][1] = d[4];  D[2][2] = d[5];  D[2][3] = d[8];  D[2][4] = d[12]; D[2][5] = d[17]; D[2][6] = d[23]; D[2][7] = d[30]; D[2][8] = d[38];
	D[3][0] = d[6];  D[3][1] = d[7];  D[3][2] = d[8];  D[3][3] = d[9];  D[3][4] = d[13]; D[3][5] = d[18]; D[3][6] = d[24]; D[3][7] = d[31]; D[3][8] = d[39];
	D[4][0] = d[10]; D[4][1] = d[11]; D[4][2] = d[12]; D[4][3] = d[13]; D[4][4] = d[14]; D[4][5] = d[19]; D[4][6] = d[25]; D[4][7] = d[32]; D[4][8] = d[40];
	D[5][0] = d[15]; D[5][1] = d[16]; D[5][2] = d[17]; D[5][3] = d[18]; D[5][4] = d[19]; D[5][5] = d[20]; D[5][6] = d[26]; D[5][7] = d[33]; D[5][8] = d[41];
	D[6][0] = d[21]; D[6][1] = d[22]; D[6][2] = d[23]; D[6][3] = d[24]; D[6][4] = d[25]; D[6][5] = d[26]; D[6][6] = d[27]; D[6][7] = d[34]; D[6][8] = d[42];
	D[7][0] = d[28]; D[7][1] = d[29]; D[7][2] = d[30]; D[7][3] = d[31]; D[7][4] = d[32]; D[7][5] = d[33]; D[7][6] = d[34]; D[7][7] = d[35]; D[7][8] = d[43];
	D[8][0] = d[36]; D[8][1] = d[37]; D[8][2] = d[38]; D[8][3] = d[39]; D[8][4] = d[40]; D[8][5] = d[41]; D[8][6] = d[42]; D[8][7] = d[43]; D[8][8] = d[44];
}

// compute the super symmetric (major and minor symmetric) component of the tensor
// Sijkl = (1/4)*(Cijkl + Cijlk + Cjikl + Cjilk)
inline tens4ds tens4dms::supersymm() const
{
	tens4ds s;

	s.d[0] = d[0]; s.d[1] = d[1]; s.d[3] = d[3]; s.d[6] = (d[6] + d[21])/2.;            s.d[10] = (d[10] + d[28])/2.;                 s.d[15] = (d[15] + d[36])/2.;
	               s.d[2] = d[2]; s.d[4] = d[4]; s.d[7] = (d[7] + d[22])/2.;            s.d[11] = (d[11] + d[29])/2.;                 s.d[16] = (d[16] + d[37])/2.;
				                  s.d[5] = d[5]; s.d[8] = (d[8] + d[23])/2.;            s.d[12] = (d[12] + d[30])/2.;                 s.d[17] = (d[17] + d[38])/2.;
								                 s.d[9] = (d[9] + 2.*d[24] + d[27])/4.; s.d[13] = (d[13] + d[31] + d[25] + d[34])/4.; s.d[18] = (d[18] + d[39] + d[26] + d[42])/4.;
												                                        s.d[14] = (d[14] + 2.*d[32] + d[35])/4.;      s.d[19] = (d[19] + d[40] + d[33] + d[43])/4.;
																			                                                          s.d[20] = (d[20] + 2.*d[41] + d[44])/4.;

	return s;
}

//-----------------------------------------------------------------------------
// (a dyad1s a)_ijkl = a_ij a_kl
inline tens4dms dyad1ms(const mat3d& a)
{
	tens4dms c;
	
	c.d[ 0] = a(0,0)*a(0,0);
	c.d[ 1] = a(0,0)*a(1,1);
	c.d[ 3] = a(0,0)*a(2,2);
	c.d[ 6] = a(0,0)*a(0,1);
	c.d[10] = a(0,0)*a(1,2);
	c.d[15] = a(0,0)*a(0,2);
	c.d[21] = a(0,0)*a(1,0);
	c.d[28] = a(0,0)*a(2,1);
	c.d[36] = a(0,0)*a(2,0);

	c.d[ 2] = a(1,1)*a(1,1);
	c.d[ 4] = a(1,1)*a(2,2);
	c.d[ 7] = a(1,1)*a(0,1);
	c.d[11] = a(1,1)*a(1,2);
	c.d[16] = a(1,1)*a(0,2);
	c.d[22] = a(1,1)*a(1,0);
	c.d[29] = a(1,1)*a(2,1);
	c.d[37] = a(1,1)*a(2,0);

	c.d[ 5] = a(2,2)*a(2,2);
	c.d[ 8] = a(2,2)*a(0,1);
	c.d[12] = a(2,2)*a(1,2);
	c.d[17] = a(2,2)*a(0,2);
	c.d[23] = a(2,2)*a(1,0);
	c.d[30] = a(2,2)*a(2,1);
	c.d[38] = a(2,2)*a(2,0);

	c.d[ 9] = a(0,1)*a(0,1);
	c.d[13] = a(0,1)*a(1,2);
	c.d[18] = a(0,1)*a(0,2);
	c.d[24] = a(0,1)*a(1,0);
	c.d[31] = a(0,1)*a(2,1);
	c.d[39] = a(0,1)*a(2,0);

	c.d[14] = a(1,2)*a(1,2);
	c.d[19] = a(1,2)*a(0,2);
	c.d[25] = a(1,2)*a(1,0);
	c.d[32] = a(1,2)*a(2,1);
	c.d[40] = a(1,2)*a(2,0);

	c.d[20] = a(0,2)*a(0,2);
	c.d[26] = a(0,2)*a(1,0);
	c.d[33] = a(0,2)*a(2,1);
	c.d[41] = a(0,2)*a(2,0);

	c.d[27] = a(1,0)*a(1,0);
	c.d[34] = a(1,0)*a(2,1);
	c.d[42] = a(1,0)*a(2,0);

	c.d[35] = a(2,1)*a(2,1);
	c.d[43] = a(2,1)*a(2,0);

	c.d[44] = a(2,0)*a(2,0);

	return c;
}

//-----------------------------------------------------------------------------
// (a dyad1s b)_ijkl = a_ij b_kl + b_ij a_kl
inline tens4dms dyad1ms(const mat3d& a, const mat3d& b)
{
	tens4dms c;
	
	c.d[0] = 2*a(0,0)*b(0,0);
	c.d[1] = a(0,0)*b(1,1) + b(0,0)*a(1,1);
	c.d[3] = a(0,0)*b(2,2) + b(0,0)*a(2,2);
	c.d[6] = a(0,0)*b(0,1) + b(0,0)*a(0,1);
	c.d[10] = a(0,0)*b(1,2) + b(0,0)*a(1,2);
	c.d[15] = a(0,0)*b(0,2) + b(0,0)*a(0,2);
	c.d[21] = a(0,0)*b(1,0) + b(0,0)*a(1,0);
	c.d[28] = a(0,0)*b(2,1) + b(0,0)*a(2,1);
	c.d[36] = a(0,0)*b(2,0) + b(0,0)*a(2,0);
	
	c.d[2] = 2*a(1,1)*b(1,1);
	c.d[4] = a(1,1)*b(2,2) + b(1,1)*a(2,2);
	c.d[7] = a(1,1)*b(0,1) + b(1,1)*a(0,1);
	c.d[11] = a(1,1)*b(1,2) + b(1,1)*a(1,2);
	c.d[16] = a(1,1)*b(0,2) + b(1,1)*a(0,2);
	c.d[22] = a(1,1)*b(1,0) + b(1,1)*a(1,0);
	c.d[29] = a(1,1)*b(2,1) + b(1,1)*a(2,1);
	c.d[37] = a(1,1)*b(2,0) + b(1,1)*a(2,0);
	
	c.d[5] = 2*a(2,2)*b(2,2);
	c.d[8] = a(2,2)*b(0,1) + b(2,2)*a(0,1);
	c.d[12] = a(2,2)*b(1,2) + b(2,2)*a(1,2);
	c.d[17] = a(2,2)*b(0,2) + b(2,2)*a(0,2);
	c.d[23] = a(2,2)*b(1,0) + b(2,2)*a(1,0);
	c.d[30] = a(2,2)*b(2,1) + b(2,2)*a(2,1);
	c.d[38] = a(2,2)*b(2,0) + b(2,2)*a(2,0);

	c.d[9] = 2*a(0,1)*b(0,1);
	c.d[13] = a(0,1)*b(1,2) + b(0,1)*a(1,2);
	c.d[18] = a(0,1)*b(0,2) + b(0,1)*a(0,2);
	c.d[24] = a(0,1)*b(1,0) + b(0,1)*a(1,0);
	c.d[31] = a(0,1)*b(2,1) + b(0,1)*a(2,1);
	c.d[39] = a(0,1)*b(2,0) + b(0,1)*a(2,0);

	c.d[14] = 2*a(1,2)*b(1,2);
	c.d[19] = a(1,2)*b(0,2) + b(1,2)*a(0,2);
	c.d[25] = a(1,2)*b(1,0) + b(1,2)*a(1,0);
	c.d[32] = a(1,2)*b(2,1) + b(1,2)*a(2,1);
	c.d[40] = a(1,2)*b(2,0) + b(1,2)*a(2,0);

	c.d[20] = 2*a(0,2)*b(0,2);
	c.d[26] = a(0,2)*b(1,0) + b(0,2)*a(1,0);
	c.d[33] = a(0,2)*b(2,1) + b(0,2)*a(2,1);
	c.d[41] = a(0,2)*b(2,0) + b(0,2)*a(2,0);

	c.d[27] = 2*a(1,0)*b(1,0);
	c.d[34] = a(1,0)*b(2,1) + b(1,0)*a(2,1);
	c.d[42] = a(1,0)*b(2,0) + b(1,0)*a(2,0);

	c.d[35] = 2*a(2,1)*b(2,1);
	c.d[43] = a(2,1)*b(2,0) + b(2,1)*a(2,0);

	c.d[44] = 2*a(2,0)*b(2,0);

	return c;
}
