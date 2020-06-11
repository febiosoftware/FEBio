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

inline tens4ds::tens4ds(const double g)
{
	d[ 0] =
	d[ 1] = d[ 2] =
	d[ 3] = d[ 4] = d[ 5] =
	d[ 6] = d[ 7] = d[ 8] = d[ 9] =
	d[10] = d[11] = d[12] = d[13] = d[14] =
	d[15] = d[16] = d[17] = d[18] = d[19] = d[20] = g;
}

inline tens4ds::tens4ds(double m[6][6])
{
	d[ 0] = m[0][0];
	d[ 1] = m[0][1]; d[ 2] = m[1][1];
	d[ 3] = m[0][2]; d[ 4] = m[1][2]; d[ 5] = m[2][2];
	d[ 6] = m[0][3]; d[ 7] = m[1][3]; d[ 8] = m[2][3]; d[ 9] = m[3][3];
	d[10] = m[0][4]; d[11] = m[1][4]; d[12] = m[2][4]; d[13] = m[3][4]; d[14] = m[4][4];
	d[15] = m[0][5]; d[16] = m[1][5]; d[17] = m[2][5]; d[18] = m[3][5]; d[19] = m[4][5]; d[20] = m[5][5];
}

inline double& tens4ds::operator () (int i, int j, int k, int l)
{
	const int m[3][3] = {{0,3,5},{3,1,4},{5,4,2}};
	tens4ds& T = (*this);
	return T(m[i][j], m[k][l]);
}

inline double tens4ds::operator () (int i, int j, int k, int l) const
{
	const int m[3][3] = {{0,3,5},{3,1,4},{5,4,2}};
	const tens4ds& T = (*this);
	return T(m[i][j], m[k][l]);
}

inline double& tens4ds::operator () (int i, int j)
{
	const int m[6] = {0, 1, 3, 6, 10, 15};
	if (i<=j) return d[m[j]+i]; else return d[m[i]+j];
}

inline double tens4ds::operator () (int i, int j) const
{
	const int m[6] = {0, 1, 3, 6, 10, 15};
	if (i<=j) return d[m[j]+i]; else return d[m[i]+j];
}

// operator +
inline tens4ds tens4ds::operator + (const tens4ds& t) const
{
	tens4ds s;
//	for (int i=0; i<NNZ; i++)
//		s.d[i] = d[i] + t.d[i];
	s.d[ 0] = d[ 0] + t.d[ 0];
	s.d[ 1] = d[ 1] + t.d[ 1];
	s.d[ 2] = d[ 2] + t.d[ 2];
	s.d[ 3] = d[ 3] + t.d[ 3];
	s.d[ 4] = d[ 4] + t.d[ 4];
	s.d[ 5] = d[ 5] + t.d[ 5];
	s.d[ 6] = d[ 6] + t.d[ 6];
	s.d[ 7] = d[ 7] + t.d[ 7];
	s.d[ 8] = d[ 8] + t.d[ 8];
	s.d[ 9] = d[ 9] + t.d[ 9];
	s.d[10] = d[10] + t.d[10];
	s.d[11] = d[11] + t.d[11];
	s.d[12] = d[12] + t.d[12];
	s.d[13] = d[13] + t.d[13];
	s.d[14] = d[14] + t.d[14];
	s.d[15] = d[15] + t.d[15];
	s.d[16] = d[16] + t.d[16];
	s.d[17] = d[17] + t.d[17];
	s.d[18] = d[18] + t.d[18];
	s.d[19] = d[19] + t.d[19];
	s.d[20] = d[20] + t.d[20];
	return s;
}

// operator -
inline tens4ds tens4ds::operator - (const tens4ds& t) const
{
	tens4ds s;
//	for (int i=0; i<NNZ; i++)
//		s.d[i] = d[i] - t.d[i];
	s.d[ 0] = d[ 0] - t.d[ 0];
	s.d[ 1] = d[ 1] - t.d[ 1];
	s.d[ 2] = d[ 2] - t.d[ 2];
	s.d[ 3] = d[ 3] - t.d[ 3];
	s.d[ 4] = d[ 4] - t.d[ 4];
	s.d[ 5] = d[ 5] - t.d[ 5];
	s.d[ 6] = d[ 6] - t.d[ 6];
	s.d[ 7] = d[ 7] - t.d[ 7];
	s.d[ 8] = d[ 8] - t.d[ 8];
	s.d[ 9] = d[ 9] - t.d[ 9];
	s.d[10] = d[10] - t.d[10];
	s.d[11] = d[11] - t.d[11];
	s.d[12] = d[12] - t.d[12];
	s.d[13] = d[13] - t.d[13];
	s.d[14] = d[14] - t.d[14];
	s.d[15] = d[15] - t.d[15];
	s.d[16] = d[16] - t.d[16];
	s.d[17] = d[17] - t.d[17];
	s.d[18] = d[18] - t.d[18];
	s.d[19] = d[19] - t.d[19];
	s.d[20] = d[20] - t.d[20];
	return s;
}

// operator *
inline tens4ds tens4ds::operator * (double g) const
{
	tens4ds s;
//	for (int i=0; i<NNZ; i++)
//		s.d[i] = g*d[i];
	s.d[ 0] = g*d[ 0];
	s.d[ 1] = g*d[ 1];
	s.d[ 2] = g*d[ 2];
	s.d[ 3] = g*d[ 3];
	s.d[ 4] = g*d[ 4];
	s.d[ 5] = g*d[ 5];
	s.d[ 6] = g*d[ 6];
	s.d[ 7] = g*d[ 7];
	s.d[ 8] = g*d[ 8];
	s.d[ 9] = g*d[ 9];
	s.d[10] = g*d[10];
	s.d[11] = g*d[11];
	s.d[12] = g*d[12];
	s.d[13] = g*d[13];
	s.d[14] = g*d[14];
	s.d[15] = g*d[15];
	s.d[16] = g*d[16];
	s.d[17] = g*d[17];
	s.d[18] = g*d[18];
	s.d[19] = g*d[19];
	s.d[20] = g*d[20];
	return s;
}

// operator /
inline tens4ds tens4ds::operator / (double g) const
{
	tens4ds s;
//	for (int i=0; i<NNZ; i++)
//		s.d[i] = d[i]/g;
	s.d[ 0] = d[ 0]/g;
	s.d[ 1] = d[ 1]/g;
	s.d[ 2] = d[ 2]/g;
	s.d[ 3] = d[ 3]/g;
	s.d[ 4] = d[ 4]/g;
	s.d[ 5] = d[ 5]/g;
	s.d[ 6] = d[ 6]/g;
	s.d[ 7] = d[ 7]/g;
	s.d[ 8] = d[ 8]/g;
	s.d[ 9] = d[ 9]/g;
	s.d[10] = d[10]/g;
	s.d[11] = d[11]/g;
	s.d[12] = d[12]/g;
	s.d[13] = d[13]/g;
	s.d[14] = d[14]/g;
	s.d[15] = d[15]/g;
	s.d[16] = d[16]/g;
	s.d[17] = d[17]/g;
	s.d[18] = d[18]/g;
	s.d[19] = d[19]/g;
	s.d[20] = d[20]/g;
	return s;
}

// assignment operator +=
inline tens4ds& tens4ds::operator += (const tens4ds& t)
{
//	for (int i=0; i<NNZ; i++)
//		d[i] += t.d[i];
	d[ 0] += t.d[ 0];
	d[ 1] += t.d[ 1];
	d[ 2] += t.d[ 2];
	d[ 3] += t.d[ 3];
	d[ 4] += t.d[ 4];
	d[ 5] += t.d[ 5];
	d[ 6] += t.d[ 6];
	d[ 7] += t.d[ 7];
	d[ 8] += t.d[ 8];
	d[ 9] += t.d[ 9];
	d[10] += t.d[10];
	d[11] += t.d[11];
	d[12] += t.d[12];
	d[13] += t.d[13];
	d[14] += t.d[14];
	d[15] += t.d[15];
	d[16] += t.d[16];
	d[17] += t.d[17];
	d[18] += t.d[18];
	d[19] += t.d[19];
	d[20] += t.d[20];
	return (*this);
}

// assignment operator -=
inline tens4ds& tens4ds::operator -= (const tens4ds& t)
{
//	for (int i=0; i<NNZ; i++)
//		d[i] -= t.d[i];
	d[ 0] -= t.d[ 0];
	d[ 1] -= t.d[ 1];
	d[ 2] -= t.d[ 2];
	d[ 3] -= t.d[ 3];
	d[ 4] -= t.d[ 4];
	d[ 5] -= t.d[ 5];
	d[ 6] -= t.d[ 6];
	d[ 7] -= t.d[ 7];
	d[ 8] -= t.d[ 8];
	d[ 9] -= t.d[ 9];
	d[10] -= t.d[10];
	d[11] -= t.d[11];
	d[12] -= t.d[12];
	d[13] -= t.d[13];
	d[14] -= t.d[14];
	d[15] -= t.d[15];
	d[16] -= t.d[16];
	d[17] -= t.d[17];
	d[18] -= t.d[18];
	d[19] -= t.d[19];
	d[20] -= t.d[20];
	return (*this);
}

// assignment operator *=
inline tens4ds& tens4ds::operator *= (double g)
{
//	for (int i=0; i<NNZ; i++)
//		d[i] *= g;
	d[ 0] *= g;
	d[ 1] *= g;
	d[ 2] *= g;
	d[ 3] *= g;
	d[ 4] *= g;
	d[ 5] *= g;
	d[ 6] *= g;
	d[ 7] *= g;
	d[ 8] *= g;
	d[ 9] *= g;
	d[10] *= g;
	d[11] *= g;
	d[12] *= g;
	d[13] *= g;
	d[14] *= g;
	d[15] *= g;
	d[16] *= g;
	d[17] *= g;
	d[18] *= g;
	d[19] *= g;
	d[20] *= g;
	return (*this);
}

// assignment operator /=
inline tens4ds& tens4ds::operator /= (double g)
{
//	for (int i=0; i<NNZ; i++)
//		d[i] /= g;
	d[ 0] /= g;
	d[ 1] /= g;
	d[ 2] /= g;
	d[ 3] /= g;
	d[ 4] /= g;
	d[ 5] /= g;
	d[ 6] /= g;
	d[ 7] /= g;
	d[ 8] /= g;
	d[ 9] /= g;
	d[10] /= g;
	d[11] /= g;
	d[12] /= g;
	d[13] /= g;
	d[14] /= g;
	d[15] /= g;
	d[16] /= g;
	d[17] /= g;
	d[18] /= g;
	d[19] /= g;
	d[20] /= g;
	return (*this);
}

// unary operator -
inline tens4ds tens4ds::operator - () const
{
	tens4ds s;
	s.d[ 0] = -d[ 0];
	s.d[ 1] = -d[ 1];
	s.d[ 2] = -d[ 2];
	s.d[ 3] = -d[ 3];
	s.d[ 4] = -d[ 4];
	s.d[ 5] = -d[ 5];
	s.d[ 6] = -d[ 6];
	s.d[ 7] = -d[ 7];
	s.d[ 8] = -d[ 8];
	s.d[ 9] = -d[ 9];
	s.d[10] = -d[10];
	s.d[11] = -d[11];
	s.d[12] = -d[12];
	s.d[13] = -d[13];
	s.d[14] = -d[14];
	s.d[15] = -d[15];
	s.d[16] = -d[16];
	s.d[17] = -d[17];
	s.d[18] = -d[18];
	s.d[19] = -d[19];
	s.d[20] = -d[20];
	return s;
}

// trace
// C.tr() = I:C:I
inline double tens4ds::tr() const
{
	return (d[0]+d[2]+d[5]+2*(d[1]+d[3]+d[4]));
}

// intialize to zero
inline void tens4ds::zero()
{
	d[0] = d[1] = d[2] = d[3] = d[4] = d[5] = d[6] = d[7] = d[8] = d[9] =
	d[10] = d[11] = d[12] = d[13] = d[14] = d[15] = d[16] = d[17] = d[18] = d[19] = d[20] = 0;
}

// extract 6x6 matrix
inline void tens4ds::extract(double D[6][6])
{
	D[0][0] = d[0];  D[0][1] = d[1];  D[0][2] = d[3];  D[0][3] = d[6];  D[0][4] = d[10]; D[0][5] = d[15];
	D[1][0] = d[1];  D[1][1] = d[2];  D[1][2] = d[4];  D[1][3] = d[7];  D[1][4] = d[11]; D[1][5] = d[16];
	D[2][0] = d[3];  D[2][1] = d[4];  D[2][2] = d[5];  D[2][3] = d[8];  D[2][4] = d[12]; D[2][5] = d[17];
	D[3][0] = d[6];  D[3][1] = d[7];  D[3][2] = d[8];  D[3][3] = d[9];  D[3][4] = d[13]; D[3][5] = d[18];
	D[4][0] = d[10]; D[4][1] = d[11]; D[4][2] = d[12]; D[4][3] = d[13]; D[4][4] = d[14]; D[4][5] = d[19];
	D[5][0] = d[15]; D[5][1] = d[16]; D[5][2] = d[17]; D[5][3] = d[18]; D[5][4] = d[19]; D[5][5] = d[20];
}

//-----------------------------------------------------------------------------
// (a dyad1s a)_ijkl = a_ij a_kl
inline tens4ds dyad1s(const mat3dd& a)
{
	tens4ds c;
    
    c.d[ 0] = a.xx()*a.xx();

    c.d[ 1] = a.xx()*a.yy();
    c.d[ 2] = a.yy()*a.yy();

    c.d[ 3] = a.xx()*a.zz();
    c.d[ 4] = a.yy()*a.zz();
    c.d[ 5] = a.zz()*a.zz();

    c.d[ 6] = 0.0;
    c.d[ 7] = 0.0;
    c.d[ 8] = 0.0;
    c.d[ 9] = 0.0;

    c.d[10] = 0.0;
    c.d[11] = 0.0;
    c.d[12] = 0.0;
    c.d[13] = 0.0;
    c.d[14] = 0.0;

    c.d[15] = 0.0;
    c.d[16] = 0.0;
    c.d[17] = 0.0;
    c.d[18] = 0.0;
    c.d[19] = 0.0;
    c.d[20] = 0.0;
    
	return c;
}

//-----------------------------------------------------------------------------
// (a dyad1s a)_ijkl = a_ij a_kl
inline tens4ds dyad1s(const mat3ds& a)
{
	tens4ds c;
    
    c.d[ 0] = a.xx()*a.xx();

    c.d[ 1] = a.xx()*a.yy();
    c.d[ 2] = a.yy()*a.yy();

    c.d[ 3] = a.xx()*a.zz();
    c.d[ 4] = a.yy()*a.zz();
    c.d[ 5] = a.zz()*a.zz();

    c.d[ 6] = a.xx()*a.xy();
    c.d[ 7] = a.xy()*a.yy();
    c.d[ 8] = a.xy()*a.zz();
    c.d[ 9] = a.xy()*a.xy();

    c.d[10] = a.xx()*a.yz();
    c.d[11] = a.yy()*a.yz();
    c.d[12] = a.yz()*a.zz();
    c.d[13] = a.xy()*a.yz();
    c.d[14] = a.yz()*a.yz();

    c.d[15] = a.xx()*a.xz();
    c.d[16] = a.xz()*a.yy();
    c.d[17] = a.xz()*a.zz();
    c.d[18] = a.xy()*a.xz();
    c.d[19] = a.xz()*a.yz();
    c.d[20] = a.xz()*a.xz();
    
	return c;
}

//-----------------------------------------------------------------------------
// (a dyad1s b)_ijkl = a_ij b_kl + b_ij a_kl
inline tens4ds dyad1s(const mat3ds& a, const mat3ds& b)
{
	tens4ds c;

    c.d[ 0] = 2*a.xx()*b.xx();
    
    c.d[ 1] = a.xx()*b.yy() + b.xx()*a.yy();
    c.d[ 2] = 2*a.yy()*b.yy();

    c.d[ 3] = a.xx()*b.zz() + b.xx()*a.zz();
    c.d[ 4] = a.yy()*b.zz() + b.yy()*a.zz();
    c.d[ 5] = 2*a.zz()*b.zz();

    c.d[ 6] = a.xx()*b.xy() + b.xx()*a.xy();
    c.d[ 7] = a.xy()*b.yy() + b.xy()*a.yy();
    c.d[ 8] = a.xy()*b.zz() + b.xy()*a.zz();
    c.d[ 9] = 2*a.xy()*b.xy();

    c.d[10] = a.xx()*b.yz() + b.xx()*a.yz();
    c.d[11] = a.yy()*b.yz() + b.yy()*a.yz();
    c.d[12] = a.yz()*b.zz() + b.yz()*a.zz();
    c.d[13] = a.xy()*b.yz() + b.xy()*a.yz();
    c.d[14] = 2*a.yz()*b.yz();

    c.d[15] = a.xx()*b.xz() + b.xx()*a.xz();
    c.d[16] = a.xz()*b.yy() + b.xz()*a.yy();
    c.d[17] = a.xz()*b.zz() + b.xz()*a.zz();
    c.d[18] = a.xy()*b.xz() + b.xy()*a.xz();
    c.d[19] = a.xz()*b.yz() + b.xz()*a.yz();
    c.d[20] = 2*a.xz()*b.xz();

	return c;
}

//-----------------------------------------------------------------------------
// (a dyad1s b)_ijkl = a_ij b_kl + b_ij a_kl
inline tens4ds dyad1s(const mat3dd& a, const mat3dd& b)
{
	tens4ds c;
    
    c.d[ 0] = 2*a.xx()*b.xx();
    
    c.d[ 1] = a.xx()*b.yy() + b.xx()*a.yy();
    c.d[ 2] = 2*a.yy()*b.yy();

    c.d[ 3] = a.xx()*b.zz() + b.xx()*a.zz();
    c.d[ 4] = a.yy()*b.zz() + b.yy()*a.zz();
    c.d[ 5] = 2*a.zz()*b.zz();

    c.d[ 6] = 0.0;
    c.d[ 7] = 0.0;
    c.d[ 8] = 0.0;
    c.d[ 9] = 0.0;

    c.d[10] = 0.0;
    c.d[11] = 0.0;
    c.d[12] = 0.0;
    c.d[13] = 0.0;
    c.d[14] = 0.0;

    c.d[15] = 0.0;
    c.d[16] = 0.0;
    c.d[17] = 0.0;
    c.d[18] = 0.0;
    c.d[19] = 0.0;
    c.d[20] = 0.0;
    
	return c;
}

//-----------------------------------------------------------------------------
// (a dyad1s b)_ijkl = a_ij b_kl + b_ij a_kl
inline tens4ds dyad1s(const mat3ds& a, const mat3dd& b)
{
	tens4ds c;
    
    c.d[ 0] = 2*a.xx()*b.xx();
    
    c.d[ 1] = a.xx()*b.yy() + b.xx()*a.yy();
    c.d[ 2] = 2*a.yy()*b.yy();

    c.d[ 3] = a.xx()*b.zz() + b.xx()*a.zz();
    c.d[ 4] = a.yy()*b.zz() + b.yy()*a.zz();
    c.d[ 5] = 2*a.zz()*b.zz();

    c.d[ 6] = b.xx()*a.xy();
    c.d[ 7] = a.xy()*b.yy();
    c.d[ 8] = a.xy()*b.zz();
    c.d[ 9] = 0.0;

    c.d[10] = b.xx()*a.yz();
    c.d[11] = b.yy()*a.yz();
    c.d[12] = a.yz()*b.zz();
    c.d[13] = 0.0;
    c.d[14] = 0.0;

    c.d[15] = b.xx()*a.xz();
    c.d[16] = a.xz()*b.yy();
    c.d[17] = a.xz()*b.zz();
    c.d[18] = 0.0;
    c.d[19] = 0.0;
    c.d[20] = 0.0;
    
	return c;
}

//-----------------------------------------------------------------------------
// (a dyad4s a)_ijkl = (a_ik a_jl + a_il a_jk)/2
inline tens4ds dyad4s(const mat3dd& a)
{
	tens4ds c;
    
    c.d[ 0] = a.xx()*a.xx();

    c.d[ 1] = 0.0;
    c.d[ 2] = a.yy()*a.yy();

    c.d[ 3] = 0.0;
    c.d[ 4] = 0.0;
    c.d[ 5] = a.zz()*a.zz();

    c.d[ 6] = 0.0;
    c.d[ 7] = 0.0;
    c.d[ 8] = 0.0;
    c.d[ 9] = a.xx()*a.yy()/2;

    c.d[10] = 0.0;
    c.d[11] = 0.0;
    c.d[12] = 0.0;
    c.d[13] = 0.0;
    c.d[14] = a.yy()*a.zz()/2;

    c.d[15] = 0.0;
    c.d[16] = 0.0;
    c.d[17] = 0.0;
    c.d[18] = 0.0;
    c.d[19] = 0.0;
    c.d[20] = a.xx()*a.zz()/2;

	return c;
}

//-----------------------------------------------------------------------------
// (a dyad4s a)_ijkl = (a_ik a_jl + a_il a_jk)/2
inline tens4ds dyad4s(const mat3ds& a)
{
	tens4ds c;
    
    c.d[ 0] = a.xx()*a.xx();

    c.d[ 1] = a.xy()*a.xy();
    c.d[ 2] = a.yy()*a.yy();

    c.d[ 3] = a.xz()*a.xz();
    c.d[ 4] = a.yz()*a.yz();
    c.d[ 5] = a.zz()*a.zz();

    c.d[ 6] = a.xx()*a.xy();
    c.d[ 7] = a.xy()*a.yy();
    c.d[ 8] = a.xz()*a.yz();
    c.d[ 9] = (a.xx()*a.yy() + a.xy()*a.xy())/2;

    c.d[10] = a.xy()*a.xz();
    c.d[11] = a.yy()*a.yz();
    c.d[12] = a.yz()*a.zz();
    c.d[13] = (a.xy()*a.yz() + a.xz()*a.yy())/2;
    c.d[14] = (a.yy()*a.zz() + a.yz()*a.yz())/2;

    c.d[15] = a.xx()*a.xz();
    c.d[16] = a.xy()*a.yz();
    c.d[17] = a.xz()*a.zz();
    c.d[18] = (a.xx()*a.yz() + a.xy()*a.xz())/2;
    c.d[19] = (a.xy()*a.zz() + a.xz()*a.yz())/2;
    c.d[20] = (a.xx()*a.zz() + a.xz()*a.xz())/2;

	return c;
}

//-----------------------------------------------------------------------------
// (a dyad4s b)_ijkl = (a_ik b_jl + a_il b_jk)/2 +  (b_ik a_jl + b_il a_jk)/2
inline tens4ds dyad4s(const mat3ds& a, const mat3dd& b)
{
	tens4ds c;

    c.d[ 0] = 2*a.xx()*b.xx();

    c.d[ 1] = 0.0;
    c.d[ 2] = 2*a.yy()*b.yy();

    c.d[ 3] = 0.0;
    c.d[ 4] = 0.0;
    c.d[ 5] = 2*a.zz()*b.zz();

    c.d[ 6] = b.xx()*a.xy();
    c.d[ 7] = a.xy()*b.yy();
    c.d[ 8] = 0.0;
    c.d[ 9] = (a.xx()*b.yy() + b.xx()*a.yy())/2;

    c.d[10] = 0.0;
    c.d[11] = b.yy()*a.yz();
    c.d[12] = a.yz()*b.zz();
    c.d[13] = a.xz()*b.yy()/2;
    c.d[14] = (a.yy()*b.zz() + b.yy()*a.zz())/2;

    c.d[15] = b.xx()*a.xz();
    c.d[16] = 0.0;
    c.d[17] = a.xz()*b.zz();
    c.d[18] = b.xx()*a.yz()/2;
    c.d[19] = a.xy()*b.zz()/2;
    c.d[20] = (a.xx()*b.zz() + b.xx()*a.zz())/2;
    
	return c;
}

//-----------------------------------------------------------------------------
// (a dyad5s b)_ijkl = (a_ik b_jl + a_il b_jk)/2 +  (a_jl b_ik + a_jk b_il)/2
inline tens4ds dyad5s(const mat3ds& a, const mat3ds& b)
{
	int L[21][4] = { { 0, 0, 0, 0},
	{ 0, 0, 1, 1},{ 1, 1, 1, 1},
	{ 0, 0, 2, 2},{ 1, 1, 2, 2},{ 2, 2, 2, 2},
	{ 0, 0, 0, 1},{ 1, 1, 0, 1},{ 2, 2, 0, 1},{ 0, 1, 0, 1},
	{ 0, 0, 1, 2},{ 1, 1, 1, 2},{ 2, 2, 1, 2},{ 0, 1, 1, 2},{ 1, 2, 1, 2 },
	{ 0, 0, 0, 2},{ 1, 1, 0, 2},{ 2, 2, 0, 2},{ 0, 1, 0, 2},{ 1, 2, 0, 2 },{ 0, 2, 0, 2 }};

	tens4ds c;
	for (int n = 0; n < 21; ++n)
	{
		int i = L[n][0];
		int j = L[n][1];
		int k = L[n][2];
		int l = L[n][3];
		c.d[n] = 0.5*(a(i, k)*b(j, l) + a(i, l)*b(j, k)) + 0.5*(a(j, l)*b(i, k) + a(j, k)*b(i, l));
	}
	return c;
}

//-----------------------------------------------------------------------------
// (a dyad4s b)_ijkl = (a_ik b_jl + a_il b_jk)/2 +  (b_ik a_jl + b_il a_jk)/2
inline tens4ds dyad4s(const mat3ds& a, const mat3ds& b)
{
	tens4ds c;
    
    c.d[ 0] = 2*a.xx()*b.xx();

    c.d[ 1] = 2*a.xy()*b.xy();
    c.d[ 2] = 2*a.yy()*b.yy();

    c.d[ 3] = 2*a.xz()*b.xz();
    c.d[ 4] = 2*a.yz()*b.yz();
    c.d[ 5] = 2*a.zz()*b.zz();

    c.d[ 6] = a.xx()*b.xy() + b.xx()*a.xy();
    c.d[ 7] = a.xy()*b.yy() + b.xy()*a.yy();
    c.d[ 8] = a.xz()*b.yz() + b.xz()*a.yz();
    c.d[ 9] = (a.xx()*b.yy() + 2*a.xy()*b.xy() + b.xx()*a.yy())/2;

    c.d[10] = a.xy()*b.xz() + b.xy()*a.xz();
    c.d[11] = a.yy()*b.yz() + b.yy()*a.yz();
    c.d[12] = a.yz()*b.zz() + b.yz()*a.zz();
    c.d[13] = (a.xy()*b.yz() + a.xz()*b.yy() + b.xy()*a.yz() + b.xz()*a.yy())/2;
    c.d[14] = (a.yy()*b.zz() + 2*a.yz()*b.yz() + b.yy()*a.zz())/2;

    c.d[15] = a.xx()*b.xz() + b.xx()*a.xz();
    c.d[16] = a.xy()*b.yz() + b.xy()*a.yz();
    c.d[17] = a.xz()*b.zz() + b.xz()*a.zz();
    c.d[18] = (a.xx()*b.yz() + a.xy()*b.xz() + b.xx()*a.yz() + b.xy()*a.xz())/2;
    c.d[19] = (a.xy()*b.zz() + a.xz()*b.yz() + b.xy()*a.zz() + b.xz()*a.yz())/2;
    c.d[20] = (a.xx()*b.zz() + 2*a.xz()*b.xz() + b.xx()*a.zz())/2;
    
	return c;

}

//-----------------------------------------------------------------------------
// (a ddots b)_ijkl = a_ijmn b_mnkl + b_ijmn a_mnkl
inline tens4ds ddots(const tens4ds& a, const tens4ds& b)
{
	tens4ds c;
	
	c.d[0] = 2*(a.d[0]*b.d[0] + a.d[1]*b.d[1] + a.d[3]*b.d[3] + 2*a.d[6]*b.d[6] 
				+ 2*a.d[10]*b.d[10] + 2*a.d[15]*b.d[15]);
	c.d[1] = a.d[0]*b.d[1] + a.d[2]*b.d[1] + a.d[1]*(b.d[0] + b.d[2]) + a.d[4]*b.d[3] 
	+ a.d[3]*b.d[4] + 2*a.d[7]*b.d[6] + 2*a.d[6]*b.d[7] + 2*a.d[11]*b.d[10] 
	+ 2*a.d[10]*b.d[11] + 2*a.d[16]*b.d[15] + 2*a.d[15]*b.d[16];
	c.d[3] = a.d[4]*b.d[1] + a.d[0]*b.d[3] + a.d[5]*b.d[3] + a.d[1]*b.d[4] 
	+ a.d[3]*(b.d[0] + b.d[5]) + 2*a.d[8]*b.d[6] + 2*a.d[6]*b.d[8] + 2*a.d[12]*b.d[10] 
	+ 2*a.d[10]*b.d[12] + 2*a.d[17]*b.d[15] + 2*a.d[15]*b.d[17];
	c.d[6] = a.d[7]*b.d[1] + a.d[8]*b.d[3] + a.d[0]*b.d[6] + 2*a.d[9]*b.d[6] + a.d[1]*b.d[7] 
	+ a.d[3]*b.d[8] + a.d[6]*(b.d[0] + 2*b.d[9]) + 2*a.d[13]*b.d[10] + 2*a.d[10]*b.d[13] 
	+ 2*a.d[18]*b.d[15] + 2*a.d[15]*b.d[18];
	c.d[10] = a.d[11]*b.d[1] + a.d[12]*b.d[3] + 2*a.d[13]*b.d[6] + a.d[0]*b.d[10] 
	+ 2*a.d[14]*b.d[10] + a.d[1]*b.d[11] + a.d[3]*b.d[12] + 2*a.d[6]*b.d[13] 
	+ a.d[10]*(b.d[0] + 2*b.d[14]) + 2*a.d[19]*b.d[15] + 2*a.d[15]*b.d[19];
	c.d[15] = a.d[16]*b.d[1] + a.d[17]*b.d[3] + 2*a.d[18]*b.d[6] + 2*a.d[19]*b.d[10] 
	+ a.d[0]*b.d[15] + 2*a.d[20]*b.d[15] + a.d[1]*b.d[16] + a.d[3]*b.d[17] 
	+ 2*a.d[6]*b.d[18] + 2*a.d[10]*b.d[19] + a.d[15]*(b.d[0] + 2*b.d[20]);
	
	c.d[2] = 2*(a.d[1]*b.d[1] + a.d[2]*b.d[2] + a.d[4]*b.d[4] + 2*a.d[7]*b.d[7] 
				+ 2*a.d[11]*b.d[11] + 2*a.d[16]*b.d[16]);
	c.d[4] = a.d[3]*b.d[1] + a.d[1]*b.d[3] + a.d[2]*b.d[4] + a.d[5]*b.d[4] 
	+ a.d[4]*(b.d[2] + b.d[5]) + 2*a.d[8]*b.d[7] + 2*a.d[7]*b.d[8] + 2*a.d[12]*b.d[11] 
	+ 2*a.d[11]*b.d[12] + 2*a.d[17]*b.d[16] + 2*a.d[16]*b.d[17];
	c.d[7] = a.d[6]*b.d[1] + a.d[8]*b.d[4] + a.d[1]*b.d[6] + a.d[2]*b.d[7] 
	+ 2*a.d[9]*b.d[7] + a.d[4]*b.d[8] + a.d[7]*(b.d[2] + 2*b.d[9]) + 2*a.d[13]*b.d[11] 
	+ 2*a.d[11]*b.d[13] + 2*a.d[18]*b.d[16] + 2*a.d[16]*b.d[18];
	c.d[11] = a.d[10]*b.d[1] + a.d[12]*b.d[4] + 2*a.d[13]*b.d[7] + a.d[1]*b.d[10] 
	+ a.d[2]*b.d[11] + 2*a.d[14]*b.d[11] + a.d[4]*b.d[12] + 2*a.d[7]*b.d[13] 
	+ a.d[11]*(b.d[2] + 2*b.d[14]) + 2*a.d[19]*b.d[16] + 2*a.d[16]*b.d[19];
	c.d[16] = a.d[15]*b.d[1] + a.d[17]*b.d[4] + 2*a.d[18]*b.d[7] + 2*a.d[19]*b.d[11] 
	+ a.d[1]*b.d[15] + a.d[2]*b.d[16] + 2*a.d[20]*b.d[16] + a.d[4]*b.d[17] 
	+ 2*a.d[7]*b.d[18] + 2*a.d[11]*b.d[19] + a.d[16]*(b.d[2] + 2*b.d[20]);
	
	c.d[5] = 2*(a.d[3]*b.d[3] + a.d[4]*b.d[4] + a.d[5]*b.d[5] + 2*a.d[8]*b.d[8] 
				+ 2*a.d[12]*b.d[12] + 2*a.d[17]*b.d[17]);
	c.d[8] = a.d[6]*b.d[3] + a.d[7]*b.d[4] + a.d[8]*b.d[5] + a.d[3]*b.d[6] 
	+ a.d[4]*b.d[7] + a.d[5]*b.d[8] + 2*a.d[9]*b.d[8] + 2*a.d[8]*b.d[9] + 2*a.d[13]*b.d[12] 
	+ 2*a.d[12]*b.d[13] + 2*a.d[18]*b.d[17] + 2*a.d[17]*b.d[18];
	c.d[12] = a.d[10]*b.d[3] + a.d[11]*b.d[4] + a.d[12]*b.d[5] + 2*a.d[13]*b.d[8] 
	+ a.d[3]*b.d[10] + a.d[4]*b.d[11] + a.d[5]*b.d[12] + 2*a.d[14]*b.d[12] + 2*a.d[8]*b.d[13] 
	+ 2*a.d[12]*b.d[14] + 2*a.d[19]*b.d[17] + 2*a.d[17]*b.d[19];
	c.d[17] = a.d[15]*b.d[3] + a.d[16]*b.d[4] + a.d[17]*b.d[5] + 2*a.d[18]*b.d[8] 
	+ 2*a.d[19]*b.d[12] + a.d[3]*b.d[15] + a.d[4]*b.d[16] + a.d[5]*b.d[17] 
	+ 2*a.d[20]*b.d[17] + 2*a.d[8]*b.d[18] + 2*a.d[12]*b.d[19] + 2*a.d[17]*b.d[20];
	
	c.d[9] = 2*(a.d[6]*b.d[6] + a.d[7]*b.d[7] + a.d[8]*b.d[8] + 2*a.d[9]*b.d[9] 
				+ 2*a.d[13]*b.d[13] + 2*a.d[18]*b.d[18]);
	c.d[13] = a.d[10]*b.d[6] + a.d[11]*b.d[7] + a.d[12]*b.d[8] + 2*a.d[13]*b.d[9] 
	+ a.d[6]*b.d[10] + a.d[7]*b.d[11] + a.d[8]*b.d[12] + 2*a.d[9]*b.d[13] 
	+ 2*a.d[14]*b.d[13] + 2*a.d[13]*b.d[14] + 2*a.d[19]*b.d[18] + 2*a.d[18]*b.d[19];
	c.d[18] = a.d[15]*b.d[6] + a.d[16]*b.d[7] + a.d[17]*b.d[8] + 2*a.d[18]*b.d[9] 
	+ 2*a.d[19]*b.d[13] + a.d[6]*b.d[15] + a.d[7]*b.d[16] + a.d[8]*b.d[17] 
	+ 2*a.d[9]*b.d[18] + 2*a.d[20]*b.d[18] + 2*a.d[13]*b.d[19] + 2*a.d[18]*b.d[20];
	
	c.d[14] = 2*(a.d[10]*b.d[10] + a.d[11]*b.d[11] + a.d[12]*b.d[12] 
				 + 2*a.d[13]*b.d[13] + 2*a.d[14]*b.d[14] + 2*a.d[19]*b.d[19]);
	c.d[19] = a.d[15]*b.d[10] + a.d[16]*b.d[11] + a.d[17]*b.d[12] 
	+ 2*a.d[18]*b.d[13] + 2*a.d[19]*b.d[14] + a.d[10]*b.d[15] + a.d[11]*b.d[16] 
	+ a.d[12]*b.d[17] + 2*a.d[13]*b.d[18] + 2*a.d[14]*b.d[19] + 2*a.d[20]*b.d[19] + 2*a.d[19]*b.d[20];
	
	c.d[20] = 2*(a.d[15]*b.d[15] + a.d[16]*b.d[16] + a.d[17]*b.d[17] 
				 + 2*a.d[18]*b.d[18] + 2*a.d[19]*b.d[19] + 2*a.d[20]*b.d[20]);
	
	return c;
	
}

//-----------------------------------------------------------------------------
// Evaluates the dyadic product C_ijkl = 0.25*(a_i * K_jk * b_l + perms.)
inline tens4ds dyad4s(const vec3d& a, const mat3d& K, const vec3d& b)
{
	tens4ds c;

	c(0,0) = a.x * K(0,0) * b.x;
	c(1,1) = a.y * K(1,1) * b.y;
	c(2,2) = a.z * K(2,2) * b.z;

	c(0,1) = 0.5*(a.x * K(0,1) * b.y + a.y * K(1,0) * b.x);
	c(0,2) = 0.5*(a.x * K(0,2) * b.z + a.z * K(2,0) * b.x);
	c(1,2) = 0.5*(a.y * K(1,2) * b.z + a.z * K(2,1) * b.y);

	c(0,3) = 0.25*(a.x * K(0,0) * b.y + a.x * K(0,1) * b.x + a.x * K(1,0) * b.x + a.y * K(0,0) * b.x);
	c(0,4) = 0.25*(a.x * K(0,1) * b.z + a.x * K(0,2) * b.y + a.y * K(2,0) * b.x + a.z * K(1,0) * b.x);
	c(0,5) = 0.25*(a.x * K(0,0) * b.z + a.x * K(0,2) * b.x + a.x * K(2,0) * b.x + a.z * K(0,0) * b.x);
	c(1,3) = 0.25*(a.y * K(1,0) * b.y + a.y * K(1,1) * b.x + a.x * K(1,1) * b.y + a.y * K(0,1) * b.y);
	c(1,4) = 0.25*(a.y * K(1,1) * b.z + a.y * K(1,2) * b.y + a.y * K(2,1) * b.y + a.z * K(1,1) * b.y);
	c(1,5) = 0.25*(a.y * K(1,0) * b.z + a.y * K(1,2) * b.x + a.x * K(2,1) * b.y + a.z * K(0,1) * b.y);
	c(2,3) = 0.25*(a.z * K(2,0) * b.y + a.z * K(2,1) * b.x + a.x * K(1,2) * b.z + a.y * K(0,2) * b.z);
	c(2,4) = 0.25*(a.z * K(2,1) * b.z + a.z * K(2,2) * b.y + a.y * K(2,2) * b.z + a.z * K(1,2) * b.z);
	c(2,5) = 0.25*(a.z * K(2,0) * b.z + a.z * K(2,2) * b.x + a.x * K(2,2) * b.z + a.z * K(0,2) * b.z);

	c(3,3) = 0.25*(a.x * K(1,0) * b.y + a.x * K(1,1) * b.x + a.x * K(1,0) * b.y + a.y * K(0,0) * b.y);
	c(3,4) = 0.25*(a.x * K(1,1) * b.z + a.x * K(1,2) * b.y + a.y * K(2,0) * b.y + a.z * K(1,0) * b.y);
	c(3,5) = 0.25*(a.x * K(1,0) * b.z + a.x * K(1,2) * b.x + a.x * K(2,0) * b.y + a.z * K(0,0) * b.y);

	c(4,4) = 0.25*(a.y * K(2,1) * b.z + a.y * K(2,2) * b.y + a.y * K(2,1) * b.z + a.z * K(1,1) * b.z);
	c(4,5) = 0.25*(a.y * K(2,0) * b.z + a.y * K(2,2) * b.x + a.x * K(2,1) * b.z + a.z * K(0,1) * b.z);

	c(5,5) = 0.25*(a.x * K(2,0) * b.z + a.x * K(2,2) * b.x + a.x * K(2,0) * b.z + a.z * K(0,0) * b.z);

	return c;
}

//-----------------------------------------------------------------------------
// double contraction of symmetric 4th-order tensor with a symmetric 2nd-order tensor
// Aij = Dijkl Mkl
inline mat3ds tens4ds::dot(const mat3ds &m) const
{
	mat3ds a;
	a.xx() = d[ 0]*m.xx() + d[ 1]*m.yy() + d[ 3]*m.zz() + 2*d[ 6]*m.xy() + 2*d[10]*m.yz() + 2*d[15]*m.xz();
	a.yy() = d[ 1]*m.xx() + d[ 2]*m.yy() + d[ 4]*m.zz() + 2*d[ 7]*m.xy() + 2*d[11]*m.yz() + 2*d[16]*m.xz();
	a.zz() = d[ 3]*m.xx() + d[ 4]*m.yy() + d[ 5]*m.zz() + 2*d[ 8]*m.xy() + 2*d[12]*m.yz() + 2*d[17]*m.xz();
	a.xy() = d[ 6]*m.xx() + d[ 7]*m.yy() + d[ 8]*m.zz() + 2*d[ 9]*m.xy() + 2*d[13]*m.yz() + 2*d[18]*m.xz();
	a.yz() = d[10]*m.xx() + d[11]*m.yy() + d[12]*m.zz() + 2*d[13]*m.xy() + 2*d[14]*m.yz() + 2*d[19]*m.xz();
	a.xz() = d[15]*m.xx() + d[16]*m.yy() + d[17]*m.zz() + 2*d[18]*m.xy() + 2*d[19]*m.yz() + 2*d[20]*m.xz();
	return a;
}

//-----------------------------------------------------------------------------
// double contraction of symmetric 4th-order tensor with a general 2nd-order tensor
// Aij = Dijkl Mkl
inline mat3ds tens4ds::dot(const mat3d &m) const
{
    mat3ds a;
    a.xx() = d[ 0]*m(0,0) + d[ 1]*m(1,1) + d[ 3]*m(2,2) + d[ 6]*(m(0,1)+m(1,0)) + d[10]*(m(1,2)+m(2,1)) + d[15]*(m(0,2)+m(2,0));
    a.yy() = d[ 1]*m(0,0) + d[ 2]*m(1,1) + d[ 4]*m(2,2) + d[ 7]*(m(0,1)+m(1,0)) + d[11]*(m(1,2)+m(2,1)) + d[16]*(m(0,2)+m(2,0));
    a.zz() = d[ 3]*m(0,0) + d[ 4]*m(1,1) + d[ 5]*m(2,2) + d[ 8]*(m(0,1)+m(1,0)) + d[12]*(m(1,2)+m(2,1)) + d[17]*(m(0,2)+m(2,0));
    a.xy() = d[ 6]*m(0,0) + d[ 7]*m(1,1) + d[ 8]*m(2,2) + d[ 9]*(m(0,1)+m(1,0)) + d[13]*(m(1,2)+m(2,1)) + d[18]*(m(0,2)+m(2,0));
    a.yz() = d[10]*m(0,0) + d[11]*m(1,1) + d[12]*m(2,2) + d[13]*(m(0,1)+m(1,0)) + d[14]*(m(1,2)+m(2,1)) + d[19]*(m(0,2)+m(2,0));
    a.xz() = d[15]*m(0,0) + d[16]*m(1,1) + d[17]*m(2,2) + d[18]*(m(0,1)+m(1,0)) + d[19]*(m(1,2)+m(2,1)) + d[20]*(m(0,2)+m(2,0));
    return a;
}

//-----------------------------------------------------------------------------
// contraction of symmetric 4th-order tensor with a vector
// Aijk = Dijkl Ml
//
//     / 0   1   3   6   10   15  \   / C0000  C0011  C0022  C0001  C0012  C0002 \
//     |     2   4   7   11   16  |   |        C1111  C1122  C1101  C1112  C1102 |
//     |         5   8   12   17  |   |               C2222  C2201  C2212  C2202 |
// A = |             9   13   18  | = |                      C0101  C0112  C0102 |
//     |                 14   19  |   |                             C1212  C1202 |
//     \                      20  /   \                                    C0202 /
//
// [G] = [G111 G112 G113 G121 G122 G123 G131 G132 G133 G221 G222 G223 G231 G232 G233 G331 G332 G333]
//     =    G0   G1   G2   G3   G4   G5   G6   G7   G8   G9  G10  G11  G12  G13  G14  G15  G16  G17
inline tens3dls tens4ds::dot(const vec3d &m) const
{
    tens3dls a;
    a.d[0] = d[0]*m.x + d[6]*m.y + d[15]*m.z;
    a.d[1] = d[6]*m.x + d[1]*m.y + d[10]*m.z;
    a.d[2] = d[15]*m.x + d[10]*m.y + d[3]*m.z;
    a.d[3] = d[6]*m.x + d[9]*m.y + d[18]*m.z;
    a.d[4] = d[9]*m.x + d[7]*m.y + d[13]*m.z;
    a.d[5] = d[18]*m.x + d[13]*m.y + d[8]*m.z;
    a.d[6] = d[15]*m.x + d[18]*m.y + d[20]*m.z;
    a.d[7] = d[18]*m.x + d[16]*m.y + d[19]*m.z;
    a.d[8] = d[20]*m.x + d[19]*m.y + d[17]*m.z;
    a.d[9] = d[1]*m.x + d[7]*m.y + d[16]*m.z;
    a.d[10] = d[7]*m.x + d[2]*m.y + d[11]*m.z;
    a.d[11] = d[16]*m.x + d[11]*m.y + d[4]*m.z;
    a.d[12] = d[10]*m.x + d[13]*m.y + d[19]*m.z;
    a.d[13] = d[13]*m.x + d[11]*m.y + d[14]*m.z;
    a.d[14] = d[19]*m.x + d[14]*m.y + d[12]*m.z;
    a.d[15] = d[3]*m.x + d[8]*m.y + d[17]*m.z;
    a.d[16] = d[8]*m.x + d[4]*m.y + d[12]*m.z;
    a.d[17] = d[17]*m.x + d[12]*m.y + d[5]*m.z;
    return a;
}

//-----------------------------------------------------------------------------
// double contraction of symmetric 4th-order tensor with a general 2nd-order tensor (2nd kind)
// Aij = Dikjl Mkl
inline mat3ds tens4ds::dot2(const mat3d &m) const
{
    mat3ds a;
    a.xx() = d[ 0]*m(0,0) + d[ 9]*m(1,1) + d[20]*m(2,2) + d[ 6]*(m(0,1)+m(1,0)) + d[18]*(m(1,2)+m(2,1)) + d[15]*(m(0,2)+m(2,0));
    a.yy() = d[ 9]*m(0,0) + d[ 2]*m(1,1) + d[14]*m(2,2) + d[ 7]*(m(0,1)+m(1,0)) + d[11]*(m(1,2)+m(2,1)) + d[13]*(m(0,2)+m(2,0));
    a.zz() = d[20]*m(0,0) + d[14]*m(1,1) + d[ 5]*m(2,2) + d[19]*(m(0,1)+m(1,0)) + d[12]*(m(1,2)+m(2,1)) + d[17]*(m(0,2)+m(2,0));
    a.xy() = d[ 6]*m(0,0) + d[ 7]*m(1,1) + d[19]*m(2,2) + d[ 1]*(m(0,1)+m(1,0)) + d[13]*(m(1,2)+m(2,1)) + d[10]*(m(0,2)+m(2,0));
    a.yz() = d[18]*m(0,0) + d[11]*m(1,1) + d[12]*m(2,2) + d[13]*(m(0,1)+m(1,0)) + d[ 4]*(m(1,2)+m(2,1)) + d[ 8]*(m(0,2)+m(2,0));
    a.xz() = d[15]*m(0,0) + d[13]*m(1,1) + d[17]*m(2,2) + d[10]*(m(0,1)+m(1,0)) + d[ 8]*(m(1,2)+m(2,1)) + d[ 3]*(m(0,2)+m(2,0));
    return a;
}

// contraction
// Aij = c_qqij = c_ijqq
inline mat3ds tens4ds::contract() const
{
	mat3ds a;
	a.xx() = d[ 0] + d[ 1] + d[ 3];
	a.yy() = d[ 1] + d[ 2] + d[ 4];
	a.zz() = d[ 3] + d[ 4] + d[ 5];
	a.xy() = d[ 6] + d[ 7] + d[ 8];
	a.yz() = d[10] + d[11] + d[12];
	a.xz() = d[15] + d[16] + d[17];
	return a;
}

//-----------------------------------------------------------------------------
// vdotTdotv_jk = a_i T_ijkl b_l
inline mat3d vdotTdotv(const vec3d& a, const tens4ds& T, const vec3d& b)
{
	return mat3d(a.x*(b.x*T.d[0] + b.y*T.d[6] + b.z*T.d[15]) + a.y*(b.x*T.d[6] + b.y*T.d[9] + b.z*T.d[18]) + a.z*(b.x*T.d[15] + b.y*T.d[18] + b.z*T.d[20]),
				 a.x*(b.y*T.d[1] + b.x*T.d[6] + b.z*T.d[10]) + a.y*(b.y*T.d[7] + b.x*T.d[9] + b.z*T.d[13]) + a.z*(b.y*T.d[16] + b.x*T.d[18] + b.z*T.d[19]),
				 a.x*(b.z*T.d[3] + b.y*T.d[10] + b.x*T.d[15]) + a.y*(b.z*T.d[8] + b.y*T.d[13] + b.x*T.d[18]) + a.z*(b.z*T.d[17] + b.y*T.d[19] + b.x*T.d[20]),
				 a.y*(b.x*T.d[1] + b.y*T.d[7] + b.z*T.d[16]) + a.x*(b.x*T.d[6] + b.y*T.d[9] + b.z*T.d[18]) + a.z*(b.x*T.d[10] + b.y*T.d[13] + b.z*T.d[19]),
				 a.y*(b.y*T.d[2] + b.x*T.d[7] + b.z*T.d[11]) + a.x*(b.y*T.d[7] + b.x*T.d[9] + b.z*T.d[13]) + a.z*(b.y*T.d[11] + b.x*T.d[13] + b.z*T.d[14]),
				 a.y*(b.z*T.d[4] + b.y*T.d[11] + b.x*T.d[16]) + a.x*(b.z*T.d[8] + b.y*T.d[13] + b.x*T.d[18]) + a.z*(b.z*T.d[12] + b.y*T.d[14] + b.x*T.d[19]),
				 a.z*(b.x*T.d[3] + b.y*T.d[8] + b.z*T.d[17]) + a.y*(b.x*T.d[10] + b.y*T.d[13] + b.z*T.d[19]) + a.x*(b.x*T.d[15] + b.y*T.d[18] + b.z*T.d[20]),
				 a.z*(b.y*T.d[4] + b.x*T.d[8] + b.z*T.d[12]) + a.y*(b.y*T.d[11] + b.x*T.d[13] + b.z*T.d[14]) + a.x*(b.y*T.d[16] + b.x*T.d[18] + b.z*T.d[19]),
				 a.z*(b.z*T.d[5] + b.y*T.d[12] + b.x*T.d[17]) + a.y*(b.z*T.d[12] + b.y*T.d[14] + b.x*T.d[19]) + a.x*(b.z*T.d[17] + b.y*T.d[19] + b.x*T.d[20]));
}


//-----------------------------------------------------------------------------
// inverse
inline tens4ds tens4ds::inverse() const
{
	matrix c(6,6);
	
	// populate c
	c(0,0) = d[0];
	c(1,1) = d[2];
	c(2,2) = d[5];
	c(3,3) = d[9];
	c(4,4) = d[14];
	c(5,5) = d[20];
	c(0,1) = c(1,0) = d[1];
	c(0,2) = c(2,0) = d[3];
	c(0,3) = c(3,0) = d[6];
	c(0,4) = c(4,0) = d[10];
	c(0,5) = c(5,0) = d[15];
	c(1,2) = c(2,1) = d[4];
	c(1,3) = c(3,1) = d[7];
	c(1,4) = c(4,1) = d[11];
	c(1,5) = c(5,1) = d[16];
	c(2,3) = c(3,2) = d[8];
	c(2,4) = c(4,2) = d[12];
	c(2,5) = c(5,2) = d[17];
	c(3,4) = c(4,3) = d[13];
	c(3,5) = c(5,3) = d[18];
	c(4,5) = c(5,4) = d[19];
	
	// invert c
	matrix s = c.inverse();
	
	// return inverse
	tens4ds S;
	S.d[ 0] = s(0,0);
	S.d[ 2] = s(1,1);
	S.d[ 5] = s(2,2);
	S.d[ 9] = s(3,3);
	S.d[14] = s(4,4);
	S.d[20] = s(5,5);
	S.d[ 1] = s(0,1);
	S.d[ 3] = s(0,2);
	S.d[ 6] = s(0,3);
	S.d[10] = s(0,4);
	S.d[15] = s(0,5);
	S.d[ 4] = s(1,2);
	S.d[ 7] = s(1,3);
	S.d[11] = s(1,4);
	S.d[16] = s(1,5);
	S.d[ 8] = s(2,3);
	S.d[12] = s(2,4);
	S.d[17] = s(2,5);
	S.d[13] = s(3,4);
	S.d[18] = s(3,5);
	S.d[19] = s(4,5);
	
	return S;
}

//-----------------------------------------------------------------------------
// evaluate push/pull operation
// c_ijpq = F_ik F_jl C_klmn F_pm F_qn
inline tens4ds tens4ds::pp(const mat3d& F)
{
    tens4ds c;

    c.d[0] = d[0]*pow(F(0,0),4) + 4*d[6]*pow(F(0,0),3)*F(0,1) +
    2*d[1]*pow(F(0,0),2)*pow(F(0,1),2) + 4*d[9]*pow(F(0,0),2)*pow(F(0,1),2) +
    4*d[7]*F(0,0)*pow(F(0,1),3) + d[2]*pow(F(0,1),4) +
    4*d[15]*pow(F(0,0),3)*F(0,2) + 4*d[10]*pow(F(0,0),2)*F(0,1)*F(0,2) +
    8*d[18]*pow(F(0,0),2)*F(0,1)*F(0,2) + 8*d[13]*F(0,0)*pow(F(0,1),2)*F(0,2) +
    4*d[16]*F(0,0)*pow(F(0,1),2)*F(0,2) + 4*d[11]*pow(F(0,1),3)*F(0,2) +
    2*d[3]*pow(F(0,0),2)*pow(F(0,2),2) + 4*d[20]*pow(F(0,0),2)*pow(F(0,2),2) +
    4*d[8]*F(0,0)*F(0,1)*pow(F(0,2),2) + 8*d[19]*F(0,0)*F(0,1)*pow(F(0,2),2) +
    2*d[4]*pow(F(0,1),2)*pow(F(0,2),2) + 4*d[14]*pow(F(0,1),2)*pow(F(0,2),2) +
    4*d[17]*F(0,0)*pow(F(0,2),3) + 4*d[12]*F(0,1)*pow(F(0,2),3) +
    d[5]*pow(F(0,2),4);
    
    c.d[1] = d[0]*pow(F(0,0),2)*pow(F(1,0),2) + d[1]*pow(F(0,1),2)*pow(F(1,0),2) +
    2*d[15]*F(0,0)*F(0,2)*pow(F(1,0),2) + 2*d[10]*F(0,1)*F(0,2)*pow(F(1,0),2) +
    d[3]*pow(F(0,2),2)*pow(F(1,0),2) + 4*d[9]*F(0,0)*F(0,1)*F(1,0)*F(1,1) +
    2*d[7]*pow(F(0,1),2)*F(1,0)*F(1,1) + 4*d[18]*F(0,0)*F(0,2)*F(1,0)*F(1,1) +
    4*d[13]*F(0,1)*F(0,2)*F(1,0)*F(1,1) + 2*d[8]*pow(F(0,2),2)*F(1,0)*F(1,1) +
    d[1]*pow(F(0,0),2)*pow(F(1,1),2) + 2*d[7]*F(0,0)*F(0,1)*pow(F(1,1),2) +
    d[2]*pow(F(0,1),2)*pow(F(1,1),2) + 2*d[16]*F(0,0)*F(0,2)*pow(F(1,1),2) +
    2*d[11]*F(0,1)*F(0,2)*pow(F(1,1),2) + d[4]*pow(F(0,2),2)*pow(F(1,1),2) +
    2*d[6]*F(0,0)*F(1,0)*(F(0,1)*F(1,0) + F(0,0)*F(1,1)) +
    2*d[15]*pow(F(0,0),2)*F(1,0)*F(1,2) + 4*d[18]*F(0,0)*F(0,1)*F(1,0)*F(1,2) +
    2*d[16]*pow(F(0,1),2)*F(1,0)*F(1,2) + 4*d[20]*F(0,0)*F(0,2)*F(1,0)*F(1,2) +
    4*d[19]*F(0,1)*F(0,2)*F(1,0)*F(1,2) + 2*d[17]*pow(F(0,2),2)*F(1,0)*F(1,2) +
    2*d[10]*pow(F(0,0),2)*F(1,1)*F(1,2) + 4*d[13]*F(0,0)*F(0,1)*F(1,1)*F(1,2) +
    2*d[11]*pow(F(0,1),2)*F(1,1)*F(1,2) + 4*d[19]*F(0,0)*F(0,2)*F(1,1)*F(1,2) +
    4*d[14]*F(0,1)*F(0,2)*F(1,1)*F(1,2) + 2*d[12]*pow(F(0,2),2)*F(1,1)*F(1,2) +
    d[3]*pow(F(0,0),2)*pow(F(1,2),2) + 2*d[8]*F(0,0)*F(0,1)*pow(F(1,2),2) +
    d[4]*pow(F(0,1),2)*pow(F(1,2),2) + 2*d[17]*F(0,0)*F(0,2)*pow(F(1,2),2) +
    2*d[12]*F(0,1)*F(0,2)*pow(F(1,2),2) + d[5]*pow(F(0,2),2)*pow(F(1,2),2);
    
    c.d[2] = d[0]*pow(F(1,0),4) + 4*d[6]*pow(F(1,0),3)*F(1,1) +
    2*d[1]*pow(F(1,0),2)*pow(F(1,1),2) + 4*d[9]*pow(F(1,0),2)*pow(F(1,1),2) +
    4*d[7]*F(1,0)*pow(F(1,1),3) + d[2]*pow(F(1,1),4) +
    4*d[15]*pow(F(1,0),3)*F(1,2) + 4*d[10]*pow(F(1,0),2)*F(1,1)*F(1,2) +
    8*d[18]*pow(F(1,0),2)*F(1,1)*F(1,2) + 8*d[13]*F(1,0)*pow(F(1,1),2)*F(1,2) +
    4*d[16]*F(1,0)*pow(F(1,1),2)*F(1,2) + 4*d[11]*pow(F(1,1),3)*F(1,2) +
    2*d[3]*pow(F(1,0),2)*pow(F(1,2),2) + 4*d[20]*pow(F(1,0),2)*pow(F(1,2),2) +
    4*d[8]*F(1,0)*F(1,1)*pow(F(1,2),2) + 8*d[19]*F(1,0)*F(1,1)*pow(F(1,2),2) +
    2*d[4]*pow(F(1,1),2)*pow(F(1,2),2) + 4*d[14]*pow(F(1,1),2)*pow(F(1,2),2) +
    4*d[17]*F(1,0)*pow(F(1,2),3) + 4*d[12]*F(1,1)*pow(F(1,2),3) +
    d[5]*pow(F(1,2),4);
    
    c.d[3] = d[0]*pow(F(0,0),2)*pow(F(2,0),2) + d[1]*pow(F(0,1),2)*pow(F(2,0),2) +
    2*d[15]*F(0,0)*F(0,2)*pow(F(2,0),2) + 2*d[10]*F(0,1)*F(0,2)*pow(F(2,0),2) +
    d[3]*pow(F(0,2),2)*pow(F(2,0),2) + 4*d[9]*F(0,0)*F(0,1)*F(2,0)*F(2,1) +
    2*d[7]*pow(F(0,1),2)*F(2,0)*F(2,1) + 4*d[18]*F(0,0)*F(0,2)*F(2,0)*F(2,1) +
    4*d[13]*F(0,1)*F(0,2)*F(2,0)*F(2,1) + 2*d[8]*pow(F(0,2),2)*F(2,0)*F(2,1) +
    d[1]*pow(F(0,0),2)*pow(F(2,1),2) + 2*d[7]*F(0,0)*F(0,1)*pow(F(2,1),2) +
    d[2]*pow(F(0,1),2)*pow(F(2,1),2) + 2*d[16]*F(0,0)*F(0,2)*pow(F(2,1),2) +
    2*d[11]*F(0,1)*F(0,2)*pow(F(2,1),2) + d[4]*pow(F(0,2),2)*pow(F(2,1),2) +
    2*d[6]*F(0,0)*F(2,0)*(F(0,1)*F(2,0) + F(0,0)*F(2,1)) +
    2*d[15]*pow(F(0,0),2)*F(2,0)*F(2,2) + 4*d[18]*F(0,0)*F(0,1)*F(2,0)*F(2,2) +
    2*d[16]*pow(F(0,1),2)*F(2,0)*F(2,2) + 4*d[20]*F(0,0)*F(0,2)*F(2,0)*F(2,2) +
    4*d[19]*F(0,1)*F(0,2)*F(2,0)*F(2,2) + 2*d[17]*pow(F(0,2),2)*F(2,0)*F(2,2) +
    2*d[10]*pow(F(0,0),2)*F(2,1)*F(2,2) + 4*d[13]*F(0,0)*F(0,1)*F(2,1)*F(2,2) +
    2*d[11]*pow(F(0,1),2)*F(2,1)*F(2,2) + 4*d[19]*F(0,0)*F(0,2)*F(2,1)*F(2,2) +
    4*d[14]*F(0,1)*F(0,2)*F(2,1)*F(2,2) + 2*d[12]*pow(F(0,2),2)*F(2,1)*F(2,2) +
    d[3]*pow(F(0,0),2)*pow(F(2,2),2) + 2*d[8]*F(0,0)*F(0,1)*pow(F(2,2),2) +
    d[4]*pow(F(0,1),2)*pow(F(2,2),2) + 2*d[17]*F(0,0)*F(0,2)*pow(F(2,2),2) +
    2*d[12]*F(0,1)*F(0,2)*pow(F(2,2),2) + d[5]*pow(F(0,2),2)*pow(F(2,2),2);
    
    c.d[4] = d[0]*pow(F(1,0),2)*pow(F(2,0),2) + d[1]*pow(F(1,1),2)*pow(F(2,0),2) +
    2*d[15]*F(1,0)*F(1,2)*pow(F(2,0),2) + 2*d[10]*F(1,1)*F(1,2)*pow(F(2,0),2) +
    d[3]*pow(F(1,2),2)*pow(F(2,0),2) + 4*d[9]*F(1,0)*F(1,1)*F(2,0)*F(2,1) +
    2*d[7]*pow(F(1,1),2)*F(2,0)*F(2,1) + 4*d[18]*F(1,0)*F(1,2)*F(2,0)*F(2,1) +
    4*d[13]*F(1,1)*F(1,2)*F(2,0)*F(2,1) + 2*d[8]*pow(F(1,2),2)*F(2,0)*F(2,1) +
    d[1]*pow(F(1,0),2)*pow(F(2,1),2) + 2*d[7]*F(1,0)*F(1,1)*pow(F(2,1),2) +
    d[2]*pow(F(1,1),2)*pow(F(2,1),2) + 2*d[16]*F(1,0)*F(1,2)*pow(F(2,1),2) +
    2*d[11]*F(1,1)*F(1,2)*pow(F(2,1),2) + d[4]*pow(F(1,2),2)*pow(F(2,1),2) +
    2*d[6]*F(1,0)*F(2,0)*(F(1,1)*F(2,0) + F(1,0)*F(2,1)) +
    2*d[15]*pow(F(1,0),2)*F(2,0)*F(2,2) + 4*d[18]*F(1,0)*F(1,1)*F(2,0)*F(2,2) +
    2*d[16]*pow(F(1,1),2)*F(2,0)*F(2,2) + 4*d[20]*F(1,0)*F(1,2)*F(2,0)*F(2,2) +
    4*d[19]*F(1,1)*F(1,2)*F(2,0)*F(2,2) + 2*d[17]*pow(F(1,2),2)*F(2,0)*F(2,2) +
    2*d[10]*pow(F(1,0),2)*F(2,1)*F(2,2) + 4*d[13]*F(1,0)*F(1,1)*F(2,1)*F(2,2) +
    2*d[11]*pow(F(1,1),2)*F(2,1)*F(2,2) + 4*d[19]*F(1,0)*F(1,2)*F(2,1)*F(2,2) +
    4*d[14]*F(1,1)*F(1,2)*F(2,1)*F(2,2) + 2*d[12]*pow(F(1,2),2)*F(2,1)*F(2,2) +
    d[3]*pow(F(1,0),2)*pow(F(2,2),2) + 2*d[8]*F(1,0)*F(1,1)*pow(F(2,2),2) +
    d[4]*pow(F(1,1),2)*pow(F(2,2),2) + 2*d[17]*F(1,0)*F(1,2)*pow(F(2,2),2) +
    2*d[12]*F(1,1)*F(1,2)*pow(F(2,2),2) + d[5]*pow(F(1,2),2)*pow(F(2,2),2);
    
    c.d[5] = d[0]*pow(F(2,0),4) + 4*d[6]*pow(F(2,0),3)*F(2,1) +
    2*d[1]*pow(F(2,0),2)*pow(F(2,1),2) + 4*d[9]*pow(F(2,0),2)*pow(F(2,1),2) +
    4*d[7]*F(2,0)*pow(F(2,1),3) + d[2]*pow(F(2,1),4) +
    4*d[15]*pow(F(2,0),3)*F(2,2) + 4*d[10]*pow(F(2,0),2)*F(2,1)*F(2,2) +
    8*d[18]*pow(F(2,0),2)*F(2,1)*F(2,2) + 8*d[13]*F(2,0)*pow(F(2,1),2)*F(2,2) +
    4*d[16]*F(2,0)*pow(F(2,1),2)*F(2,2) + 4*d[11]*pow(F(2,1),3)*F(2,2) +
    2*d[3]*pow(F(2,0),2)*pow(F(2,2),2) + 4*d[20]*pow(F(2,0),2)*pow(F(2,2),2) +
    4*d[8]*F(2,0)*F(2,1)*pow(F(2,2),2) + 8*d[19]*F(2,0)*F(2,1)*pow(F(2,2),2) +
    2*d[4]*pow(F(2,1),2)*pow(F(2,2),2) + 4*d[14]*pow(F(2,1),2)*pow(F(2,2),2) +
    4*d[17]*F(2,0)*pow(F(2,2),3) + 4*d[12]*F(2,1)*pow(F(2,2),3) +
    d[5]*pow(F(2,2),4);
    
    c.d[6] = d[0]*pow(F(0,0),3)*F(1,0) + d[1]*F(0,0)*pow(F(0,1),2)*F(1,0) +
    2*d[9]*F(0,0)*pow(F(0,1),2)*F(1,0) + d[7]*pow(F(0,1),3)*F(1,0) +
    3*d[15]*pow(F(0,0),2)*F(0,2)*F(1,0) + 2*d[10]*F(0,0)*F(0,1)*F(0,2)*F(1,0) +
    4*d[18]*F(0,0)*F(0,1)*F(0,2)*F(1,0) + 2*d[13]*pow(F(0,1),2)*F(0,2)*F(1,0) +
    d[16]*pow(F(0,1),2)*F(0,2)*F(1,0) + d[3]*F(0,0)*pow(F(0,2),2)*F(1,0) +
    2*d[20]*F(0,0)*pow(F(0,2),2)*F(1,0) + d[8]*F(0,1)*pow(F(0,2),2)*F(1,0) +
    2*d[19]*F(0,1)*pow(F(0,2),2)*F(1,0) + d[17]*pow(F(0,2),3)*F(1,0) +
    d[1]*pow(F(0,0),2)*F(0,1)*F(1,1) + 2*d[9]*pow(F(0,0),2)*F(0,1)*F(1,1) +
    3*d[7]*F(0,0)*pow(F(0,1),2)*F(1,1) + d[2]*pow(F(0,1),3)*F(1,1) +
    d[10]*pow(F(0,0),2)*F(0,2)*F(1,1) + 2*d[18]*pow(F(0,0),2)*F(0,2)*F(1,1) +
    4*d[13]*F(0,0)*F(0,1)*F(0,2)*F(1,1) + 2*d[16]*F(0,0)*F(0,1)*F(0,2)*F(1,1) +
    3*d[11]*pow(F(0,1),2)*F(0,2)*F(1,1) + d[8]*F(0,0)*pow(F(0,2),2)*F(1,1) +
    2*d[19]*F(0,0)*pow(F(0,2),2)*F(1,1) + d[4]*F(0,1)*pow(F(0,2),2)*F(1,1) +
    2*d[14]*F(0,1)*pow(F(0,2),2)*F(1,1) + d[12]*pow(F(0,2),3)*F(1,1) +
    d[6]*pow(F(0,0),2)*(3*F(0,1)*F(1,0) + F(0,0)*F(1,1)) +
    d[15]*pow(F(0,0),3)*F(1,2) + d[10]*pow(F(0,0),2)*F(0,1)*F(1,2) +
    2*d[18]*pow(F(0,0),2)*F(0,1)*F(1,2) + 2*d[13]*F(0,0)*pow(F(0,1),2)*F(1,2) +
    d[16]*F(0,0)*pow(F(0,1),2)*F(1,2) + d[11]*pow(F(0,1),3)*F(1,2) +
    d[3]*pow(F(0,0),2)*F(0,2)*F(1,2) + 2*d[20]*pow(F(0,0),2)*F(0,2)*F(1,2) +
    2*d[8]*F(0,0)*F(0,1)*F(0,2)*F(1,2) + 4*d[19]*F(0,0)*F(0,1)*F(0,2)*F(1,2) +
    d[4]*pow(F(0,1),2)*F(0,2)*F(1,2) + 2*d[14]*pow(F(0,1),2)*F(0,2)*F(1,2) +
    3*d[17]*F(0,0)*pow(F(0,2),2)*F(1,2) + 3*d[12]*F(0,1)*pow(F(0,2),2)*F(1,2) +
    d[5]*pow(F(0,2),3)*F(1,2);
    
    c.d[7] = d[0]*F(0,0)*pow(F(1,0),3) + d[15]*F(0,2)*pow(F(1,0),3) +
    d[1]*F(0,1)*pow(F(1,0),2)*F(1,1) + 2*d[9]*F(0,1)*pow(F(1,0),2)*F(1,1) +
    d[10]*F(0,2)*pow(F(1,0),2)*F(1,1) + 2*d[18]*F(0,2)*pow(F(1,0),2)*F(1,1) +
    d[1]*F(0,0)*F(1,0)*pow(F(1,1),2) + 2*d[9]*F(0,0)*F(1,0)*pow(F(1,1),2) +
    3*d[7]*F(0,1)*F(1,0)*pow(F(1,1),2) + 2*d[13]*F(0,2)*F(1,0)*pow(F(1,1),2) +
    d[16]*F(0,2)*F(1,0)*pow(F(1,1),2) + d[7]*F(0,0)*pow(F(1,1),3) +
    d[2]*F(0,1)*pow(F(1,1),3) + d[11]*F(0,2)*pow(F(1,1),3) +
    d[6]*pow(F(1,0),2)*(F(0,1)*F(1,0) + 3*F(0,0)*F(1,1)) +
    3*d[15]*F(0,0)*pow(F(1,0),2)*F(1,2) + d[10]*F(0,1)*pow(F(1,0),2)*F(1,2) +
    2*d[18]*F(0,1)*pow(F(1,0),2)*F(1,2) + d[3]*F(0,2)*pow(F(1,0),2)*F(1,2) +
    2*d[20]*F(0,2)*pow(F(1,0),2)*F(1,2) + 2*d[10]*F(0,0)*F(1,0)*F(1,1)*F(1,2) +
    4*d[18]*F(0,0)*F(1,0)*F(1,1)*F(1,2) + 4*d[13]*F(0,1)*F(1,0)*F(1,1)*F(1,2) +
    2*d[16]*F(0,1)*F(1,0)*F(1,1)*F(1,2) + 2*d[8]*F(0,2)*F(1,0)*F(1,1)*F(1,2) +
    4*d[19]*F(0,2)*F(1,0)*F(1,1)*F(1,2) + 2*d[13]*F(0,0)*pow(F(1,1),2)*F(1,2) +
    d[16]*F(0,0)*pow(F(1,1),2)*F(1,2) + 3*d[11]*F(0,1)*pow(F(1,1),2)*F(1,2) +
    d[4]*F(0,2)*pow(F(1,1),2)*F(1,2) + 2*d[14]*F(0,2)*pow(F(1,1),2)*F(1,2) +
    d[3]*F(0,0)*F(1,0)*pow(F(1,2),2) + 2*d[20]*F(0,0)*F(1,0)*pow(F(1,2),2) +
    d[8]*F(0,1)*F(1,0)*pow(F(1,2),2) + 2*d[19]*F(0,1)*F(1,0)*pow(F(1,2),2) +
    3*d[17]*F(0,2)*F(1,0)*pow(F(1,2),2) + d[8]*F(0,0)*F(1,1)*pow(F(1,2),2) +
    2*d[19]*F(0,0)*F(1,1)*pow(F(1,2),2) + d[4]*F(0,1)*F(1,1)*pow(F(1,2),2) +
    2*d[14]*F(0,1)*F(1,1)*pow(F(1,2),2) + 3*d[12]*F(0,2)*F(1,1)*pow(F(1,2),2) +
    d[17]*F(0,0)*pow(F(1,2),3) + d[12]*F(0,1)*pow(F(1,2),3) +
    d[5]*F(0,2)*pow(F(1,2),3);
    
    c.d[8] = d[0]*F(0,0)*F(1,0)*pow(F(2,0),2) + d[15]*F(0,2)*F(1,0)*pow(F(2,0),2) +
    d[1]*F(0,1)*F(1,1)*pow(F(2,0),2) + d[10]*F(0,2)*F(1,1)*pow(F(2,0),2) +
    d[15]*F(0,0)*F(1,2)*pow(F(2,0),2) + d[10]*F(0,1)*F(1,2)*pow(F(2,0),2) +
    d[3]*F(0,2)*F(1,2)*pow(F(2,0),2) + 2*d[9]*F(0,1)*F(1,0)*F(2,0)*F(2,1) +
    2*d[18]*F(0,2)*F(1,0)*F(2,0)*F(2,1) + 2*d[9]*F(0,0)*F(1,1)*F(2,0)*F(2,1) +
    2*d[7]*F(0,1)*F(1,1)*F(2,0)*F(2,1) + 2*d[13]*F(0,2)*F(1,1)*F(2,0)*F(2,1) +
    2*d[18]*F(0,0)*F(1,2)*F(2,0)*F(2,1) + 2*d[13]*F(0,1)*F(1,2)*F(2,0)*F(2,1) +
    2*d[8]*F(0,2)*F(1,2)*F(2,0)*F(2,1) + d[1]*F(0,0)*F(1,0)*pow(F(2,1),2) +
    d[7]*F(0,1)*F(1,0)*pow(F(2,1),2) + d[16]*F(0,2)*F(1,0)*pow(F(2,1),2) +
    d[7]*F(0,0)*F(1,1)*pow(F(2,1),2) + d[2]*F(0,1)*F(1,1)*pow(F(2,1),2) +
    d[11]*F(0,2)*F(1,1)*pow(F(2,1),2) + d[16]*F(0,0)*F(1,2)*pow(F(2,1),2) +
    d[11]*F(0,1)*F(1,2)*pow(F(2,1),2) + d[4]*F(0,2)*F(1,2)*pow(F(2,1),2) +
    d[6]*F(2,0)*(F(0,1)*F(1,0)*F(2,0) + F(0,0)*(F(1,1)*F(2,0) + 2*F(1,0)*F(2,1))) +
    2*d[15]*F(0,0)*F(1,0)*F(2,0)*F(2,2) + 2*d[18]*F(0,1)*F(1,0)*F(2,0)*F(2,2) +
    2*d[20]*F(0,2)*F(1,0)*F(2,0)*F(2,2) + 2*d[18]*F(0,0)*F(1,1)*F(2,0)*F(2,2) +
    2*d[16]*F(0,1)*F(1,1)*F(2,0)*F(2,2) + 2*d[19]*F(0,2)*F(1,1)*F(2,0)*F(2,2) +
    2*d[20]*F(0,0)*F(1,2)*F(2,0)*F(2,2) + 2*d[19]*F(0,1)*F(1,2)*F(2,0)*F(2,2) +
    2*d[17]*F(0,2)*F(1,2)*F(2,0)*F(2,2) + 2*d[10]*F(0,0)*F(1,0)*F(2,1)*F(2,2) +
    2*d[13]*F(0,1)*F(1,0)*F(2,1)*F(2,2) + 2*d[19]*F(0,2)*F(1,0)*F(2,1)*F(2,2) +
    2*d[13]*F(0,0)*F(1,1)*F(2,1)*F(2,2) + 2*d[11]*F(0,1)*F(1,1)*F(2,1)*F(2,2) +
    2*d[14]*F(0,2)*F(1,1)*F(2,1)*F(2,2) + 2*d[19]*F(0,0)*F(1,2)*F(2,1)*F(2,2) +
    2*d[14]*F(0,1)*F(1,2)*F(2,1)*F(2,2) + 2*d[12]*F(0,2)*F(1,2)*F(2,1)*F(2,2) +
    d[3]*F(0,0)*F(1,0)*pow(F(2,2),2) + d[8]*F(0,1)*F(1,0)*pow(F(2,2),2) +
    d[17]*F(0,2)*F(1,0)*pow(F(2,2),2) + d[8]*F(0,0)*F(1,1)*pow(F(2,2),2) +
    d[4]*F(0,1)*F(1,1)*pow(F(2,2),2) + d[12]*F(0,2)*F(1,1)*pow(F(2,2),2) +
    d[17]*F(0,0)*F(1,2)*pow(F(2,2),2) + d[12]*F(0,1)*F(1,2)*pow(F(2,2),2) +
    d[5]*F(0,2)*F(1,2)*pow(F(2,2),2);
    
    c.d[9] = d[0]*pow(F(0,0),2)*pow(F(1,0),2) + d[9]*pow(F(0,1),2)*pow(F(1,0),2) +
    2*d[15]*F(0,0)*F(0,2)*pow(F(1,0),2) + 2*d[18]*F(0,1)*F(0,2)*pow(F(1,0),2) +
    d[20]*pow(F(0,2),2)*pow(F(1,0),2) + 2*d[1]*F(0,0)*F(0,1)*F(1,0)*F(1,1) +
    2*d[9]*F(0,0)*F(0,1)*F(1,0)*F(1,1) + 2*d[7]*pow(F(0,1),2)*F(1,0)*F(1,1) +
    2*d[10]*F(0,0)*F(0,2)*F(1,0)*F(1,1) + 2*d[18]*F(0,0)*F(0,2)*F(1,0)*F(1,1) +
    2*d[13]*F(0,1)*F(0,2)*F(1,0)*F(1,1) + 2*d[16]*F(0,1)*F(0,2)*F(1,0)*F(1,1) +
    2*d[19]*pow(F(0,2),2)*F(1,0)*F(1,1) + d[9]*pow(F(0,0),2)*pow(F(1,1),2) +
    2*d[7]*F(0,0)*F(0,1)*pow(F(1,1),2) + d[2]*pow(F(0,1),2)*pow(F(1,1),2) +
    2*d[13]*F(0,0)*F(0,2)*pow(F(1,1),2) + 2*d[11]*F(0,1)*F(0,2)*pow(F(1,1),2) +
    d[14]*pow(F(0,2),2)*pow(F(1,1),2) +
    2*d[6]*F(0,0)*F(1,0)*(F(0,1)*F(1,0) + F(0,0)*F(1,1)) +
    2*d[15]*pow(F(0,0),2)*F(1,0)*F(1,2) + 2*d[10]*F(0,0)*F(0,1)*F(1,0)*F(1,2) +
    2*d[18]*F(0,0)*F(0,1)*F(1,0)*F(1,2) + 2*d[13]*pow(F(0,1),2)*F(1,0)*F(1,2) +
    2*d[3]*F(0,0)*F(0,2)*F(1,0)*F(1,2) + 2*d[20]*F(0,0)*F(0,2)*F(1,0)*F(1,2) +
    2*d[8]*F(0,1)*F(0,2)*F(1,0)*F(1,2) + 2*d[19]*F(0,1)*F(0,2)*F(1,0)*F(1,2) +
    2*d[17]*pow(F(0,2),2)*F(1,0)*F(1,2) + 2*d[18]*pow(F(0,0),2)*F(1,1)*F(1,2) +
    2*d[13]*F(0,0)*F(0,1)*F(1,1)*F(1,2) + 2*d[16]*F(0,0)*F(0,1)*F(1,1)*F(1,2) +
    2*d[11]*pow(F(0,1),2)*F(1,1)*F(1,2) + 2*d[8]*F(0,0)*F(0,2)*F(1,1)*F(1,2) +
    2*d[19]*F(0,0)*F(0,2)*F(1,1)*F(1,2) + 2*d[4]*F(0,1)*F(0,2)*F(1,1)*F(1,2) +
    2*d[14]*F(0,1)*F(0,2)*F(1,1)*F(1,2) + 2*d[12]*pow(F(0,2),2)*F(1,1)*F(1,2) +
    d[20]*pow(F(0,0),2)*pow(F(1,2),2) + 2*d[19]*F(0,0)*F(0,1)*pow(F(1,2),2) +
    d[14]*pow(F(0,1),2)*pow(F(1,2),2) + 2*d[17]*F(0,0)*F(0,2)*pow(F(1,2),2) +
    2*d[12]*F(0,1)*F(0,2)*pow(F(1,2),2) + d[5]*pow(F(0,2),2)*pow(F(1,2),2);
    
    c.d[10] = d[0]*pow(F(0,0),2)*F(1,0)*F(2,0) + d[1]*pow(F(0,1),2)*F(1,0)*F(2,0) +
    2*d[15]*F(0,0)*F(0,2)*F(1,0)*F(2,0) + 2*d[10]*F(0,1)*F(0,2)*F(1,0)*F(2,0) +
    d[3]*pow(F(0,2),2)*F(1,0)*F(2,0) + 2*d[9]*F(0,0)*F(0,1)*F(1,1)*F(2,0) +
    d[7]*pow(F(0,1),2)*F(1,1)*F(2,0) + 2*d[18]*F(0,0)*F(0,2)*F(1,1)*F(2,0) +
    2*d[13]*F(0,1)*F(0,2)*F(1,1)*F(2,0) + d[8]*pow(F(0,2),2)*F(1,1)*F(2,0) +
    d[15]*pow(F(0,0),2)*F(1,2)*F(2,0) + 2*d[18]*F(0,0)*F(0,1)*F(1,2)*F(2,0) +
    d[16]*pow(F(0,1),2)*F(1,2)*F(2,0) + 2*d[20]*F(0,0)*F(0,2)*F(1,2)*F(2,0) +
    2*d[19]*F(0,1)*F(0,2)*F(1,2)*F(2,0) + d[17]*pow(F(0,2),2)*F(1,2)*F(2,0) +
    2*d[9]*F(0,0)*F(0,1)*F(1,0)*F(2,1) + d[7]*pow(F(0,1),2)*F(1,0)*F(2,1) +
    2*d[18]*F(0,0)*F(0,2)*F(1,0)*F(2,1) + 2*d[13]*F(0,1)*F(0,2)*F(1,0)*F(2,1) +
    d[8]*pow(F(0,2),2)*F(1,0)*F(2,1) + d[1]*pow(F(0,0),2)*F(1,1)*F(2,1) +
    2*d[7]*F(0,0)*F(0,1)*F(1,1)*F(2,1) + d[2]*pow(F(0,1),2)*F(1,1)*F(2,1) +
    2*d[16]*F(0,0)*F(0,2)*F(1,1)*F(2,1) + 2*d[11]*F(0,1)*F(0,2)*F(1,1)*F(2,1) +
    d[4]*pow(F(0,2),2)*F(1,1)*F(2,1) + d[10]*pow(F(0,0),2)*F(1,2)*F(2,1) +
    2*d[13]*F(0,0)*F(0,1)*F(1,2)*F(2,1) + d[11]*pow(F(0,1),2)*F(1,2)*F(2,1) +
    2*d[19]*F(0,0)*F(0,2)*F(1,2)*F(2,1) + 2*d[14]*F(0,1)*F(0,2)*F(1,2)*F(2,1) +
    d[12]*pow(F(0,2),2)*F(1,2)*F(2,1) +
    d[6]*F(0,0)*(2*F(0,1)*F(1,0)*F(2,0) + F(0,0)*(F(1,1)*F(2,0) + F(1,0)*F(2,1))) +
    d[15]*pow(F(0,0),2)*F(1,0)*F(2,2) + 2*d[18]*F(0,0)*F(0,1)*F(1,0)*F(2,2) +
    d[16]*pow(F(0,1),2)*F(1,0)*F(2,2) + 2*d[20]*F(0,0)*F(0,2)*F(1,0)*F(2,2) +
    2*d[19]*F(0,1)*F(0,2)*F(1,0)*F(2,2) + d[17]*pow(F(0,2),2)*F(1,0)*F(2,2) +
    d[10]*pow(F(0,0),2)*F(1,1)*F(2,2) + 2*d[13]*F(0,0)*F(0,1)*F(1,1)*F(2,2) +
    d[11]*pow(F(0,1),2)*F(1,1)*F(2,2) + 2*d[19]*F(0,0)*F(0,2)*F(1,1)*F(2,2) +
    2*d[14]*F(0,1)*F(0,2)*F(1,1)*F(2,2) + d[12]*pow(F(0,2),2)*F(1,1)*F(2,2) +
    d[3]*pow(F(0,0),2)*F(1,2)*F(2,2) + 2*d[8]*F(0,0)*F(0,1)*F(1,2)*F(2,2) +
    d[4]*pow(F(0,1),2)*F(1,2)*F(2,2) + 2*d[17]*F(0,0)*F(0,2)*F(1,2)*F(2,2) +
    2*d[12]*F(0,1)*F(0,2)*F(1,2)*F(2,2) + d[5]*pow(F(0,2),2)*F(1,2)*F(2,2);
    
    c.d[11] = d[0]*pow(F(1,0),3)*F(2,0) + d[1]*F(1,0)*pow(F(1,1),2)*F(2,0) +
    2*d[9]*F(1,0)*pow(F(1,1),2)*F(2,0) + d[7]*pow(F(1,1),3)*F(2,0) +
    3*d[15]*pow(F(1,0),2)*F(1,2)*F(2,0) + 2*d[10]*F(1,0)*F(1,1)*F(1,2)*F(2,0) +
    4*d[18]*F(1,0)*F(1,1)*F(1,2)*F(2,0) + 2*d[13]*pow(F(1,1),2)*F(1,2)*F(2,0) +
    d[16]*pow(F(1,1),2)*F(1,2)*F(2,0) + d[3]*F(1,0)*pow(F(1,2),2)*F(2,0) +
    2*d[20]*F(1,0)*pow(F(1,2),2)*F(2,0) + d[8]*F(1,1)*pow(F(1,2),2)*F(2,0) +
    2*d[19]*F(1,1)*pow(F(1,2),2)*F(2,0) + d[17]*pow(F(1,2),3)*F(2,0) +
    d[1]*pow(F(1,0),2)*F(1,1)*F(2,1) + 2*d[9]*pow(F(1,0),2)*F(1,1)*F(2,1) +
    3*d[7]*F(1,0)*pow(F(1,1),2)*F(2,1) + d[2]*pow(F(1,1),3)*F(2,1) +
    d[10]*pow(F(1,0),2)*F(1,2)*F(2,1) + 2*d[18]*pow(F(1,0),2)*F(1,2)*F(2,1) +
    4*d[13]*F(1,0)*F(1,1)*F(1,2)*F(2,1) + 2*d[16]*F(1,0)*F(1,1)*F(1,2)*F(2,1) +
    3*d[11]*pow(F(1,1),2)*F(1,2)*F(2,1) + d[8]*F(1,0)*pow(F(1,2),2)*F(2,1) +
    2*d[19]*F(1,0)*pow(F(1,2),2)*F(2,1) + d[4]*F(1,1)*pow(F(1,2),2)*F(2,1) +
    2*d[14]*F(1,1)*pow(F(1,2),2)*F(2,1) + d[12]*pow(F(1,2),3)*F(2,1) +
    d[6]*pow(F(1,0),2)*(3*F(1,1)*F(2,0) + F(1,0)*F(2,1)) +
    d[15]*pow(F(1,0),3)*F(2,2) + d[10]*pow(F(1,0),2)*F(1,1)*F(2,2) +
    2*d[18]*pow(F(1,0),2)*F(1,1)*F(2,2) + 2*d[13]*F(1,0)*pow(F(1,1),2)*F(2,2) +
    d[16]*F(1,0)*pow(F(1,1),2)*F(2,2) + d[11]*pow(F(1,1),3)*F(2,2) +
    d[3]*pow(F(1,0),2)*F(1,2)*F(2,2) + 2*d[20]*pow(F(1,0),2)*F(1,2)*F(2,2) +
    2*d[8]*F(1,0)*F(1,1)*F(1,2)*F(2,2) + 4*d[19]*F(1,0)*F(1,1)*F(1,2)*F(2,2) +
    d[4]*pow(F(1,1),2)*F(1,2)*F(2,2) + 2*d[14]*pow(F(1,1),2)*F(1,2)*F(2,2) +
    3*d[17]*F(1,0)*pow(F(1,2),2)*F(2,2) + 3*d[12]*F(1,1)*pow(F(1,2),2)*F(2,2) +
    d[5]*pow(F(1,2),3)*F(2,2);
    
    c.d[12] = d[0]*F(1,0)*pow(F(2,0),3) + d[15]*F(1,2)*pow(F(2,0),3) +
    d[1]*F(1,1)*pow(F(2,0),2)*F(2,1) + 2*d[9]*F(1,1)*pow(F(2,0),2)*F(2,1) +
    d[10]*F(1,2)*pow(F(2,0),2)*F(2,1) + 2*d[18]*F(1,2)*pow(F(2,0),2)*F(2,1) +
    d[1]*F(1,0)*F(2,0)*pow(F(2,1),2) + 2*d[9]*F(1,0)*F(2,0)*pow(F(2,1),2) +
    3*d[7]*F(1,1)*F(2,0)*pow(F(2,1),2) + 2*d[13]*F(1,2)*F(2,0)*pow(F(2,1),2) +
    d[16]*F(1,2)*F(2,0)*pow(F(2,1),2) + d[7]*F(1,0)*pow(F(2,1),3) +
    d[2]*F(1,1)*pow(F(2,1),3) + d[11]*F(1,2)*pow(F(2,1),3) +
    d[6]*pow(F(2,0),2)*(F(1,1)*F(2,0) + 3*F(1,0)*F(2,1)) +
    3*d[15]*F(1,0)*pow(F(2,0),2)*F(2,2) + d[10]*F(1,1)*pow(F(2,0),2)*F(2,2) +
    2*d[18]*F(1,1)*pow(F(2,0),2)*F(2,2) + d[3]*F(1,2)*pow(F(2,0),2)*F(2,2) +
    2*d[20]*F(1,2)*pow(F(2,0),2)*F(2,2) + 2*d[10]*F(1,0)*F(2,0)*F(2,1)*F(2,2) +
    4*d[18]*F(1,0)*F(2,0)*F(2,1)*F(2,2) + 4*d[13]*F(1,1)*F(2,0)*F(2,1)*F(2,2) +
    2*d[16]*F(1,1)*F(2,0)*F(2,1)*F(2,2) + 2*d[8]*F(1,2)*F(2,0)*F(2,1)*F(2,2) +
    4*d[19]*F(1,2)*F(2,0)*F(2,1)*F(2,2) + 2*d[13]*F(1,0)*pow(F(2,1),2)*F(2,2) +
    d[16]*F(1,0)*pow(F(2,1),2)*F(2,2) + 3*d[11]*F(1,1)*pow(F(2,1),2)*F(2,2) +
    d[4]*F(1,2)*pow(F(2,1),2)*F(2,2) + 2*d[14]*F(1,2)*pow(F(2,1),2)*F(2,2) +
    d[3]*F(1,0)*F(2,0)*pow(F(2,2),2) + 2*d[20]*F(1,0)*F(2,0)*pow(F(2,2),2) +
    d[8]*F(1,1)*F(2,0)*pow(F(2,2),2) + 2*d[19]*F(1,1)*F(2,0)*pow(F(2,2),2) +
    3*d[17]*F(1,2)*F(2,0)*pow(F(2,2),2) + d[8]*F(1,0)*F(2,1)*pow(F(2,2),2) +
    2*d[19]*F(1,0)*F(2,1)*pow(F(2,2),2) + d[4]*F(1,1)*F(2,1)*pow(F(2,2),2) +
    2*d[14]*F(1,1)*F(2,1)*pow(F(2,2),2) + 3*d[12]*F(1,2)*F(2,1)*pow(F(2,2),2) +
    d[17]*F(1,0)*pow(F(2,2),3) + d[12]*F(1,1)*pow(F(2,2),3) +
    d[5]*F(1,2)*pow(F(2,2),3);
    
    c.d[13] = d[0]*F(0,0)*pow(F(1,0),2)*F(2,0) + d[15]*F(0,2)*pow(F(1,0),2)*F(2,0) +
    d[1]*F(0,1)*F(1,0)*F(1,1)*F(2,0) + d[9]*F(0,1)*F(1,0)*F(1,1)*F(2,0) +
    d[10]*F(0,2)*F(1,0)*F(1,1)*F(2,0) + d[18]*F(0,2)*F(1,0)*F(1,1)*F(2,0) +
    d[9]*F(0,0)*pow(F(1,1),2)*F(2,0) + d[7]*F(0,1)*pow(F(1,1),2)*F(2,0) +
    d[13]*F(0,2)*pow(F(1,1),2)*F(2,0) + 2*d[15]*F(0,0)*F(1,0)*F(1,2)*F(2,0) +
    d[10]*F(0,1)*F(1,0)*F(1,2)*F(2,0) + d[18]*F(0,1)*F(1,0)*F(1,2)*F(2,0) +
    d[3]*F(0,2)*F(1,0)*F(1,2)*F(2,0) + d[20]*F(0,2)*F(1,0)*F(1,2)*F(2,0) +
    2*d[18]*F(0,0)*F(1,1)*F(1,2)*F(2,0) + d[13]*F(0,1)*F(1,1)*F(1,2)*F(2,0) +
    d[16]*F(0,1)*F(1,1)*F(1,2)*F(2,0) + d[8]*F(0,2)*F(1,1)*F(1,2)*F(2,0) +
    d[19]*F(0,2)*F(1,1)*F(1,2)*F(2,0) + d[20]*F(0,0)*pow(F(1,2),2)*F(2,0) +
    d[19]*F(0,1)*pow(F(1,2),2)*F(2,0) + d[17]*F(0,2)*pow(F(1,2),2)*F(2,0) +
    d[9]*F(0,1)*pow(F(1,0),2)*F(2,1) + d[18]*F(0,2)*pow(F(1,0),2)*F(2,1) +
    d[1]*F(0,0)*F(1,0)*F(1,1)*F(2,1) + d[9]*F(0,0)*F(1,0)*F(1,1)*F(2,1) +
    2*d[7]*F(0,1)*F(1,0)*F(1,1)*F(2,1) + d[13]*F(0,2)*F(1,0)*F(1,1)*F(2,1) +
    d[16]*F(0,2)*F(1,0)*F(1,1)*F(2,1) + d[7]*F(0,0)*pow(F(1,1),2)*F(2,1) +
    d[2]*F(0,1)*pow(F(1,1),2)*F(2,1) + d[11]*F(0,2)*pow(F(1,1),2)*F(2,1) +
    d[10]*F(0,0)*F(1,0)*F(1,2)*F(2,1) + d[18]*F(0,0)*F(1,0)*F(1,2)*F(2,1) +
    2*d[13]*F(0,1)*F(1,0)*F(1,2)*F(2,1) + d[8]*F(0,2)*F(1,0)*F(1,2)*F(2,1) +
    d[19]*F(0,2)*F(1,0)*F(1,2)*F(2,1) + d[13]*F(0,0)*F(1,1)*F(1,2)*F(2,1) +
    d[16]*F(0,0)*F(1,1)*F(1,2)*F(2,1) + 2*d[11]*F(0,1)*F(1,1)*F(1,2)*F(2,1) +
    d[4]*F(0,2)*F(1,1)*F(1,2)*F(2,1) + d[14]*F(0,2)*F(1,1)*F(1,2)*F(2,1) +
    d[19]*F(0,0)*pow(F(1,2),2)*F(2,1) + d[14]*F(0,1)*pow(F(1,2),2)*F(2,1) +
    d[12]*F(0,2)*pow(F(1,2),2)*F(2,1) +
    d[6]*F(1,0)*(F(0,1)*F(1,0)*F(2,0) + F(0,0)*(2*F(1,1)*F(2,0) + F(1,0)*F(2,1))) +
    d[15]*F(0,0)*pow(F(1,0),2)*F(2,2) + d[18]*F(0,1)*pow(F(1,0),2)*F(2,2) +
    d[20]*F(0,2)*pow(F(1,0),2)*F(2,2) + d[10]*F(0,0)*F(1,0)*F(1,1)*F(2,2) +
    d[18]*F(0,0)*F(1,0)*F(1,1)*F(2,2) + d[13]*F(0,1)*F(1,0)*F(1,1)*F(2,2) +
    d[16]*F(0,1)*F(1,0)*F(1,1)*F(2,2) + 2*d[19]*F(0,2)*F(1,0)*F(1,1)*F(2,2) +
    d[13]*F(0,0)*pow(F(1,1),2)*F(2,2) + d[11]*F(0,1)*pow(F(1,1),2)*F(2,2) +
    d[14]*F(0,2)*pow(F(1,1),2)*F(2,2) + d[3]*F(0,0)*F(1,0)*F(1,2)*F(2,2) +
    d[20]*F(0,0)*F(1,0)*F(1,2)*F(2,2) + d[8]*F(0,1)*F(1,0)*F(1,2)*F(2,2) +
    d[19]*F(0,1)*F(1,0)*F(1,2)*F(2,2) + 2*d[17]*F(0,2)*F(1,0)*F(1,2)*F(2,2) +
    d[8]*F(0,0)*F(1,1)*F(1,2)*F(2,2) + d[19]*F(0,0)*F(1,1)*F(1,2)*F(2,2) +
    d[4]*F(0,1)*F(1,1)*F(1,2)*F(2,2) + d[14]*F(0,1)*F(1,1)*F(1,2)*F(2,2) +
    2*d[12]*F(0,2)*F(1,1)*F(1,2)*F(2,2) + d[17]*F(0,0)*pow(F(1,2),2)*F(2,2) +
    d[12]*F(0,1)*pow(F(1,2),2)*F(2,2) + d[5]*F(0,2)*pow(F(1,2),2)*F(2,2);
    
    c.d[14] = d[0]*pow(F(1,0),2)*pow(F(2,0),2) + d[9]*pow(F(1,1),2)*pow(F(2,0),2) +
    2*d[15]*F(1,0)*F(1,2)*pow(F(2,0),2) + 2*d[18]*F(1,1)*F(1,2)*pow(F(2,0),2) +
    d[20]*pow(F(1,2),2)*pow(F(2,0),2) + 2*d[1]*F(1,0)*F(1,1)*F(2,0)*F(2,1) +
    2*d[9]*F(1,0)*F(1,1)*F(2,0)*F(2,1) + 2*d[7]*pow(F(1,1),2)*F(2,0)*F(2,1) +
    2*d[10]*F(1,0)*F(1,2)*F(2,0)*F(2,1) + 2*d[18]*F(1,0)*F(1,2)*F(2,0)*F(2,1) +
    2*d[13]*F(1,1)*F(1,2)*F(2,0)*F(2,1) + 2*d[16]*F(1,1)*F(1,2)*F(2,0)*F(2,1) +
    2*d[19]*pow(F(1,2),2)*F(2,0)*F(2,1) + d[9]*pow(F(1,0),2)*pow(F(2,1),2) +
    2*d[7]*F(1,0)*F(1,1)*pow(F(2,1),2) + d[2]*pow(F(1,1),2)*pow(F(2,1),2) +
    2*d[13]*F(1,0)*F(1,2)*pow(F(2,1),2) + 2*d[11]*F(1,1)*F(1,2)*pow(F(2,1),2) +
    d[14]*pow(F(1,2),2)*pow(F(2,1),2) +
    2*d[6]*F(1,0)*F(2,0)*(F(1,1)*F(2,0) + F(1,0)*F(2,1)) +
    2*d[15]*pow(F(1,0),2)*F(2,0)*F(2,2) + 2*d[10]*F(1,0)*F(1,1)*F(2,0)*F(2,2) +
    2*d[18]*F(1,0)*F(1,1)*F(2,0)*F(2,2) + 2*d[13]*pow(F(1,1),2)*F(2,0)*F(2,2) +
    2*d[3]*F(1,0)*F(1,2)*F(2,0)*F(2,2) + 2*d[20]*F(1,0)*F(1,2)*F(2,0)*F(2,2) +
    2*d[8]*F(1,1)*F(1,2)*F(2,0)*F(2,2) + 2*d[19]*F(1,1)*F(1,2)*F(2,0)*F(2,2) +
    2*d[17]*pow(F(1,2),2)*F(2,0)*F(2,2) + 2*d[18]*pow(F(1,0),2)*F(2,1)*F(2,2) +
    2*d[13]*F(1,0)*F(1,1)*F(2,1)*F(2,2) + 2*d[16]*F(1,0)*F(1,1)*F(2,1)*F(2,2) +
    2*d[11]*pow(F(1,1),2)*F(2,1)*F(2,2) + 2*d[8]*F(1,0)*F(1,2)*F(2,1)*F(2,2) +
    2*d[19]*F(1,0)*F(1,2)*F(2,1)*F(2,2) + 2*d[4]*F(1,1)*F(1,2)*F(2,1)*F(2,2) +
    2*d[14]*F(1,1)*F(1,2)*F(2,1)*F(2,2) + 2*d[12]*pow(F(1,2),2)*F(2,1)*F(2,2) +
    d[20]*pow(F(1,0),2)*pow(F(2,2),2) + 2*d[19]*F(1,0)*F(1,1)*pow(F(2,2),2) +
    d[14]*pow(F(1,1),2)*pow(F(2,2),2) + 2*d[17]*F(1,0)*F(1,2)*pow(F(2,2),2) +
    2*d[12]*F(1,1)*F(1,2)*pow(F(2,2),2) + d[5]*pow(F(1,2),2)*pow(F(2,2),2);
    
    c.d[15] = d[0]*pow(F(0,0),3)*F(2,0) + d[1]*F(0,0)*pow(F(0,1),2)*F(2,0) +
    2*d[9]*F(0,0)*pow(F(0,1),2)*F(2,0) + d[7]*pow(F(0,1),3)*F(2,0) +
    3*d[15]*pow(F(0,0),2)*F(0,2)*F(2,0) + 2*d[10]*F(0,0)*F(0,1)*F(0,2)*F(2,0) +
    4*d[18]*F(0,0)*F(0,1)*F(0,2)*F(2,0) + 2*d[13]*pow(F(0,1),2)*F(0,2)*F(2,0) +
    d[16]*pow(F(0,1),2)*F(0,2)*F(2,0) + d[3]*F(0,0)*pow(F(0,2),2)*F(2,0) +
    2*d[20]*F(0,0)*pow(F(0,2),2)*F(2,0) + d[8]*F(0,1)*pow(F(0,2),2)*F(2,0) +
    2*d[19]*F(0,1)*pow(F(0,2),2)*F(2,0) + d[17]*pow(F(0,2),3)*F(2,0) +
    d[1]*pow(F(0,0),2)*F(0,1)*F(2,1) + 2*d[9]*pow(F(0,0),2)*F(0,1)*F(2,1) +
    3*d[7]*F(0,0)*pow(F(0,1),2)*F(2,1) + d[2]*pow(F(0,1),3)*F(2,1) +
    d[10]*pow(F(0,0),2)*F(0,2)*F(2,1) + 2*d[18]*pow(F(0,0),2)*F(0,2)*F(2,1) +
    4*d[13]*F(0,0)*F(0,1)*F(0,2)*F(2,1) + 2*d[16]*F(0,0)*F(0,1)*F(0,2)*F(2,1) +
    3*d[11]*pow(F(0,1),2)*F(0,2)*F(2,1) + d[8]*F(0,0)*pow(F(0,2),2)*F(2,1) +
    2*d[19]*F(0,0)*pow(F(0,2),2)*F(2,1) + d[4]*F(0,1)*pow(F(0,2),2)*F(2,1) +
    2*d[14]*F(0,1)*pow(F(0,2),2)*F(2,1) + d[12]*pow(F(0,2),3)*F(2,1) +
    d[6]*pow(F(0,0),2)*(3*F(0,1)*F(2,0) + F(0,0)*F(2,1)) +
    d[15]*pow(F(0,0),3)*F(2,2) + d[10]*pow(F(0,0),2)*F(0,1)*F(2,2) +
    2*d[18]*pow(F(0,0),2)*F(0,1)*F(2,2) + 2*d[13]*F(0,0)*pow(F(0,1),2)*F(2,2) +
    d[16]*F(0,0)*pow(F(0,1),2)*F(2,2) + d[11]*pow(F(0,1),3)*F(2,2) +
    d[3]*pow(F(0,0),2)*F(0,2)*F(2,2) + 2*d[20]*pow(F(0,0),2)*F(0,2)*F(2,2) +
    2*d[8]*F(0,0)*F(0,1)*F(0,2)*F(2,2) + 4*d[19]*F(0,0)*F(0,1)*F(0,2)*F(2,2) +
    d[4]*pow(F(0,1),2)*F(0,2)*F(2,2) + 2*d[14]*pow(F(0,1),2)*F(0,2)*F(2,2) +
    3*d[17]*F(0,0)*pow(F(0,2),2)*F(2,2) + 3*d[12]*F(0,1)*pow(F(0,2),2)*F(2,2) +
    d[5]*pow(F(0,2),3)*F(2,2);
    
    c.d[16] = d[0]*F(0,0)*pow(F(1,0),2)*F(2,0) + d[15]*F(0,2)*pow(F(1,0),2)*F(2,0) +
    2*d[9]*F(0,1)*F(1,0)*F(1,1)*F(2,0) + 2*d[18]*F(0,2)*F(1,0)*F(1,1)*F(2,0) +
    d[1]*F(0,0)*pow(F(1,1),2)*F(2,0) + d[7]*F(0,1)*pow(F(1,1),2)*F(2,0) +
    d[16]*F(0,2)*pow(F(1,1),2)*F(2,0) + 2*d[15]*F(0,0)*F(1,0)*F(1,2)*F(2,0) +
    2*d[18]*F(0,1)*F(1,0)*F(1,2)*F(2,0) + 2*d[20]*F(0,2)*F(1,0)*F(1,2)*F(2,0) +
    2*d[10]*F(0,0)*F(1,1)*F(1,2)*F(2,0) + 2*d[13]*F(0,1)*F(1,1)*F(1,2)*F(2,0) +
    2*d[19]*F(0,2)*F(1,1)*F(1,2)*F(2,0) + d[3]*F(0,0)*pow(F(1,2),2)*F(2,0) +
    d[8]*F(0,1)*pow(F(1,2),2)*F(2,0) + d[17]*F(0,2)*pow(F(1,2),2)*F(2,0) +
    d[1]*F(0,1)*pow(F(1,0),2)*F(2,1) + d[10]*F(0,2)*pow(F(1,0),2)*F(2,1) +
    2*d[9]*F(0,0)*F(1,0)*F(1,1)*F(2,1) + 2*d[7]*F(0,1)*F(1,0)*F(1,1)*F(2,1) +
    2*d[13]*F(0,2)*F(1,0)*F(1,1)*F(2,1) + d[7]*F(0,0)*pow(F(1,1),2)*F(2,1) +
    d[2]*F(0,1)*pow(F(1,1),2)*F(2,1) + d[11]*F(0,2)*pow(F(1,1),2)*F(2,1) +
    2*d[18]*F(0,0)*F(1,0)*F(1,2)*F(2,1) + 2*d[16]*F(0,1)*F(1,0)*F(1,2)*F(2,1) +
    2*d[19]*F(0,2)*F(1,0)*F(1,2)*F(2,1) + 2*d[13]*F(0,0)*F(1,1)*F(1,2)*F(2,1) +
    2*d[11]*F(0,1)*F(1,1)*F(1,2)*F(2,1) + 2*d[14]*F(0,2)*F(1,1)*F(1,2)*F(2,1) +
    d[8]*F(0,0)*pow(F(1,2),2)*F(2,1) + d[4]*F(0,1)*pow(F(1,2),2)*F(2,1) +
    d[12]*F(0,2)*pow(F(1,2),2)*F(2,1) +
    d[6]*F(1,0)*(F(0,1)*F(1,0)*F(2,0) + F(0,0)*(2*F(1,1)*F(2,0) + F(1,0)*F(2,1))) +
    d[15]*F(0,0)*pow(F(1,0),2)*F(2,2) + d[10]*F(0,1)*pow(F(1,0),2)*F(2,2) +
    d[3]*F(0,2)*pow(F(1,0),2)*F(2,2) + 2*d[18]*F(0,0)*F(1,0)*F(1,1)*F(2,2) +
    2*d[13]*F(0,1)*F(1,0)*F(1,1)*F(2,2) + 2*d[8]*F(0,2)*F(1,0)*F(1,1)*F(2,2) +
    d[16]*F(0,0)*pow(F(1,1),2)*F(2,2) + d[11]*F(0,1)*pow(F(1,1),2)*F(2,2) +
    d[4]*F(0,2)*pow(F(1,1),2)*F(2,2) + 2*d[20]*F(0,0)*F(1,0)*F(1,2)*F(2,2) +
    2*d[19]*F(0,1)*F(1,0)*F(1,2)*F(2,2) + 2*d[17]*F(0,2)*F(1,0)*F(1,2)*F(2,2) +
    2*d[19]*F(0,0)*F(1,1)*F(1,2)*F(2,2) + 2*d[14]*F(0,1)*F(1,1)*F(1,2)*F(2,2) +
    2*d[12]*F(0,2)*F(1,1)*F(1,2)*F(2,2) + d[17]*F(0,0)*pow(F(1,2),2)*F(2,2) +
    d[12]*F(0,1)*pow(F(1,2),2)*F(2,2) + d[5]*F(0,2)*pow(F(1,2),2)*F(2,2);
    
    c.d[17] = d[0]*F(0,0)*pow(F(2,0),3) + d[15]*F(0,2)*pow(F(2,0),3) +
    d[1]*F(0,1)*pow(F(2,0),2)*F(2,1) + 2*d[9]*F(0,1)*pow(F(2,0),2)*F(2,1) +
    d[10]*F(0,2)*pow(F(2,0),2)*F(2,1) + 2*d[18]*F(0,2)*pow(F(2,0),2)*F(2,1) +
    d[1]*F(0,0)*F(2,0)*pow(F(2,1),2) + 2*d[9]*F(0,0)*F(2,0)*pow(F(2,1),2) +
    3*d[7]*F(0,1)*F(2,0)*pow(F(2,1),2) + 2*d[13]*F(0,2)*F(2,0)*pow(F(2,1),2) +
    d[16]*F(0,2)*F(2,0)*pow(F(2,1),2) + d[7]*F(0,0)*pow(F(2,1),3) +
    d[2]*F(0,1)*pow(F(2,1),3) + d[11]*F(0,2)*pow(F(2,1),3) +
    d[6]*pow(F(2,0),2)*(F(0,1)*F(2,0) + 3*F(0,0)*F(2,1)) +
    3*d[15]*F(0,0)*pow(F(2,0),2)*F(2,2) + d[10]*F(0,1)*pow(F(2,0),2)*F(2,2) +
    2*d[18]*F(0,1)*pow(F(2,0),2)*F(2,2) + d[3]*F(0,2)*pow(F(2,0),2)*F(2,2) +
    2*d[20]*F(0,2)*pow(F(2,0),2)*F(2,2) + 2*d[10]*F(0,0)*F(2,0)*F(2,1)*F(2,2) +
    4*d[18]*F(0,0)*F(2,0)*F(2,1)*F(2,2) + 4*d[13]*F(0,1)*F(2,0)*F(2,1)*F(2,2) +
    2*d[16]*F(0,1)*F(2,0)*F(2,1)*F(2,2) + 2*d[8]*F(0,2)*F(2,0)*F(2,1)*F(2,2) +
    4*d[19]*F(0,2)*F(2,0)*F(2,1)*F(2,2) + 2*d[13]*F(0,0)*pow(F(2,1),2)*F(2,2) +
    d[16]*F(0,0)*pow(F(2,1),2)*F(2,2) + 3*d[11]*F(0,1)*pow(F(2,1),2)*F(2,2) +
    d[4]*F(0,2)*pow(F(2,1),2)*F(2,2) + 2*d[14]*F(0,2)*pow(F(2,1),2)*F(2,2) +
    d[3]*F(0,0)*F(2,0)*pow(F(2,2),2) + 2*d[20]*F(0,0)*F(2,0)*pow(F(2,2),2) +
    d[8]*F(0,1)*F(2,0)*pow(F(2,2),2) + 2*d[19]*F(0,1)*F(2,0)*pow(F(2,2),2) +
    3*d[17]*F(0,2)*F(2,0)*pow(F(2,2),2) + d[8]*F(0,0)*F(2,1)*pow(F(2,2),2) +
    2*d[19]*F(0,0)*F(2,1)*pow(F(2,2),2) + d[4]*F(0,1)*F(2,1)*pow(F(2,2),2) +
    2*d[14]*F(0,1)*F(2,1)*pow(F(2,2),2) + 3*d[12]*F(0,2)*F(2,1)*pow(F(2,2),2) +
    d[17]*F(0,0)*pow(F(2,2),3) + d[12]*F(0,1)*pow(F(2,2),3) +
    d[5]*F(0,2)*pow(F(2,2),3);
    
    c.d[18] = d[0]*pow(F(0,0),2)*F(1,0)*F(2,0) + d[9]*pow(F(0,1),2)*F(1,0)*F(2,0) +
    2*d[15]*F(0,0)*F(0,2)*F(1,0)*F(2,0) + 2*d[18]*F(0,1)*F(0,2)*F(1,0)*F(2,0) +
    d[20]*pow(F(0,2),2)*F(1,0)*F(2,0) + d[1]*F(0,0)*F(0,1)*F(1,1)*F(2,0) +
    d[9]*F(0,0)*F(0,1)*F(1,1)*F(2,0) + d[7]*pow(F(0,1),2)*F(1,1)*F(2,0) +
    d[10]*F(0,0)*F(0,2)*F(1,1)*F(2,0) + d[18]*F(0,0)*F(0,2)*F(1,1)*F(2,0) +
    d[13]*F(0,1)*F(0,2)*F(1,1)*F(2,0) + d[16]*F(0,1)*F(0,2)*F(1,1)*F(2,0) +
    d[19]*pow(F(0,2),2)*F(1,1)*F(2,0) + d[15]*pow(F(0,0),2)*F(1,2)*F(2,0) +
    d[10]*F(0,0)*F(0,1)*F(1,2)*F(2,0) + d[18]*F(0,0)*F(0,1)*F(1,2)*F(2,0) +
    d[13]*pow(F(0,1),2)*F(1,2)*F(2,0) + d[3]*F(0,0)*F(0,2)*F(1,2)*F(2,0) +
    d[20]*F(0,0)*F(0,2)*F(1,2)*F(2,0) + d[8]*F(0,1)*F(0,2)*F(1,2)*F(2,0) +
    d[19]*F(0,1)*F(0,2)*F(1,2)*F(2,0) + d[17]*pow(F(0,2),2)*F(1,2)*F(2,0) +
    d[1]*F(0,0)*F(0,1)*F(1,0)*F(2,1) + d[9]*F(0,0)*F(0,1)*F(1,0)*F(2,1) +
    d[7]*pow(F(0,1),2)*F(1,0)*F(2,1) + d[10]*F(0,0)*F(0,2)*F(1,0)*F(2,1) +
    d[18]*F(0,0)*F(0,2)*F(1,0)*F(2,1) + d[13]*F(0,1)*F(0,2)*F(1,0)*F(2,1) +
    d[16]*F(0,1)*F(0,2)*F(1,0)*F(2,1) + d[19]*pow(F(0,2),2)*F(1,0)*F(2,1) +
    d[9]*pow(F(0,0),2)*F(1,1)*F(2,1) + 2*d[7]*F(0,0)*F(0,1)*F(1,1)*F(2,1) +
    d[2]*pow(F(0,1),2)*F(1,1)*F(2,1) + 2*d[13]*F(0,0)*F(0,2)*F(1,1)*F(2,1) +
    2*d[11]*F(0,1)*F(0,2)*F(1,1)*F(2,1) + d[14]*pow(F(0,2),2)*F(1,1)*F(2,1) +
    d[18]*pow(F(0,0),2)*F(1,2)*F(2,1) + d[13]*F(0,0)*F(0,1)*F(1,2)*F(2,1) +
    d[16]*F(0,0)*F(0,1)*F(1,2)*F(2,1) + d[11]*pow(F(0,1),2)*F(1,2)*F(2,1) +
    d[8]*F(0,0)*F(0,2)*F(1,2)*F(2,1) + d[19]*F(0,0)*F(0,2)*F(1,2)*F(2,1) +
    d[4]*F(0,1)*F(0,2)*F(1,2)*F(2,1) + d[14]*F(0,1)*F(0,2)*F(1,2)*F(2,1) +
    d[12]*pow(F(0,2),2)*F(1,2)*F(2,1) +
    d[6]*F(0,0)*(2*F(0,1)*F(1,0)*F(2,0) + F(0,0)*(F(1,1)*F(2,0) + F(1,0)*F(2,1))) +
    d[15]*pow(F(0,0),2)*F(1,0)*F(2,2) + d[10]*F(0,0)*F(0,1)*F(1,0)*F(2,2) +
    d[18]*F(0,0)*F(0,1)*F(1,0)*F(2,2) + d[13]*pow(F(0,1),2)*F(1,0)*F(2,2) +
    d[3]*F(0,0)*F(0,2)*F(1,0)*F(2,2) + d[20]*F(0,0)*F(0,2)*F(1,0)*F(2,2) +
    d[8]*F(0,1)*F(0,2)*F(1,0)*F(2,2) + d[19]*F(0,1)*F(0,2)*F(1,0)*F(2,2) +
    d[17]*pow(F(0,2),2)*F(1,0)*F(2,2) + d[18]*pow(F(0,0),2)*F(1,1)*F(2,2) +
    d[13]*F(0,0)*F(0,1)*F(1,1)*F(2,2) + d[16]*F(0,0)*F(0,1)*F(1,1)*F(2,2) +
    d[11]*pow(F(0,1),2)*F(1,1)*F(2,2) + d[8]*F(0,0)*F(0,2)*F(1,1)*F(2,2) +
    d[19]*F(0,0)*F(0,2)*F(1,1)*F(2,2) + d[4]*F(0,1)*F(0,2)*F(1,1)*F(2,2) +
    d[14]*F(0,1)*F(0,2)*F(1,1)*F(2,2) + d[12]*pow(F(0,2),2)*F(1,1)*F(2,2) +
    d[20]*pow(F(0,0),2)*F(1,2)*F(2,2) + 2*d[19]*F(0,0)*F(0,1)*F(1,2)*F(2,2) +
    d[14]*pow(F(0,1),2)*F(1,2)*F(2,2) + 2*d[17]*F(0,0)*F(0,2)*F(1,2)*F(2,2) +
    2*d[12]*F(0,1)*F(0,2)*F(1,2)*F(2,2) + d[5]*pow(F(0,2),2)*F(1,2)*F(2,2);
    
    c.d[19] = d[0]*F(0,0)*F(1,0)*pow(F(2,0),2) + d[15]*F(0,2)*F(1,0)*pow(F(2,0),2) +
    d[9]*F(0,1)*F(1,1)*pow(F(2,0),2) + d[18]*F(0,2)*F(1,1)*pow(F(2,0),2) +
    d[15]*F(0,0)*F(1,2)*pow(F(2,0),2) + d[18]*F(0,1)*F(1,2)*pow(F(2,0),2) +
    d[20]*F(0,2)*F(1,2)*pow(F(2,0),2) + d[1]*F(0,1)*F(1,0)*F(2,0)*F(2,1) +
    d[9]*F(0,1)*F(1,0)*F(2,0)*F(2,1) + d[10]*F(0,2)*F(1,0)*F(2,0)*F(2,1) +
    d[18]*F(0,2)*F(1,0)*F(2,0)*F(2,1) + d[1]*F(0,0)*F(1,1)*F(2,0)*F(2,1) +
    d[9]*F(0,0)*F(1,1)*F(2,0)*F(2,1) + 2*d[7]*F(0,1)*F(1,1)*F(2,0)*F(2,1) +
    d[13]*F(0,2)*F(1,1)*F(2,0)*F(2,1) + d[16]*F(0,2)*F(1,1)*F(2,0)*F(2,1) +
    d[10]*F(0,0)*F(1,2)*F(2,0)*F(2,1) + d[18]*F(0,0)*F(1,2)*F(2,0)*F(2,1) +
    d[13]*F(0,1)*F(1,2)*F(2,0)*F(2,1) + d[16]*F(0,1)*F(1,2)*F(2,0)*F(2,1) +
    2*d[19]*F(0,2)*F(1,2)*F(2,0)*F(2,1) + d[9]*F(0,0)*F(1,0)*pow(F(2,1),2) +
    d[7]*F(0,1)*F(1,0)*pow(F(2,1),2) + d[13]*F(0,2)*F(1,0)*pow(F(2,1),2) +
    d[7]*F(0,0)*F(1,1)*pow(F(2,1),2) + d[2]*F(0,1)*F(1,1)*pow(F(2,1),2) +
    d[11]*F(0,2)*F(1,1)*pow(F(2,1),2) + d[13]*F(0,0)*F(1,2)*pow(F(2,1),2) +
    d[11]*F(0,1)*F(1,2)*pow(F(2,1),2) + d[14]*F(0,2)*F(1,2)*pow(F(2,1),2) +
    d[6]*F(2,0)*(F(0,1)*F(1,0)*F(2,0) + F(0,0)*(F(1,1)*F(2,0) + 2*F(1,0)*F(2,1))) +
    2*d[15]*F(0,0)*F(1,0)*F(2,0)*F(2,2) + d[10]*F(0,1)*F(1,0)*F(2,0)*F(2,2) +
    d[18]*F(0,1)*F(1,0)*F(2,0)*F(2,2) + d[3]*F(0,2)*F(1,0)*F(2,0)*F(2,2) +
    d[20]*F(0,2)*F(1,0)*F(2,0)*F(2,2) + d[10]*F(0,0)*F(1,1)*F(2,0)*F(2,2) +
    d[18]*F(0,0)*F(1,1)*F(2,0)*F(2,2) + 2*d[13]*F(0,1)*F(1,1)*F(2,0)*F(2,2) +
    d[8]*F(0,2)*F(1,1)*F(2,0)*F(2,2) + d[19]*F(0,2)*F(1,1)*F(2,0)*F(2,2) +
    d[3]*F(0,0)*F(1,2)*F(2,0)*F(2,2) + d[20]*F(0,0)*F(1,2)*F(2,0)*F(2,2) +
    d[8]*F(0,1)*F(1,2)*F(2,0)*F(2,2) + d[19]*F(0,1)*F(1,2)*F(2,0)*F(2,2) +
    2*d[17]*F(0,2)*F(1,2)*F(2,0)*F(2,2) + 2*d[18]*F(0,0)*F(1,0)*F(2,1)*F(2,2) +
    d[13]*F(0,1)*F(1,0)*F(2,1)*F(2,2) + d[16]*F(0,1)*F(1,0)*F(2,1)*F(2,2) +
    d[8]*F(0,2)*F(1,0)*F(2,1)*F(2,2) + d[19]*F(0,2)*F(1,0)*F(2,1)*F(2,2) +
    d[13]*F(0,0)*F(1,1)*F(2,1)*F(2,2) + d[16]*F(0,0)*F(1,1)*F(2,1)*F(2,2) +
    2*d[11]*F(0,1)*F(1,1)*F(2,1)*F(2,2) + d[4]*F(0,2)*F(1,1)*F(2,1)*F(2,2) +
    d[14]*F(0,2)*F(1,1)*F(2,1)*F(2,2) + d[8]*F(0,0)*F(1,2)*F(2,1)*F(2,2) +
    d[19]*F(0,0)*F(1,2)*F(2,1)*F(2,2) + d[4]*F(0,1)*F(1,2)*F(2,1)*F(2,2) +
    d[14]*F(0,1)*F(1,2)*F(2,1)*F(2,2) + 2*d[12]*F(0,2)*F(1,2)*F(2,1)*F(2,2) +
    d[20]*F(0,0)*F(1,0)*pow(F(2,2),2) + d[19]*F(0,1)*F(1,0)*pow(F(2,2),2) +
    d[17]*F(0,2)*F(1,0)*pow(F(2,2),2) + d[19]*F(0,0)*F(1,1)*pow(F(2,2),2) +
    d[14]*F(0,1)*F(1,1)*pow(F(2,2),2) + d[12]*F(0,2)*F(1,1)*pow(F(2,2),2) +
    d[17]*F(0,0)*F(1,2)*pow(F(2,2),2) + d[12]*F(0,1)*F(1,2)*pow(F(2,2),2) +
    d[5]*F(0,2)*F(1,2)*pow(F(2,2),2);
    
    c.d[20] = d[0]*pow(F(0,0),2)*pow(F(2,0),2) + d[9]*pow(F(0,1),2)*pow(F(2,0),2) +
    2*d[15]*F(0,0)*F(0,2)*pow(F(2,0),2) + 2*d[18]*F(0,1)*F(0,2)*pow(F(2,0),2) +
    d[20]*pow(F(0,2),2)*pow(F(2,0),2) + 2*d[1]*F(0,0)*F(0,1)*F(2,0)*F(2,1) +
    2*d[9]*F(0,0)*F(0,1)*F(2,0)*F(2,1) + 2*d[7]*pow(F(0,1),2)*F(2,0)*F(2,1) +
    2*d[10]*F(0,0)*F(0,2)*F(2,0)*F(2,1) + 2*d[18]*F(0,0)*F(0,2)*F(2,0)*F(2,1) +
    2*d[13]*F(0,1)*F(0,2)*F(2,0)*F(2,1) + 2*d[16]*F(0,1)*F(0,2)*F(2,0)*F(2,1) +
    2*d[19]*pow(F(0,2),2)*F(2,0)*F(2,1) + d[9]*pow(F(0,0),2)*pow(F(2,1),2) +
    2*d[7]*F(0,0)*F(0,1)*pow(F(2,1),2) + d[2]*pow(F(0,1),2)*pow(F(2,1),2) +
    2*d[13]*F(0,0)*F(0,2)*pow(F(2,1),2) + 2*d[11]*F(0,1)*F(0,2)*pow(F(2,1),2) +
    d[14]*pow(F(0,2),2)*pow(F(2,1),2) +
    2*d[6]*F(0,0)*F(2,0)*(F(0,1)*F(2,0) + F(0,0)*F(2,1)) +
    2*d[15]*pow(F(0,0),2)*F(2,0)*F(2,2) + 2*d[10]*F(0,0)*F(0,1)*F(2,0)*F(2,2) +
    2*d[18]*F(0,0)*F(0,1)*F(2,0)*F(2,2) + 2*d[13]*pow(F(0,1),2)*F(2,0)*F(2,2) +
    2*d[3]*F(0,0)*F(0,2)*F(2,0)*F(2,2) + 2*d[20]*F(0,0)*F(0,2)*F(2,0)*F(2,2) +
    2*d[8]*F(0,1)*F(0,2)*F(2,0)*F(2,2) + 2*d[19]*F(0,1)*F(0,2)*F(2,0)*F(2,2) +
    2*d[17]*pow(F(0,2),2)*F(2,0)*F(2,2) + 2*d[18]*pow(F(0,0),2)*F(2,1)*F(2,2) +
    2*d[13]*F(0,0)*F(0,1)*F(2,1)*F(2,2) + 2*d[16]*F(0,0)*F(0,1)*F(2,1)*F(2,2) +
    2*d[11]*pow(F(0,1),2)*F(2,1)*F(2,2) + 2*d[8]*F(0,0)*F(0,2)*F(2,1)*F(2,2) +
    2*d[19]*F(0,0)*F(0,2)*F(2,1)*F(2,2) + 2*d[4]*F(0,1)*F(0,2)*F(2,1)*F(2,2) +
    2*d[14]*F(0,1)*F(0,2)*F(2,1)*F(2,2) + 2*d[12]*pow(F(0,2),2)*F(2,1)*F(2,2) +
    d[20]*pow(F(0,0),2)*pow(F(2,2),2) + 2*d[19]*F(0,0)*F(0,1)*pow(F(2,2),2) +
    d[14]*pow(F(0,1),2)*pow(F(2,2),2) + 2*d[17]*F(0,0)*F(0,2)*pow(F(2,2),2) +
    2*d[12]*F(0,1)*F(0,2)*pow(F(2,2),2) + d[5]*pow(F(0,2),2)*pow(F(2,2),2);
    
    return c;
}
