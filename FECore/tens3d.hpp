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
// NOTE: This file is automatically included from tens3drs.h
// Users should not include this file manually!

// access operator
inline double tens3d::operator()(int i, int j, int k) const
{
	return d[i*9 + j*3 + k];
}

// access operator
inline double& tens3d::operator()(int i, int j, int k)
{
	return d[i*9 + j*3 + k];
}

inline tens3d::tens3d(double a)
{
    int nnz = tensor_traits<tens3d>::NNZ;
    for (int i = 0; i < nnz; ++i) d[i] = a;
}

// symmetrize a general 3o tensor
inline tens3ds tens3d::symm()
{
	tens3ds t;

	t.d[0] =  d[0]; 
	t.d[1] = (d[1] + d[3]  + d[9])/3.; 
	t.d[2] = (d[2] + d[6]  + d[18])/3.;
	t.d[3] = (d[4] + d[10] + d[12])/3.; 
	t.d[4] = (d[5] + d[11] + d[21] + d[7] + d[19] + d[15])/6.; 
	t.d[5] = (d[8] + d[20] + d[24])/3.;
	t.d[6] =  d[13]; 
	t.d[7] = (d[14] + d[16] + d[22])/3.;
	t.d[8] = (d[17] + d[23] + d[25])/3.;
	t.d[9] =  d[26]; 

	return t;
}

// return transpose of right side
// [T] = [T111 T112 T113 T121 T122 T123 T131 T132 T133 T211 T212 T213 T221 T222 T223 T231 T232 T233 T311 T312 T313 T321 T322 T323 T331 T332 T333
//     =    T0   T1   T2   T3   T4   T5   T6   T7   T8   T9  T10  T11  T12  T13  T14  T15  T16  T17  T18  T19  T20  T21  T22  T23  T24  T25  T26
inline tens3d tens3d::transposer()
{
    tens3d s;
    s.d[ 0] = d[ 0];
    s.d[ 1] = d[ 3];
    s.d[ 2] = d[ 6];
    s.d[ 3] = d[ 1];
    s.d[ 4] = d[ 4];
    s.d[ 5] = d[ 7];
    s.d[ 6] = d[ 2];
    s.d[ 7] = d[ 5];
    s.d[ 8] = d[ 8];
    s.d[ 9] = d[ 9];
    s.d[10] = d[12];
    s.d[11] = d[15];
    s.d[12] = d[10];
    s.d[13] = d[13];
    s.d[14] = d[16];
    s.d[15] = d[11];
    s.d[16] = d[14];
    s.d[17] = d[17];
    s.d[18] = d[18];
    s.d[19] = d[21];
    s.d[20] = d[24];
    s.d[21] = d[19];
    s.d[22] = d[22];
    s.d[23] = d[25];
    s.d[24] = d[20];
    s.d[25] = d[23];
    s.d[26] = d[26];
    return s;
}

// contract the right two legs by a 2o tensor  xi = Gijk*Sjk
// [T] = [T111 T112 T113 T121 T122 T123 T131 T132 T133 T211 T212 T213 T221 T222 T223 T231 T232 T233 T311 T312 T313 T321 T322 T323 T331 T332 T333
//     =    T0   T1   T2   T3   T4   T5   T6   T7   T8   T9  T10  T11  T12  T13  T14  T15  T16  T17  T18  T19  T20  T21  T22  T23  T24  T25  T26
inline vec3d tens3d::contract2(const mat3d& s) const
{
    vec3d x;
    x.x = d[0]*s[0][0] + d[1]*s[0][1] + d[2]*s[0][2] + d[3]*s[1][0] + d[4]*s[1][1] + d[5]*s[1][2] + d[6]*s[2][0] + d[7]*s[2][1] + d[8]*s[2][2];
    x.y = d[9]*s[0][0] + d[10]*s[0][1] + d[11]*s[0][2] + d[12]*s[1][0] + d[13]*s[1][1] + d[14]*s[1][2] + d[15]*s[2][0] + d[16]*s[2][1] + d[17]*s[2][2];
    x.z = d[18]*s[0][0] + d[19]*s[0][1] + d[20]*s[0][2] + d[21]*s[1][0] + d[22]*s[1][1] + d[23]*s[1][2] + d[24]*s[2][0] + d[25]*s[2][1] + d[26]*s[2][2];
    
    return x;
}

// contract the right leg by a vector  xij = Gijk*vk
// [T] = [T111 T112 T113 T121 T122 T123 T131 T132 T133 T211 T212 T213 T221 T222 T223 T231 T232 T233 T311 T312 T313 T321 T322 T323 T331 T332 T333
//     =    T0   T1   T2   T3   T4   T5   T6   T7   T8   T9  T10  T11  T12  T13  T14  T15  T16  T17  T18  T19  T20  T21  T22  T23  T24  T25  T26
inline mat3d tens3d::contract1(const vec3d& v) const
{
    mat3d x;
    
    x[0][0] = d[0]*v.x + d[1]*v.y + d[2]*v.z;
    x[0][1] = d[3]*v.x + d[4]*v.y + d[5]*v.z;
    x[0][2] = d[6]*v.x + d[7]*v.y + d[8]*v.z;
    
    x[1][0] = d[9]*v.x + d[10]*v.y + d[11]*v.z;
    x[1][1] = d[12]*v.x + d[13]*v.y + d[14]*v.z;
    x[1][2] = d[15]*v.x + d[16]*v.y + d[17]*v.z;
    
    x[2][0] = d[18]*v.x + d[19]*v.y + d[20]*v.z;
    x[2][1] = d[21]*v.x + d[22]*v.y + d[23]*v.z;
    x[2][2] = d[24]*v.x + d[25]*v.y + d[26]*v.z;
    
    return x;
}

inline tens3d operator + (const tens3dls& l, const tens3drs& r)
{
	tens3d s;
	s.d[ 0] = l.d[ 0] + r.d[ 0];	// S111 = L111 + R111
	s.d[ 1] = l.d[ 1] + r.d[ 1];	// S112 = L112 + R112
	s.d[ 2] = l.d[ 2] + r.d[ 2];	// S113 = L113 + R113
	s.d[ 3] = l.d[ 3] + r.d[ 1];	// S121 = L121 + R112
	s.d[ 4] = l.d[ 4] + r.d[ 3];	// S122 = L122 + R122
	s.d[ 5] = l.d[ 5] + r.d[ 4];	// S123 = L123 + R123
	s.d[ 6] = l.d[ 6] + r.d[ 2];	// S131 = L131 + R113
	s.d[ 7] = l.d[ 7] + r.d[ 4];	// S132 = L132 + R123
	s.d[ 8] = l.d[ 8] + r.d[ 5];	// S133 = L133 + R133
	s.d[ 9] = l.d[ 3] + r.d[ 6];	// S211 = L121 + R211
	s.d[10] = l.d[ 4] + r.d[ 7];	// S212 = L122 + R212
	s.d[11] = l.d[ 5] + r.d[ 8];	// S213 = L123 + R213
	s.d[12] = l.d[ 9] + r.d[ 7];	// S221 = L221 + R212
	s.d[13] = l.d[10] + r.d[ 9];	// S222 = L222 + R222
	s.d[14] = l.d[11] + r.d[10];	// S223 = L223 + R223
	s.d[15] = l.d[12] + r.d[ 8];	// S231 = L231 + R213
	s.d[16] = l.d[13] + r.d[10];	// S232 = L232 + R223
	s.d[17] = l.d[14] + r.d[11];	// S233 = L233 + R233
	s.d[18] = l.d[ 6] + r.d[12];	// S311 = L131 + R311
	s.d[19] = l.d[ 7] + r.d[13];	// S312 = L132 + R312
	s.d[20] = l.d[ 8] + r.d[14];	// S313 = L133 + R313
	s.d[21] = l.d[12] + r.d[13];	// S321 = L231 + R312
	s.d[22] = l.d[13] + r.d[15];	// S322 = L232 + R322
	s.d[23] = l.d[14] + r.d[16];	// S323 = L233 + R323
	s.d[24] = l.d[15] + r.d[14];	// S331 = L331 + R313
	s.d[25] = l.d[16] + r.d[16];	// S332 = L332 + R323
	s.d[26] = l.d[17] + r.d[17];	// S333 = L333 + R333
	return s;
}
