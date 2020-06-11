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

#include "mat3d.h"

// access operator
inline double tens3drs::operator () (int i, int j, int k) const
{
	const int m[3][3] = {{0,1,2},{1,3,4},{2,4,5}};
	return d[6*i + m[j][k] ];
}

// access operator
inline double& tens3drs::operator () (int i, int j, int k)
{
	const int m[3][3] = {{0,1,2},{1,3,4},{2,4,5}};
	return d[6*i + m[j][k] ];
}

inline tens3drs::tens3drs(double a)
{
	int nnz = tensor_traits<tens3drs>::NNZ;
	for (int i = 0; i < nnz; ++i) d[i] = a;
}

// contract the right two legs by the dyad formed by a vector  xi = Gijk*Xj*Xk
inline vec3d tens3drs::contractdyad1(const vec3d& v) const
{
    vec3d x;
	x.x = d[ 0]*v.x*v.x + 2*d[ 1]*v.x*v.y + 2*d[ 2]*v.x*v.z + d[ 3]*v.y*v.y + 2*d[ 4]*v.y*v.z + d[ 5]*v.z*v.z;
	x.y = d[ 6]*v.x*v.x + 2*d[ 7]*v.x*v.y + 2*d[ 8]*v.x*v.z + d[ 9]*v.y*v.y + 2*d[10]*v.y*v.z + d[11]*v.z*v.z;
	x.z = d[12]*v.x*v.x + 2*d[13]*v.x*v.y + 2*d[14]*v.x*v.z + d[15]*v.y*v.y + 2*d[16]*v.y*v.z + d[17]*v.z*v.z;

	return x;
}

// contract the right two legs by a symmetric 2o tensor  xi = Gijk*Sjk
inline vec3d tens3drs::contract2s(const mat3ds& s) const
{
    vec3d x;
	x.x = d[ 0]*s.xx() + 2*d[ 1]*s.xy() + 2*d[ 2]*s.xz() + d[ 3]*s.yy() + 2*d[ 4]*s.yz() + d[ 5]*s.zz();
	x.y = d[ 6]*s.xx() + 2*d[ 7]*s.xy() + 2*d[ 8]*s.xz() + d[ 9]*s.yy() + 2*d[10]*s.yz() + d[11]*s.zz();
	x.z = d[12]*s.xx() + 2*d[13]*s.xy() + 2*d[14]*s.xz() + d[15]*s.yy() + 2*d[16]*s.yz() + d[17]*s.zz();

	return x;
}

// triple contraction by a similar 3o tensor m = Gijk*Hijk
inline double tens3drs::tripledot(const tens3drs& H) const
{
	const double* h = H.d;
	return	d[ 0]*h[ 0] + 2*d[ 1]*h[ 1] + 2*d[ 2]*h[ 2] + d[ 3]*h[ 3] + 2*d[ 4]*h[ 4] + d[ 5]*h[ 5] 
		  + d[ 6]*h[ 6] + 2*d[ 7]*h[ 7] + 2*d[ 8]*h[ 8] + d[ 9]*h[ 9] + 2*d[10]*h[10] + d[11]*h[11]
		  + d[12]*h[12] + 2*d[13]*h[13] + 2*d[14]*h[14] + d[15]*h[15] + 2*d[16]*h[16] + d[17]*h[17];  
}

// contract the right two legs by the dyad formed by a vector  xi = Gijk*Vj*Wk
inline vec3d tens3drs::contractdyad2(const vec3d& v, const vec3d& w)
{
    vec3d x;
	x.x = d[ 0]*v.x*w.x + d[ 1]*(v.x*w.y + v.y*w.x) + d[ 2]*(v.x*w.z + v.z*w.x) + d[ 3]*v.y*w.y + d[ 4]*(v.y*w.z + v.z*w.y) + d[ 5]*v.z*w.z;
	x.y = d[ 6]*v.x*w.x + d[ 7]*(v.x*w.y + v.y*w.x) + d[ 8]*(v.x*w.z + v.z*w.x) + d[ 9]*v.y*w.y + d[10]*(v.y*w.z + v.z*w.y) + d[11]*v.z*w.z;
	x.z = d[12]*v.x*w.x + d[13]*(v.x*w.y + v.y*w.x) + d[14]*(v.x*w.z + v.z*w.x) + d[15]*v.y*w.y + d[16]*(v.y*w.z + v.z*w.y) + d[17]*v.z*w.z;

	return x;
}

// calculates the dyadic product T_ijk = l_i*r_j*r_k
inline tens3drs dyad3rs(const vec3d& l, const vec3d& r)
{
	tens3drs a;
	a.d[ 0] = l.x*r.x*r.x; 
	a.d[ 1] = l.x*r.x*r.y;
	a.d[ 2] = l.x*r.x*r.z;
	a.d[ 3] = l.x*r.y*r.y;
	a.d[ 4] = l.x*r.y*r.z;
	a.d[ 5] = l.x*r.z*r.z;
	a.d[ 6] = l.y*r.x*r.x;
	a.d[ 7] = l.y*r.x*r.y;
	a.d[ 8] = l.y*r.x*r.z;
	a.d[ 9] = l.y*r.y*r.y;
	a.d[10] = l.y*r.y*r.z;
	a.d[11] = l.y*r.z*r.z;
	a.d[12] = l.z*r.x*r.x;
	a.d[13] = l.z*r.x*r.y;
	a.d[14] = l.z*r.x*r.z;
	a.d[15] = l.z*r.y*r.y;
	a.d[16] = l.z*r.y*r.z;
	a.d[17] = l.z*r.z*r.z;
	return a;
}

// calculates the dyadic product T_ijk = 1/2*(L_ij*r_k + L_ik*r_j)
inline tens3drs dyad3rs(const mat3d& L, const vec3d& r)
{
	tens3drs a;
	a.d[0] =  L(0, 0)*r.x;
	a.d[1] = (L(0, 0)*r.y + L(0, 1)*r.x)*0.5;
	a.d[2] = (L(0, 0)*r.z + L(0, 2)*r.x)*0.5;
	a.d[3] =  L(0, 1)*r.y;
	a.d[4] = (L(0, 1)*r.z + L(0, 2)*r.y)*0.5;
	a.d[5] =  L(0, 2)*r.z;

	a.d[ 6] =  L(1, 0)*r.x;
	a.d[ 7] = (L(1, 0)*r.y + L(1, 1)*r.x)*0.5;
	a.d[ 8] = (L(1, 0)*r.z + L(1, 2)*r.x)*0.5;
	a.d[ 9] =  L(1, 1)*r.y;
	a.d[10] = (L(1, 1)*r.z + L(1, 2)*r.y)*0.5;
	a.d[11] =  L(1, 2)*r.z;

	a.d[12] =  L(2, 0)*r.x;
	a.d[13] = (L(2, 0)*r.y + L(2, 1)*r.x)*0.5;
	a.d[14] = (L(2, 0)*r.z + L(2, 2)*r.x)*0.5;
	a.d[15] =  L(2, 1)*r.y;
	a.d[16] = (L(2, 1)*r.z + L(2, 2)*r.y)*0.5;
	a.d[17] =  L(2, 2)*r.z;

	return a;
}

// calculate the transpose ((G_iJK)T = G_KJi)
inline tens3dls tens3drs::transpose()
{
	tens3dls GLC;

	GLC.d[ 0] = d[ 0];
	GLC.d[ 3] = d[ 1];
	GLC.d[ 6] = d[ 2];
	GLC.d[ 9] = d[ 3];
	GLC.d[12] = d[ 4];
	GLC.d[15] = d[ 5];
	GLC.d[ 1] = d[ 6];
	GLC.d[ 4] = d[ 7];
	GLC.d[ 7] = d[ 8];
	GLC.d[10] = d[ 9];
	GLC.d[13] = d[10];
	GLC.d[16] = d[11];
	GLC.d[ 2] = d[12];
	GLC.d[ 5] = d[13];
	GLC.d[ 8] = d[14];
	GLC.d[11] = d[15];
	GLC.d[14] = d[16];
	GLC.d[17] = d[17];

	return GLC;
}

// contract each leg by a 2o tensor (intended to calculate the inverse deformation hessian according to Finv_Ii * G_iJK * Finv_Jj * Fin_Kk)
inline void tens3drs::contractleg2(const mat3d& F, int leg)
{
	tens3drs G = *this;
	
	if (leg == 1)
	{
		d[0] = F(0,0)*G.d[0] + F(0,1)*G.d[6] + F(0,2)*G.d[12];
		d[1] = F(0,0)*G.d[1] + F(0,1)*G.d[7] + F(0,2)*G.d[13];
		d[2] = F(0,0)*G.d[2] + F(0,1)*G.d[8] + F(0,2)*G.d[14];
		d[3] = F(0,0)*G.d[3] + F(0,1)*G.d[9] + F(0,2)*G.d[15];
		d[4] = F(0,0)*G.d[4] + F(0,1)*G.d[10] + F(0,2)*G.d[16];
		d[5] = F(0,0)*G.d[5] + F(0,1)*G.d[11] + F(0,2)*G.d[17];
		d[6] = F(1,0)*G.d[0] + F(1,1)*G.d[6] + F(1,2)*G.d[12];
		d[7] = F(1,0)*G.d[1] + F(1,1)*G.d[7] + F(1,2)*G.d[13];
		d[8] = F(1,0)*G.d[2] + F(1,1)*G.d[8] + F(1,2)*G.d[14];
		d[9] = F(1,0)*G.d[3] + F(1,1)*G.d[9] + F(1,2)*G.d[15];
		d[10] = F(1,0)*G.d[4] + F(1,1)*G.d[10] + F(1,2)*G.d[16];
		d[11] = F(1,0)*G.d[5] + F(1,1)*G.d[11] + F(1,2)*G.d[17];
		d[12] = F(2,0)*G.d[0] + F(2,1)*G.d[6] + F(2,2)*G.d[12];
		d[13] = F(2,0)*G.d[1] + F(2,1)*G.d[7] + F(2,2)*G.d[13];
		d[14] = F(2,0)*G.d[2] + F(2,1)*G.d[8] + F(2,2)*G.d[14];
		d[15] = F(2,0)*G.d[3] + F(2,1)*G.d[9] + F(2,2)*G.d[15];
		d[16] = F(2,0)*G.d[4] + F(2,1)*G.d[10] + F(2,2)*G.d[16];
		d[17] = F(2,0)*G.d[5] + F(2,1)*G.d[11] + F(2,2)*G.d[17];
	}
	else if (leg == 2)
	{
		d[0] = G.d[0]*F(0,0) + G.d[1]*F(1,0) + G.d[2]*F(2,0);
		d[1] = G.d[1]*F(0,0) + G.d[3]*F(1,0) + G.d[4]*F(2,0);
		d[2] = G.d[2]*F(0,0) + G.d[4]*F(1,0) + G.d[5]*F(2,0);
		d[3] = G.d[1]*F(0,1) + G.d[3]*F(1,1) + G.d[4]*F(2,1);
		d[4] = G.d[2]*F(0,1) + G.d[4]*F(1,1) + G.d[5]*F(2,1);
		d[5] = G.d[2]*F(0,2) + G.d[4]*F(1,2) + G.d[5]*F(2,2);
		d[6] = G.d[6]*F(0,0) + G.d[7]*F(1,0) + G.d[8]*F(2,0);
		d[7] = G.d[7]*F(0,0) + G.d[9]*F(1,0) + G.d[10]*F(2,0);
		d[8] = G.d[8]*F(0,0) + G.d[10]*F(1,0) + G.d[11]*F(2,0);
		d[9] = G.d[7]*F(0,1) + G.d[9]*F(1,1) + G.d[10]*F(2,1);
		d[10] = G.d[8]*F(0,1) + G.d[10]*F(1,1) + G.d[11]*F(2,1);
		d[11] = G.d[8]*F(0,2) + G.d[10]*F(1,2) + G.d[11]*F(2,2);
		d[12] = G.d[12]*F(0,0) + G.d[13]*F(1,0) + G.d[4]*F(2,0);
		d[13] = G.d[13]*F(0,0) + G.d[14]*F(1,0) + G.d[16]*F(2,0);
		d[14] = G.d[14]*F(0,0) + G.d[16]*F(1,0) + G.d[17]*F(2,0);
		d[15] = G.d[13]*F(0,1) + G.d[15]*F(1,1) + G.d[16]*F(2,1);
		d[16] = G.d[14]*F(0,1) + G.d[16]*F(1,1) + G.d[17]*F(2,1);
		d[17] = G.d[14]*F(0,2) + G.d[16]*F(1,2) + G.d[17]*F(2,2);
	}
	else if (leg == 3)
	{
		d[0] = G.d[0]*F(0,0) + G.d[1]*F(1,0) + G.d[2]*F(2,0);
		d[1] = G.d[0]*F(0,1) + G.d[1]*F(1,1) + G.d[2]*F(2,1);
		d[2] = G.d[0]*F(0,2) + G.d[1]*F(1,2) + G.d[2]*F(2,2);
		d[3] = G.d[1]*F(0,1) + G.d[3]*F(1,1) + G.d[4]*F(2,1);
		d[4] = G.d[1]*F(0,2) + G.d[3]*F(1,2) + G.d[4]*F(2,2);
		d[5] = G.d[2]*F(0,2) + G.d[4]*F(1,2) + G.d[5]*F(2,2);
		d[6] = G.d[6]*F(0,0) + G.d[7]*F(1,0) + G.d[8]*F(2,0);
		d[7] = G.d[6]*F(0,1) + G.d[7]*F(1,1) + G.d[8]*F(2,1);
		d[8] = G.d[6]*F(0,2) + G.d[7]*F(1,2) + G.d[8]*F(2,2);
		d[9] = G.d[7]*F(0,1) + G.d[9]*F(1,1) + G.d[10]*F(2,1);
		d[10] = G.d[7]*F(0,2) + G.d[9]*F(1,2) + G.d[10]*F(2,2);
		d[11] = G.d[8]*F(0,2) + G.d[10]*F(1,2) + G.d[11]*F(2,2);
		d[12] = G.d[12]*F(0,0) + G.d[13]*F(1,0) + G.d[14]*F(2,0);
		d[13] = G.d[12]*F(0,1) + G.d[13]*F(1,1) + G.d[14]*F(2,1);
		d[14] = G.d[12]*F(0,2) + G.d[13]*F(1,2) + G.d[14]*F(2,2);
		d[15] = G.d[13]*F(0,1) + G.d[15]*F(1,1) + G.d[16]*F(2,1);
		d[16] = G.d[13]*F(0,2) + G.d[15]*F(1,2) + G.d[16]*F(2,2);
		d[17] = G.d[14]*F(0,2) + G.d[16]*F(1,2) + G.d[17]*F(2,2);
	}
}

// multiply by a 2o tensor on the left (F_Ii * G_iJK)
inline tens3drs operator * (const mat3d& F, const tens3drs& t)
{
	tens3drs G;

	const double* d = t.d;
	G.d[ 0] = F(0,0)*d[0] + F(0,1)*d[ 6] + F(0,2)*d[12];
	G.d[ 1] = F(0,0)*d[1] + F(0,1)*d[ 7] + F(0,2)*d[13];
	G.d[ 2] = F(0,0)*d[2] + F(0,1)*d[ 8] + F(0,2)*d[14];
	G.d[ 3] = F(0,0)*d[3] + F(0,1)*d[ 9] + F(0,2)*d[15];
	G.d[ 4] = F(0,0)*d[4] + F(0,1)*d[10] + F(0,2)*d[16];
	G.d[ 5] = F(0,0)*d[5] + F(0,1)*d[11] + F(0,2)*d[17];

	G.d[ 6] = F(1,0)*d[0] + F(1,1)*d[ 6] + F(1,2)*d[12];
	G.d[ 7] = F(1,0)*d[1] + F(1,1)*d[ 7] + F(1,2)*d[13];
	G.d[ 8] = F(1,0)*d[2] + F(1,1)*d[ 8] + F(1,2)*d[14];
	G.d[ 9] = F(1,0)*d[3] + F(1,1)*d[ 9] + F(1,2)*d[15];
	G.d[10] = F(1,0)*d[4] + F(1,1)*d[10] + F(1,2)*d[16];
	G.d[11] = F(1,0)*d[5] + F(1,1)*d[11] + F(1,2)*d[17];

	G.d[12] = F(2,0)*d[0] + F(2,1)*d[ 6] + F(2,2)*d[12];
	G.d[13] = F(2,0)*d[1] + F(2,1)*d[ 7] + F(2,2)*d[13];
	G.d[14] = F(2,0)*d[2] + F(2,1)*d[ 8] + F(2,2)*d[14];
	G.d[15] = F(2,0)*d[3] + F(2,1)*d[ 9] + F(2,2)*d[15];
	G.d[16] = F(2,0)*d[4] + F(2,1)*d[10] + F(2,2)*d[16];
	G.d[17] = F(2,0)*d[5] + F(2,1)*d[11] + F(2,2)*d[17];

	return G;
}
