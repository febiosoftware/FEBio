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
inline double tens3ds::operator()(int i, int j, int k) const
{
	const int LUT[3][3][3] = {
		{{0,1,2},{1,3,4},{2,4,5}},
		{{1,3,4},{3,6,7},{4,7,8}},
		{{2,4,5},{4,7,8},{5,8,9}}};

	return d[LUT[i][j][k]];
}

// contract the right two legs by the dyad formed by a vector  xi = Tijk*Xi*Xk
inline vec3d tens3ds::contractdyad1(const vec3d& v)
{
    vec3d x;
	x.x = d[0]*v.x*v.x + 2*d[1]*v.x*v.y + 2*d[2]*v.x*v.z + d[3]*v.y*v.y + 2*d[4]*v.y*v.z + d[5]*v.z*v.z;
	x.y = d[1]*v.x*v.x + 2*d[3]*v.x*v.y + 2*d[4]*v.x*v.z + d[6]*v.y*v.y + 2*d[7]*v.y*v.z + d[8]*v.z*v.z;
	x.z = d[2]*v.x*v.x + 2*d[4]*v.x*v.y + 2*d[5]*v.x*v.z + d[7]*v.y*v.y + 2*d[8]*v.y*v.z + d[9]*v.z*v.z;

	return x;
}

// triple contraction by a similar 3o tensor m = Tijk*Hijk
inline double tens3ds::tripledot(const tens3ds& H)
{
	const double* h = H.d;
	return d[0]*h[0] + 3*d[1]*h[1] + 3*d[2]*h[2] + 3*d[3]*h[3] + 6*d[4]*h[4] + 3*d[5]*h[5] + d[6]*h[6] + 3*d[7]*h[7] + 3*d[8]*h[8] + d[9]*h[9];
}

// calculates the symmetric tensor A_ijk = (l_i*m_j*r_k + perm(i,j,k))/6
inline tens3ds dyad3s(const vec3d& l, const vec3d& m, const vec3d& r)
{
	tens3ds a;
	a.d[0] = (l.x*m.x*r.x); 
	a.d[1] = (l.x*m.x*r.y + l.x*m.y*r.x + l.y*m.x*r.x)/3.0; 
	a.d[2] = (l.x*m.x*r.z + l.x*m.z*r.x + l.z*m.x*r.x)/3.0;
	a.d[3] = (l.x*m.y*r.y + l.y*m.x*r.y + l.y*m.y*r.x)/3.0; 
	a.d[4] = (l.x*m.y*r.z + l.y*m.x*r.z + l.z*m.y*r.x + l.x*m.z*r.y + l.z*m.x*r.y + l.y*m.z*r.x)/6.0; 
	a.d[5] = (l.x*m.z*r.z + l.z*m.x*r.z + l.z*m.z*r.x)/3.0;
	a.d[6] = (l.y*m.y*r.y); 
	a.d[7] = (l.y*m.y*r.z + l.y*m.z*r.y + l.z*m.y*r.y)/3.0;
	a.d[8] = (l.y*m.z*r.z + l.z*m.y*r.z + l.z*m.z*r.y)/3.0;
	a.d[9] = (l.z*m.z*r.z);
	return a;
}
