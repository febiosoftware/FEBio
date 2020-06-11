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
// NOTE: This file is automatically included from tens3d.h
// Users should not include this file manually!

// calculate the transpose ((G_KJi)T = G_iJK)
inline tens3drs tens3dls::transpose()
{
	tens3drs GRC;

	GRC.d[0] =  d[ 0];
	GRC.d[1] =  d[ 3];
	GRC.d[2] =  d[ 6];
	GRC.d[3] =  d[ 9];
	GRC.d[4] =  d[12];
	GRC.d[5] =  d[15];
	GRC.d[6] =  d[ 1];
	GRC.d[7] =  d[ 4];
	GRC.d[8] =  d[ 7];
	GRC.d[9] =  d[10];
	GRC.d[10] = d[13];
	GRC.d[11] = d[16];
	GRC.d[12] = d[ 2];
	GRC.d[13] = d[ 5];
	GRC.d[14] = d[ 8];
	GRC.d[15] = d[11];
	GRC.d[16] = d[14];
	GRC.d[17] = d[17];

	return GRC;
}

// multiply by a double on right (G_ifk * d)
inline tens3dls tens3dls::operator * (const double& f) const
{
    tens3dls G;
    
    G.d[0] = d[0]*f;
    G.d[1] = d[1]*f;
    G.d[2] = d[2]*f;
    
    G.d[3] = d[3]*f;
    G.d[4] = d[4]*f;
    G.d[5] = d[5]*f;
    
    G.d[6] = d[6]*f;
    G.d[7] = d[7]*f;
    G.d[8] = d[8]*f;
    
    G.d[ 9] = d[9]*f;
    G.d[10] = d[10]*f;
    G.d[11] = d[11]*f;
    
    G.d[12] = d[12]*f;
    G.d[13] = d[13]*f;
    G.d[14] = d[14]*f;
    
    G.d[15] = d[15]*f;
    G.d[16] = d[16]*f;
    G.d[17] = d[17]*f;
    
    return G;
}

// multiply by a 2o tensor on the right (G_KJi * F_iI)
inline tens3dls tens3dls::operator * (const mat3d& F) const
{
	tens3dls G;

	G.d[0] = d[0]*F(0,0) + d[1]*F(1,0) + d[2]*F(2,0);
	G.d[1] = d[0]*F(0,1) + d[1]*F(1,1) + d[2]*F(2,1);
	G.d[2] = d[0]*F(0,2) + d[1]*F(1,2) + d[2]*F(2,2);
	
	G.d[3] = d[3]*F(0,0) + d[4]*F(1,0) + d[5]*F(2,0);
	G.d[4] = d[3]*F(0,1) + d[4]*F(1,1) + d[5]*F(2,1);
	G.d[5] = d[3]*F(0,2) + d[4]*F(1,2) + d[5]*F(2,2);

	G.d[6] = d[6]*F(0,0) + d[7]*F(1,0) + d[8]*F(2,0);
	G.d[7] = d[6]*F(0,1) + d[7]*F(1,1) + d[8]*F(2,1);
	G.d[8] = d[6]*F(0,2) + d[7]*F(1,2) + d[8]*F(2,2);

	G.d[ 9] = d[9]*F(0,0) + d[10]*F(1,0) + d[11]*F(2,0);
	G.d[10] = d[9]*F(0,1) + d[10]*F(1,1) + d[11]*F(2,1);
	G.d[11] = d[9]*F(0,2) + d[10]*F(1,2) + d[11]*F(2,2);

	G.d[12] = d[12]*F(0,0) + d[13]*F(1,0) + d[14]*F(2,0);
	G.d[13] = d[12]*F(0,1) + d[13]*F(1,1) + d[14]*F(2,1);
	G.d[14] = d[12]*F(0,2) + d[13]*F(1,2) + d[14]*F(2,2);

	G.d[15] = d[15]*F(0,0) + d[16]*F(1,0) + d[17]*F(2,0);
	G.d[16] = d[15]*F(0,1) + d[16]*F(1,1) + d[17]*F(2,1);
	G.d[17] = d[15]*F(0,2) + d[16]*F(1,2) + d[17]*F(2,2);

	return G;
}

// calculate the trace Tijk -> Tijj
inline vec3d tens3dls::trace()
{
    double a = d[0] + d[4] + d[8];
    double b = d[3] + d[10] + d[14];
    double c = d[6] + d[13] + d[17];
    vec3d v = vec3d(a,b,c);
    
    return v;
}

// generalize tensor from tens3dls to tens3d
// [G] = [G111 G112 G113 G121 G122 G123 G131 G132 G133 G221 G222 G223 G231 G232 G233 G331 G332 G333]
//     =    G0   G1   G2   G3   G4   G5   G6   G7   G8   G9  G10  G11  G12  G13  G14  G15  G16  G17
// [T] = [T111 T112 T113 T121 T122 T123 T131 T132 T133 T211 T212 T213 T221 T222 T223 T231 T232 T233 T311 T312 T313 T321 T322 T323 T331 T332 T333
//     =    T0   T1   T2   T3   T4   T5   T6   T7   T8   T9  T10  T11  T12  T13  T14  T15  T16  T17  T18  T19  T20  T21  T22  T23  T24  T25  T26
inline tens3d tens3dls::generalize()
{
    tens3d G;
    
    G.d[0] = d[0];
    G.d[1] = d[1];
    G.d[2] = d[2];
    G.d[3] = d[3];
    G.d[4] = d[4];
    G.d[5] = d[5];
    G.d[6] = d[6];
    G.d[7] = d[7];
    G.d[8] = d[8];
    G.d[9] = d[3];
    G.d[10] = d[4];
    G.d[11] = d[5];
    G.d[12] = d[9];
    G.d[13] = d[10];
    G.d[14] = d[11];
    G.d[15] = d[12];
    G.d[16] = d[13];
    G.d[17] = d[14];
    G.d[18] = d[6];
    G.d[19] = d[7];
    G.d[20] = d[8];
    G.d[21] = d[12];
    G.d[22] = d[13];
    G.d[23] = d[14];
    G.d[24] = d[15];
    G.d[25] = d[16];
    G.d[26] = d[17];
    
    return G;
}

// calculates the dyadic product T_ijk = 1/2*(L_ij*r_k + L_ji*r_k)
// [G] = [G111 G112 G113 G121 G122 G123 G131 G132 G133 G221 G222 G223 G231 G232 G233 G331 G332 G333]
//     =    G0   G1   G2   G3   G4   G5   G6   G7   G8   G9  G10  G11  G12  G13  G14  G15  G16  G17
inline tens3dls dyad3ls(const mat3ds& L, const vec3d& r)
{
    tens3dls a;
    a.d[ 0] = L(0, 0)*r.x;
    a.d[ 1] = L(0, 0)*r.y;
    a.d[ 2] = L(0, 0)*r.z;
    a.d[ 3] = L(0, 1)*r.x;
    a.d[ 4] = L(0, 1)*r.y;
    a.d[ 5] = L(0, 1)*r.z;
    a.d[ 6] = L(0, 2)*r.x;
    a.d[ 7] = L(0, 2)*r.y;
    a.d[ 8] = L(0, 2)*r.z;
    a.d[ 9] = L(1, 1)*r.x;
    a.d[10] = L(1, 1)*r.y;
    a.d[11] = L(1, 1)*r.z;
    a.d[12] = L(1, 2)*r.x;
    a.d[13] = L(1, 2)*r.y;
    a.d[14] = L(1, 2)*r.z;
    a.d[15] = L(2, 2)*r.x;
    a.d[16] = L(2, 2)*r.y;
    a.d[17] = L(2, 2)*r.z;
    
    return a;
}
