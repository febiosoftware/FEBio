/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, Columbia University, and others.

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
