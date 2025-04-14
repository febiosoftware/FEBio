/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
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

inline tens4dmm::tens4dmm(const double g)
{
    d[ 0] = d[ 1] = d[ 2] = d[ 3] = d[ 4] = d[ 5] =
    d[ 6] = d[ 7] = d[ 8] = d[ 9] = d[10] = d[11] =
    d[12] = d[13] = d[14] = d[15] = d[16] = d[17] =
    d[18] = d[19] = d[20] = d[21] = d[22] = d[23] =
    d[24] = d[25] = d[26] = d[27] = d[28] = d[29] =
    d[30] = d[31] = d[32] = d[33] = d[34] = d[35] = g;
}

inline tens4dmm::tens4dmm(const tens4ds t)
{
    d[ 0] = t.d[ 0];    d[ 6] = t.d[ 1];    d[12] = t.d[ 3];    d[18] = t.d[ 6];    d[24] = t.d[10];    d[30] = t.d[15];
    d[ 1] = t.d[ 1];    d[ 7] = t.d[ 2];    d[13] = t.d[ 4];    d[19] = t.d[ 7];    d[25] = t.d[11];    d[31] = t.d[16];
    d[ 2] = t.d[ 3];    d[ 8] = t.d[ 4];    d[14] = t.d[ 5];    d[20] = t.d[ 8];    d[26] = t.d[12];    d[32] = t.d[17];
    d[ 3] = t.d[ 6];    d[ 9] = t.d[ 7];    d[15] = t.d[ 8];    d[21] = t.d[ 9];    d[27] = t.d[13];    d[33] = t.d[18];
    d[ 4] = t.d[10];    d[10] = t.d[11];    d[16] = t.d[12];    d[22] = t.d[13];    d[28] = t.d[14];    d[34] = t.d[19];
    d[ 5] = t.d[15];    d[11] = t.d[16];    d[17] = t.d[17];    d[23] = t.d[18];    d[29] = t.d[19];    d[35] = t.d[20];
}

inline tens4dmm::tens4dmm(double m[6][6])
{
    d[ 0]=m[0][0]; d[ 6]=m[0][1]; d[12]=m[0][2]; d[18]=m[0][3]; d[24]=m[0][4]; d[30]=m[0][5];
    d[ 1]=m[1][0]; d[ 7]=m[1][1]; d[13]=m[1][2]; d[19]=m[1][3]; d[25]=m[1][4]; d[31]=m[1][5];
    d[ 2]=m[2][0]; d[ 8]=m[2][1]; d[14]=m[2][2]; d[20]=m[2][3]; d[26]=m[2][4]; d[32]=m[2][5];
    d[ 3]=m[3][0]; d[ 9]=m[3][1]; d[15]=m[3][2]; d[21]=m[3][3]; d[27]=m[3][4]; d[33]=m[3][5];
    d[ 4]=m[4][0]; d[10]=m[4][1]; d[16]=m[4][2]; d[22]=m[4][3]; d[28]=m[4][4]; d[34]=m[4][5];
    d[ 5]=m[5][0]; d[11]=m[5][1]; d[17]=m[5][2]; d[23]=m[5][3]; d[29]=m[5][4]; d[35]=m[5][5];
}

inline double& tens4dmm::operator () (int i, int j, int k, int l)
{
    const int m[3][3] = {{0,3,5},{3,1,4},{5,4,2}};
    tens4dmm& T = (*this);
    return T(m[i][j], m[k][l]);
}

inline double tens4dmm::operator () (int i, int j, int k, int l) const
{
    const int m[3][3] = {{0,3,5},{3,1,4},{5,4,2}};
    const tens4dmm& T = (*this);
    return T(m[i][j], m[k][l]);
}

inline double& tens4dmm::operator () (int i, int j)
{
    const int m[6] = {0, 6, 12, 18, 24, 30};
    return d[m[j]+i];
}

inline double tens4dmm::operator () (int i, int j) const
{
    const int m[6] = {0, 6, 12, 18, 24, 30};
    return d[m[j]+i];
}

//-----------------------------------------------------------------------------
// operator +
inline tens4dmm tens4dmm::operator + (const tens4dmm& t) const
{
    tens4dmm s;
//    for (int i=0; i<NNZ; i++)
//        s.d[i] = d[i] + t.d[i];
    s.d[ 0] = d[ 0] + t.d[ 0];  s.d[12] = d[12] + t.d[12];  s.d[24] = d[24] + t.d[24];
    s.d[ 1] = d[ 1] + t.d[ 1];  s.d[13] = d[13] + t.d[13];  s.d[25] = d[25] + t.d[25];
    s.d[ 2] = d[ 2] + t.d[ 2];  s.d[14] = d[14] + t.d[14];  s.d[26] = d[26] + t.d[26];
    s.d[ 3] = d[ 3] + t.d[ 3];  s.d[15] = d[15] + t.d[15];  s.d[27] = d[27] + t.d[27];
    s.d[ 4] = d[ 4] + t.d[ 4];  s.d[16] = d[16] + t.d[16];  s.d[28] = d[28] + t.d[28];
    s.d[ 5] = d[ 5] + t.d[ 5];  s.d[17] = d[17] + t.d[17];  s.d[29] = d[29] + t.d[29];
    s.d[ 6] = d[ 6] + t.d[ 6];  s.d[18] = d[18] + t.d[18];  s.d[30] = d[30] + t.d[30];
    s.d[ 7] = d[ 7] + t.d[ 7];  s.d[19] = d[19] + t.d[19];  s.d[31] = d[31] + t.d[31];
    s.d[ 8] = d[ 8] + t.d[ 8];  s.d[20] = d[20] + t.d[20];  s.d[32] = d[32] + t.d[32];
    s.d[ 9] = d[ 9] + t.d[ 9];  s.d[21] = d[21] + t.d[21];  s.d[33] = d[33] + t.d[33];
    s.d[10] = d[10] + t.d[10];  s.d[22] = d[22] + t.d[22];  s.d[34] = d[34] + t.d[34];
    s.d[11] = d[11] + t.d[11];  s.d[23] = d[23] + t.d[23];  s.d[35] = d[35] + t.d[35];
    
    return s;
}

//-----------------------------------------------------------------------------
// operator -
inline tens4dmm tens4dmm::operator - (const tens4dmm& t) const
{
    tens4dmm s;
//    for (int i=0; i<NNZ; i++)
//        s.d[i] = d[i] - t.d[i];
    s.d[ 0] = d[ 0] - t.d[ 0];  s.d[12] = d[12] - t.d[12];  s.d[24] = d[24] - t.d[24];
    s.d[ 1] = d[ 1] - t.d[ 1];  s.d[13] = d[13] - t.d[13];  s.d[25] = d[25] - t.d[25];
    s.d[ 2] = d[ 2] - t.d[ 2];  s.d[14] = d[14] - t.d[14];  s.d[26] = d[26] - t.d[26];
    s.d[ 3] = d[ 3] - t.d[ 3];  s.d[15] = d[15] - t.d[15];  s.d[27] = d[27] - t.d[27];
    s.d[ 4] = d[ 4] - t.d[ 4];  s.d[16] = d[16] - t.d[16];  s.d[28] = d[28] - t.d[28];
    s.d[ 5] = d[ 5] - t.d[ 5];  s.d[17] = d[17] - t.d[17];  s.d[29] = d[29] - t.d[29];
    s.d[ 6] = d[ 6] - t.d[ 6];  s.d[18] = d[18] - t.d[18];  s.d[30] = d[30] - t.d[30];
    s.d[ 7] = d[ 7] - t.d[ 7];  s.d[19] = d[19] - t.d[19];  s.d[31] = d[31] - t.d[31];
    s.d[ 8] = d[ 8] - t.d[ 8];  s.d[20] = d[20] - t.d[20];  s.d[32] = d[32] - t.d[32];
    s.d[ 9] = d[ 9] - t.d[ 9];  s.d[21] = d[21] - t.d[21];  s.d[33] = d[33] - t.d[33];
    s.d[10] = d[10] - t.d[10];  s.d[22] = d[22] - t.d[22];  s.d[34] = d[34] - t.d[34];
    s.d[11] = d[11] - t.d[11];  s.d[23] = d[23] - t.d[23];  s.d[35] = d[35] - t.d[35];
    
    return s;
}

//-----------------------------------------------------------------------------
// operator + tens4ds
inline tens4dmm tens4dmm::operator + (const tens4ds& t) const
{
    tens4dmm s;
    
    s.d[ 0] = d[ 0] + t.d[ 0];    s.d[12] = d[12] + t.d[ 3];    s.d[24] = d[24] + t.d[10];
    s.d[ 1] = d[ 1] + t.d[ 1];    s.d[13] = d[13] + t.d[ 4];    s.d[25] = d[25] + t.d[11];
    s.d[ 2] = d[ 2] + t.d[ 3];    s.d[14] = d[14] + t.d[ 5];    s.d[26] = d[26] + t.d[12];
    s.d[ 3] = d[ 3] + t.d[ 6];    s.d[15] = d[15] + t.d[ 8];    s.d[27] = d[27] + t.d[13];
    s.d[ 4] = d[ 4] + t.d[10];    s.d[16] = d[16] + t.d[12];    s.d[28] = d[28] + t.d[14];
    s.d[ 5] = d[ 5] + t.d[15];    s.d[17] = d[17] + t.d[17];    s.d[29] = d[29] + t.d[19];
    s.d[ 6] = d[ 6] + t.d[ 1];    s.d[18] = d[18] + t.d[ 6];    s.d[30] = d[30] + t.d[15];
    s.d[ 7] = d[ 7] + t.d[ 2];    s.d[19] = d[19] + t.d[ 7];    s.d[31] = d[31] + t.d[16];
    s.d[ 8] = d[ 8] + t.d[ 4];    s.d[20] = d[20] + t.d[ 8];    s.d[32] = d[32] + t.d[17];
    s.d[ 9] = d[ 9] + t.d[ 7];    s.d[21] = d[21] + t.d[ 9];    s.d[33] = d[33] + t.d[18];
    s.d[10] = d[10] + t.d[11];    s.d[22] = d[22] + t.d[13];    s.d[34] = d[34] + t.d[19];
    s.d[11] = d[11] + t.d[16];    s.d[23] = d[23] + t.d[18];    s.d[35] = d[35] + t.d[20];
    
    return s;
}

//-----------------------------------------------------------------------------
// operator - tens4ds
inline tens4dmm tens4dmm::operator - (const tens4ds& t) const
{
    tens4dmm s;
    
    s.d[ 0] = d[ 0] - t.d[ 0];    s.d[12] = d[12] - t.d[ 3];    s.d[24] = d[24] - t.d[10];
    s.d[ 1] = d[ 1] - t.d[ 1];    s.d[13] = d[13] - t.d[ 4];    s.d[25] = d[25] - t.d[11];
    s.d[ 2] = d[ 2] - t.d[ 3];    s.d[14] = d[14] - t.d[ 5];    s.d[26] = d[26] - t.d[12];
    s.d[ 3] = d[ 3] - t.d[ 6];    s.d[15] = d[15] - t.d[ 8];    s.d[27] = d[27] - t.d[13];
    s.d[ 4] = d[ 4] - t.d[10];    s.d[16] = d[16] - t.d[12];    s.d[28] = d[28] - t.d[14];
    s.d[ 5] = d[ 5] - t.d[15];    s.d[17] = d[17] - t.d[17];    s.d[29] = d[29] - t.d[19];
    s.d[ 6] = d[ 6] - t.d[ 1];    s.d[18] = d[18] - t.d[ 6];    s.d[30] = d[30] - t.d[15];
    s.d[ 7] = d[ 7] - t.d[ 2];    s.d[19] = d[19] - t.d[ 7];    s.d[31] = d[31] - t.d[16];
    s.d[ 8] = d[ 8] - t.d[ 4];    s.d[20] = d[20] - t.d[ 8];    s.d[32] = d[32] - t.d[17];
    s.d[ 9] = d[ 9] - t.d[ 7];    s.d[21] = d[21] - t.d[ 9];    s.d[33] = d[33] - t.d[18];
    s.d[10] = d[10] - t.d[11];    s.d[22] = d[22] - t.d[13];    s.d[34] = d[34] - t.d[19];
    s.d[11] = d[11] - t.d[16];    s.d[23] = d[23] - t.d[18];    s.d[35] = d[35] - t.d[20];

    return s;
}

//-----------------------------------------------------------------------------
// operator *
inline tens4dmm tens4dmm::operator * (double g) const
{
    tens4dmm s;
//    for (int i=0; i<NNZ; i++)
//        s.d[i] = g*d[i];
    s.d[ 0] = g*d[ 0];    s.d[12] = g*d[12];    s.d[24] = g*d[24];
    s.d[ 1] = g*d[ 1];    s.d[13] = g*d[13];    s.d[25] = g*d[25];
    s.d[ 2] = g*d[ 2];    s.d[14] = g*d[14];    s.d[26] = g*d[26];
    s.d[ 3] = g*d[ 3];    s.d[15] = g*d[15];    s.d[27] = g*d[27];
    s.d[ 4] = g*d[ 4];    s.d[16] = g*d[16];    s.d[28] = g*d[28];
    s.d[ 5] = g*d[ 5];    s.d[17] = g*d[17];    s.d[29] = g*d[29];
    s.d[ 6] = g*d[ 6];    s.d[18] = g*d[18];    s.d[30] = g*d[30];
    s.d[ 7] = g*d[ 7];    s.d[19] = g*d[19];    s.d[31] = g*d[31];
    s.d[ 8] = g*d[ 8];    s.d[20] = g*d[20];    s.d[32] = g*d[32];
    s.d[ 9] = g*d[ 9];    s.d[21] = g*d[21];    s.d[33] = g*d[33];
    s.d[10] = g*d[10];    s.d[22] = g*d[22];    s.d[34] = g*d[34];
    s.d[11] = g*d[11];    s.d[23] = g*d[23];    s.d[35] = g*d[35];
    
    return s;
}

//-----------------------------------------------------------------------------
// operator /
inline tens4dmm tens4dmm::operator / (double g) const
{
    tens4dmm s;
//    for (int i=0; i<NNZ; i++)
//        s.d[i] = d[i]/g;
    s.d[ 0] = d[ 0]/g;    s.d[12] = d[12]/g;    s.d[24] = d[24]/g;
    s.d[ 1] = d[ 1]/g;    s.d[13] = d[13]/g;    s.d[25] = d[25]/g;
    s.d[ 2] = d[ 2]/g;    s.d[14] = d[14]/g;    s.d[26] = d[26]/g;
    s.d[ 3] = d[ 3]/g;    s.d[15] = d[15]/g;    s.d[27] = d[27]/g;
    s.d[ 4] = d[ 4]/g;    s.d[16] = d[16]/g;    s.d[28] = d[28]/g;
    s.d[ 5] = d[ 5]/g;    s.d[17] = d[17]/g;    s.d[29] = d[29]/g;
    s.d[ 6] = d[ 6]/g;    s.d[18] = d[18]/g;    s.d[30] = d[30]/g;
    s.d[ 7] = d[ 7]/g;    s.d[19] = d[19]/g;    s.d[31] = d[31]/g;
    s.d[ 8] = d[ 8]/g;    s.d[20] = d[20]/g;    s.d[32] = d[32]/g;
    s.d[ 9] = d[ 9]/g;    s.d[21] = d[21]/g;    s.d[33] = d[33]/g;
    s.d[10] = d[10]/g;    s.d[22] = d[22]/g;    s.d[34] = d[34]/g;
    s.d[11] = d[11]/g;    s.d[23] = d[23]/g;    s.d[35] = d[35]/g;
    
    return s;
}

//-----------------------------------------------------------------------------
// assignment operator +=
inline tens4dmm& tens4dmm::operator += (const tens4dmm& t)
{
//    for (int i=0; i<NNZ; i++)
//        d[i] += t.d[i];
    d[ 0] += t.d[ 0];    d[12] += t.d[12];    d[24] += t.d[24];
    d[ 1] += t.d[ 1];    d[13] += t.d[13];    d[25] += t.d[25];
    d[ 2] += t.d[ 2];    d[14] += t.d[14];    d[26] += t.d[26];
    d[ 3] += t.d[ 3];    d[15] += t.d[15];    d[27] += t.d[27];
    d[ 4] += t.d[ 4];    d[16] += t.d[16];    d[28] += t.d[28];
    d[ 5] += t.d[ 5];    d[17] += t.d[17];    d[29] += t.d[29];
    d[ 6] += t.d[ 6];    d[18] += t.d[18];    d[30] += t.d[30];
    d[ 7] += t.d[ 7];    d[19] += t.d[19];    d[31] += t.d[31];
    d[ 8] += t.d[ 8];    d[20] += t.d[20];    d[32] += t.d[32];
    d[ 9] += t.d[ 9];    d[21] += t.d[21];    d[33] += t.d[33];
    d[10] += t.d[10];    d[22] += t.d[22];    d[34] += t.d[34];
    d[11] += t.d[11];    d[23] += t.d[23];    d[35] += t.d[35];
    
    return (*this);
}

//-----------------------------------------------------------------------------
// assignment operator -=
inline tens4dmm& tens4dmm::operator -= (const tens4dmm& t)
{
//    for (int i=0; i<NNZ; i++)
//        d[i] -= t.d[i];
    d[ 0] -= t.d[ 0];    d[12] -= t.d[12];    d[24] -= t.d[24];
    d[ 1] -= t.d[ 1];    d[13] -= t.d[13];    d[25] -= t.d[25];
    d[ 2] -= t.d[ 2];    d[14] -= t.d[14];    d[26] -= t.d[26];
    d[ 3] -= t.d[ 3];    d[15] -= t.d[15];    d[27] -= t.d[27];
    d[ 4] -= t.d[ 4];    d[16] -= t.d[16];    d[28] -= t.d[28];
    d[ 5] -= t.d[ 5];    d[17] -= t.d[17];    d[29] -= t.d[29];
    d[ 6] -= t.d[ 6];    d[18] -= t.d[18];    d[30] -= t.d[30];
    d[ 7] -= t.d[ 7];    d[19] -= t.d[19];    d[31] -= t.d[31];
    d[ 8] -= t.d[ 8];    d[20] -= t.d[20];    d[32] -= t.d[32];
    d[ 9] -= t.d[ 9];    d[21] -= t.d[21];    d[33] -= t.d[33];
    d[10] -= t.d[10];    d[22] -= t.d[22];    d[34] -= t.d[34];
    d[11] -= t.d[11];    d[23] -= t.d[23];    d[35] -= t.d[35];
    
    return (*this);
}

//-----------------------------------------------------------------------------
// assignment operator += tens4ds
inline tens4dmm& tens4dmm::operator += (const tens4ds& t)
{
    d[ 0] += t.d[ 0];    d[12] += t.d[ 3];    d[24] += t.d[10];
    d[ 1] += t.d[ 1];    d[13] += t.d[ 4];    d[25] += t.d[11];
    d[ 2] += t.d[ 3];    d[14] += t.d[ 5];    d[26] += t.d[12];
    d[ 3] += t.d[ 6];    d[15] += t.d[ 8];    d[27] += t.d[13];
    d[ 4] += t.d[10];    d[16] += t.d[12];    d[28] += t.d[14];
    d[ 5] += t.d[15];    d[17] += t.d[17];    d[29] += t.d[19];
    d[ 6] += t.d[ 1];    d[18] += t.d[ 6];    d[30] += t.d[15];
    d[ 7] += t.d[ 2];    d[19] += t.d[ 7];    d[31] += t.d[16];
    d[ 8] += t.d[ 4];    d[20] += t.d[ 8];    d[32] += t.d[17];
    d[ 9] += t.d[ 7];    d[21] += t.d[ 9];    d[33] += t.d[18];
    d[10] += t.d[11];    d[22] += t.d[13];    d[34] += t.d[19];
    d[11] += t.d[16];    d[23] += t.d[18];    d[35] += t.d[20];
    
    return (*this);
}

//-----------------------------------------------------------------------------
// assignment operator -= tens4ds
inline tens4dmm& tens4dmm::operator -= (const tens4ds& t)
{
    d[ 0] -= t.d[ 0];    d[12] -= t.d[ 3];    d[24] -= t.d[10];
    d[ 1] -= t.d[ 1];    d[13] -= t.d[ 4];    d[25] -= t.d[11];
    d[ 2] -= t.d[ 3];    d[14] -= t.d[ 5];    d[26] -= t.d[12];
    d[ 3] -= t.d[ 6];    d[15] -= t.d[ 8];    d[27] -= t.d[13];
    d[ 4] -= t.d[10];    d[16] -= t.d[12];    d[28] -= t.d[14];
    d[ 5] -= t.d[15];    d[17] -= t.d[17];    d[29] -= t.d[19];
    d[ 6] -= t.d[ 1];    d[18] -= t.d[ 6];    d[30] -= t.d[15];
    d[ 7] -= t.d[ 2];    d[19] -= t.d[ 7];    d[31] -= t.d[16];
    d[ 8] -= t.d[ 4];    d[20] -= t.d[ 8];    d[32] -= t.d[17];
    d[ 9] -= t.d[ 7];    d[21] -= t.d[ 9];    d[33] -= t.d[18];
    d[10] -= t.d[11];    d[22] -= t.d[13];    d[34] -= t.d[19];
    d[11] -= t.d[16];    d[23] -= t.d[18];    d[35] -= t.d[20];
    
    return (*this);
}

//-----------------------------------------------------------------------------
// assignment operator *=
inline tens4dmm& tens4dmm::operator *= (double g)
{
//    for (int i=0; i<NNZ; i++)
//        d[i] *= g;
    d[ 0] *= g;    d[12] *= g;    d[24] *= g;
    d[ 1] *= g;    d[13] *= g;    d[25] *= g;
    d[ 2] *= g;    d[14] *= g;    d[26] *= g;
    d[ 3] *= g;    d[15] *= g;    d[27] *= g;
    d[ 4] *= g;    d[16] *= g;    d[28] *= g;
    d[ 5] *= g;    d[17] *= g;    d[29] *= g;
    d[ 6] *= g;    d[18] *= g;    d[30] *= g;
    d[ 7] *= g;    d[19] *= g;    d[31] *= g;
    d[ 8] *= g;    d[20] *= g;    d[32] *= g;
    d[ 9] *= g;    d[21] *= g;    d[33] *= g;
    d[10] *= g;    d[22] *= g;    d[34] *= g;
    d[11] *= g;    d[23] *= g;    d[35] *= g;
    
    return (*this);
}

//-----------------------------------------------------------------------------
// assignment operator /=
inline tens4dmm& tens4dmm::operator /= (double g)
{
//    for (int i=0; i<NNZ; i++)
//        d[i] /= g;
    d[ 0] /= g;    d[12] /= g;    d[24] /= g;
    d[ 1] /= g;    d[13] /= g;    d[25] /= g;
    d[ 2] /= g;    d[14] /= g;    d[26] /= g;
    d[ 3] /= g;    d[15] /= g;    d[27] /= g;
    d[ 4] /= g;    d[16] /= g;    d[28] /= g;
    d[ 5] /= g;    d[17] /= g;    d[29] /= g;
    d[ 6] /= g;    d[18] /= g;    d[30] /= g;
    d[ 7] /= g;    d[19] /= g;    d[31] /= g;
    d[ 8] /= g;    d[20] /= g;    d[32] /= g;
    d[ 9] /= g;    d[21] /= g;    d[33] /= g;
    d[10] /= g;    d[22] /= g;    d[34] /= g;
    d[11] /= g;    d[23] /= g;    d[35] /= g;
    
    return (*this);
}

//-----------------------------------------------------------------------------
// unary operator -
inline tens4dmm tens4dmm::operator - () const
{
    tens4dmm s;
    s.d[ 0] = -d[ 0];    s.d[12] = -d[12];    s.d[24] = -d[24];
    s.d[ 1] = -d[ 1];    s.d[13] = -d[13];    s.d[25] = -d[25];
    s.d[ 2] = -d[ 2];    s.d[14] = -d[14];    s.d[26] = -d[26];
    s.d[ 3] = -d[ 3];    s.d[15] = -d[15];    s.d[27] = -d[27];
    s.d[ 4] = -d[ 4];    s.d[16] = -d[16];    s.d[28] = -d[28];
    s.d[ 5] = -d[ 5];    s.d[17] = -d[17];    s.d[29] = -d[29];
    s.d[ 6] = -d[ 6];    s.d[18] = -d[18];    s.d[30] = -d[30];
    s.d[ 7] = -d[ 7];    s.d[19] = -d[19];    s.d[31] = -d[31];
    s.d[ 8] = -d[ 8];    s.d[20] = -d[20];    s.d[32] = -d[32];
    s.d[ 9] = -d[ 9];    s.d[21] = -d[21];    s.d[33] = -d[33];
    s.d[10] = -d[10];    s.d[22] = -d[22];    s.d[34] = -d[34];
    s.d[11] = -d[11];    s.d[23] = -d[23];    s.d[35] = -d[35];
    
    return s;
}

//-----------------------------------------------------------------------------
// vdotTdotv_jk = a_i T_ijkl b_l
inline mat3d vdotTdotv(const vec3d& a, const tens4dmm& T, const vec3d& b)
{
    return mat3d(a.x*b.x*T.d[0] + a.y*b.x*T.d[3] + a.z*b.x*T.d[5] + a.x*b.y*T.d[18] + a.y*b.y*T.d[21] +
                 a.z*b.y*T.d[23] + a.x*b.z*T.d[30] + a.y*b.z*T.d[33] + a.z*b.z*T.d[35],
                 a.x*b.y*T.d[6] + a.y*b.y*T.d[9] + a.z*b.y*T.d[11] + a.x*b.x*T.d[18] + a.y*b.x*T.d[21] + a.z*b.x*T.d[23] +
                 a.x*b.z*T.d[24] + a.y*b.z*T.d[27] + a.z*b.z*T.d[29],
                 a.x*b.z*T.d[12] + a.y*b.z*T.d[15] + a.z*b.z*T.d[17] + a.x*b.y*T.d[24] + a.y*b.y*T.d[27] + a.z*b.y*T.d[29] +
                 a.x*b.x*T.d[30] + a.y*b.x*T.d[33] + a.z*b.x*T.d[35],
                 a.y*b.x*T.d[1] + a.x*b.x*T.d[3] + a.z*b.x*T.d[4] + a.y*b.y*T.d[19] + a.x*b.y*T.d[21] + a.z*b.y*T.d[22] +
                 a.y*b.z*T.d[31] + a.x*b.z*T.d[33] + a.z*b.z*T.d[34],
                 a.y*b.y*T.d[7] + a.x*b.y*T.d[9] + a.z*b.y*T.d[10] + a.y*b.x*T.d[19] + a.x*b.x*T.d[21] + a.z*b.x*T.d[22] +
                 a.y*b.z*T.d[25] + a.x*b.z*T.d[27] + a.z*b.z*T.d[28],
                 a.y*b.z*T.d[13] + a.x*b.z*T.d[15] + a.z*b.z*T.d[16] + a.y*b.y*T.d[25] + a.x*b.y*T.d[27] + a.z*b.y*T.d[28] +
                 a.y*b.x*T.d[31] + a.x*b.x*T.d[33] + a.z*b.x*T.d[34],
                 a.z*b.x*T.d[2] + a.y*b.x*T.d[4] + a.x*b.x*T.d[5] + a.z*b.y*T.d[20] + a.y*b.y*T.d[22] + a.x*b.y*T.d[23] +
                 a.z*b.z*T.d[32] + a.y*b.z*T.d[34] + a.x*b.z*T.d[35],
                 a.z*b.y*T.d[8] + a.y*b.y*T.d[10] + a.x*b.y*T.d[11] + a.z*b.x*T.d[20] + a.y*b.x*T.d[22] + a.x*b.x*T.d[23] +
                 a.z*b.z*T.d[26] + a.y*b.z*T.d[28] + a.x*b.z*T.d[29],
                 a.z*b.z*T.d[14] + a.y*b.z*T.d[16] + a.x*b.z*T.d[17] + a.z*b.y*T.d[26] + a.y*b.y*T.d[28] + a.x*b.y*T.d[29] +
                 a.z*b.x*T.d[32] + a.y*b.x*T.d[34] + a.x*b.x*T.d[35]);
}

//-----------------------------------------------------------------------------
// double contraction of general 4th-order tensor with a general 2nd-order tensor
// Aij = Dijkl Mkl
inline mat3ds tens4dmm::dot(const mat3ds &m) const
{
    mat3ds a;
    
    a.xx() = d[0]*m.xx() + d[18]*m.xy() + d[30]*m.xz() + d[6]*m.yy() + d[24]*m.yz() + d[12]*m.zz();
    a.yy() = d[1]*m.xx() + d[19]*m.xy() + d[31]*m.xz() + d[7]*m.yy() + d[25]*m.yz() + d[13]*m.zz();
    a.zz() = d[2]*m.xx() + d[20]*m.xy() + d[32]*m.xz() + d[8]*m.yy() + d[26]*m.yz() + d[14]*m.zz();
    a.xy() = d[3]*m.xx() + d[21]*m.xy() + d[33]*m.xz() + d[9]*m.yy() + d[27]*m.yz() + d[15]*m.zz();
    a.yz() = d[4]*m.xx() + d[22]*m.xy() + d[34]*m.xz() + d[10]*m.yy() + d[28]*m.yz() + d[16]*m.zz();
    a.xz() = d[5]*m.xx() + d[23]*m.xy() + d[35]*m.xz() + d[11]*m.yy() + d[29]*m.yz() + d[17]*m.zz();
    
    return a;
}

//-----------------------------------------------------------------------------
// trace
// C.tr() = I:C:I
inline double tens4dmm::tr() const
{
    return (d[0]+d[1]+d[2]+d[6]+d[7]+d[8]+d[12]+d[13]+d[14]);
}

//-----------------------------------------------------------------------------
// intialize to zero
inline void tens4dmm::zero()
{
    d[ 0] = d[ 1] = d[ 2] = d[ 3] = d[ 4] = d[ 5] = d[ 6] = d[ 7] = d[ 8] =
    d[ 9] = d[10] = d[11] = d[12] = d[13] = d[14] = d[15] = d[16] = d[17] =
    d[18] = d[19] = d[20] = d[21] = d[22] = d[23] = d[24] = d[25] = d[26] =
    d[27] = d[28] = d[29] = d[30] = d[31] = d[32] = d[33] = d[34] = d[35] = 0;
}

//-----------------------------------------------------------------------------
// extract 6x6 matrix
inline void tens4dmm::extract(double D[6][6])
{
    D[0][0] = d[0]; D[0][1] = d[ 6]; D[0][2] = d[12]; D[0][3] = d[18]; D[0][4] = d[24]; D[0][5] = d[30];
    D[1][0] = d[1]; D[1][1] = d[ 7]; D[1][2] = d[13]; D[1][3] = d[19]; D[1][4] = d[25]; D[1][5] = d[31];
    D[2][0] = d[2]; D[2][1] = d[ 8]; D[2][2] = d[14]; D[2][3] = d[20]; D[2][4] = d[26]; D[2][5] = d[32];
    D[3][0] = d[3]; D[3][1] = d[ 9]; D[3][2] = d[15]; D[3][3] = d[21]; D[3][4] = d[27]; D[3][5] = d[33];
    D[4][0] = d[4]; D[4][1] = d[10]; D[4][2] = d[16]; D[4][3] = d[22]; D[4][4] = d[28]; D[4][5] = d[34];
    D[5][0] = d[5]; D[5][1] = d[11]; D[5][2] = d[17]; D[5][3] = d[23]; D[5][4] = d[29]; D[5][5] = d[35];
}


//-----------------------------------------------------------------------------
// compute the super symmetric (major and minor symmetric) component of the tensor
// Sijkl = (1/8)*(Cijkl + Cjikl + Cijlk + Cjilk + Cklij + Clkij + Cklji + Clkji)
inline tens4ds tens4dmm::supersymm() const
{
    tens4ds s;
    
    s.d[ 0] = d[0];
    s.d[ 1] = (d[1]+d[6])/2;
    s.d[ 2] = d[7];
    s.d[ 3] = (d[2]+d[12])/2;
    s.d[ 4] = (d[8]+d[13])/2;
    s.d[ 5] = d[14];
    s.d[ 6] = (d[3]+d[18])/2;
    s.d[ 7] = (d[9]+d[19])/2;
    s.d[ 8] = (d[15]+d[20])/2;
    s.d[ 9] = d[21];
    s.d[10] = (d[4]+d[24])/2;
    s.d[11] = (d[10]+d[25])/2;
    s.d[12] = (d[16]+d[26])/2;
    s.d[13] = (d[22]+d[27])/2;
    s.d[14] = d[28];
    s.d[15] = (d[5]+d[30])/2;
    s.d[16] = (d[11]+d[31])/2;
    s.d[17] = (d[17]+d[32])/2;
    s.d[18] = (d[23]+d[33])/2;
    s.d[19] = (d[29]+d[34])/2;
    s.d[20] = d[35];
    
    return s;
}

//-----------------------------------------------------------------------------
// compute the major transpose of the tensor
// Sijkl -> Sklij
inline tens4dmm tens4dmm::transpose() const
{
    tens4dmm s;
    
    s.d[ 0] = d[ 0];    s.d[ 2] = d[12];    s.d[ 4] = d[24];
    s.d[ 6] = d[ 1];    s.d[ 8] = d[13];    s.d[10] = d[25];
    s.d[12] = d[ 2];    s.d[14] = d[14];    s.d[16] = d[26];
    s.d[18] = d[ 3];    s.d[20] = d[15];    s.d[22] = d[27];
    s.d[24] = d[ 4];    s.d[26] = d[16];    s.d[28] = d[28];
    s.d[30] = d[ 5];    s.d[32] = d[17];    s.d[34] = d[29];
    s.d[ 1] = d[ 6];    s.d[ 3] = d[18];    s.d[ 5] = d[30];
    s.d[ 7] = d[ 7];    s.d[ 9] = d[19];    s.d[11] = d[31];
    s.d[13] = d[ 8];    s.d[15] = d[20];    s.d[17] = d[32];
    s.d[19] = d[ 9];    s.d[21] = d[21];    s.d[23] = d[33];
    s.d[25] = d[10];    s.d[27] = d[22];    s.d[29] = d[34];
    s.d[31] = d[11];    s.d[33] = d[23];    s.d[35] = d[35];
    return s;
}

//-----------------------------------------------------------------------------
// (a dyad1 b)_ijkl = a_ij b_kl
inline tens4dmm dyad1mm(const mat3ds& a, const mat3ds& b)
{
    tens4dmm c;
    
    c.d[ 0] = a.xx()*b.xx();    c.d[12] = a.xx()*b.zz();    c.d[24] = a.xx()*b.yz();
    c.d[ 1] = a.yy()*b.xx();    c.d[13] = a.yy()*b.zz();    c.d[25] = a.yy()*b.yz();
    c.d[ 2] = a.zz()*b.xx();    c.d[14] = a.zz()*b.zz();    c.d[26] = a.zz()*b.yz();
    c.d[ 3] = a.xy()*b.xx();    c.d[15] = a.xy()*b.zz();    c.d[27] = a.xy()*b.yz();
    c.d[ 4] = a.yz()*b.xx();    c.d[16] = a.yz()*b.zz();    c.d[28] = a.yz()*b.yz();
    c.d[ 5] = a.xz()*b.xx();    c.d[17] = a.xz()*b.zz();    c.d[29] = a.xz()*b.yz();
    c.d[ 6] = a.xx()*b.yy();    c.d[18] = a.xx()*b.xy();    c.d[30] = a.xx()*b.xz();
    c.d[ 7] = a.yy()*b.yy();    c.d[19] = a.yy()*b.xy();    c.d[31] = a.yy()*b.xz();
    c.d[ 8] = a.zz()*b.yy();    c.d[20] = a.zz()*b.xy();    c.d[32] = a.zz()*b.xz();
    c.d[ 9] = a.xy()*b.yy();    c.d[21] = a.xy()*b.xy();    c.d[33] = a.xy()*b.xz();
    c.d[10] = a.yz()*b.yy();    c.d[22] = a.yz()*b.xy();    c.d[34] = a.yz()*b.xz();
    c.d[11] = a.xz()*b.yy();    c.d[23] = a.xz()*b.xy();    c.d[35] = a.xz()*b.xz();
    
    return c;
}

//-----------------------------------------------------------------------------
// (a dyad2 b)_ijkl = a_ik b_jl
inline tens4dmm dyad2mm(const mat3ds& a, const mat3ds& b)
{
    tens4dmm c;
    
    c.d[ 0] = a.xx()*b.xx();    c.d[12] = a.xz()*b.xz();    c.d[24] = a.xy()*b.xz();
    c.d[ 1] = a.xy()*b.xy();    c.d[13] = a.yz()*b.yz();    c.d[25] = a.yy()*b.yz();
    c.d[ 2] = a.xz()*b.xz();    c.d[14] = a.zz()*b.zz();    c.d[26] = a.yz()*b.zz();
    c.d[ 3] = a.xx()*b.xy();    c.d[15] = a.xz()*b.yz();    c.d[27] = a.xy()*b.yz();
    c.d[ 4] = a.xy()*b.xz();    c.d[16] = a.yz()*b.zz();    c.d[28] = a.yy()*b.zz();
    c.d[ 5] = a.xx()*b.xz();    c.d[17] = a.xz()*b.zz();    c.d[29] = a.xy()*b.zz();
    c.d[ 6] = a.xy()*b.xy();    c.d[18] = a.xx()*b.xy();    c.d[30] = a.xx()*b.xz();
    c.d[ 7] = a.yy()*b.yy();    c.d[19] = a.xy()*b.yy();    c.d[31] = a.xy()*b.yz();
    c.d[ 8] = a.yz()*b.yz();    c.d[20] = a.xz()*b.yz();    c.d[32] = a.xz()*b.zz();
    c.d[ 9] = a.xy()*b.yy();    c.d[21] = a.xx()*b.yy();    c.d[33] = a.xx()*b.yz();
    c.d[10] = a.yy()*b.yz();    c.d[22] = a.xy()*b.yz();    c.d[34] = a.xy()*b.zz();
    c.d[11] = a.xy()*b.yz();    c.d[23] = a.xx()*b.yz();    c.d[35] = a.xx()*b.zz();
    
    return c;
}

//-----------------------------------------------------------------------------
// (a dyad3 b)_ijkl = a_il b_jk
inline tens4dmm dyad3mm(const mat3ds& a, const mat3ds& b)
{
    tens4dmm c;

    c.d[ 0] = a.xx()*b.xx();    c.d[12] = a.xz()*b.xz();    c.d[24] = a.xz()*b.xy();
    c.d[ 1] = a.xy()*b.xy();    c.d[13] = a.yz()*b.yz();    c.d[25] = a.yz()*b.yy();
    c.d[ 2] = a.xz()*b.xz();    c.d[14] = a.zz()*b.zz();    c.d[26] = a.zz()*b.yz();
    c.d[ 3] = a.xx()*b.xy();    c.d[15] = a.xz()*b.yz();    c.d[27] = a.xz()*b.yy();
    c.d[ 4] = a.xy()*b.xz();    c.d[16] = a.yz()*b.zz();    c.d[28] = a.yz()*b.yz();
    c.d[ 5] = a.xx()*b.xz();    c.d[17] = a.xz()*b.zz();    c.d[29] = a.xz()*b.yz();
    c.d[ 6] = a.xy()*b.xy();    c.d[18] = a.xy()*b.xx();    c.d[30] = a.xz()*b.xx();
    c.d[ 7] = a.yy()*b.yy();    c.d[19] = a.yy()*b.xy();    c.d[31] = a.yz()*b.xy();
    c.d[ 8] = a.yz()*b.yz();    c.d[20] = a.yz()*b.xz();    c.d[32] = a.zz()*b.xz();
    c.d[ 9] = a.xy()*b.yy();    c.d[21] = a.xy()*b.xy();    c.d[33] = a.xz()*b.xy();
    c.d[10] = a.yy()*b.yz();    c.d[22] = a.yy()*b.xz();    c.d[34] = a.yz()*b.xz();
    c.d[11] = a.xy()*b.yz();    c.d[23] = a.xy()*b.xz();    c.d[35] = a.xz()*b.xz();

    return c;
}

//-----------------------------------------------------------------------------
// (a dyad4 b)_ijkl = 0.5(a_ik b_jl + a_il b_jk)
inline tens4dmm dyad4mm(const mat3ds& a, const mat3ds& b)
{
    return (dyad2mm(a,b) + dyad3mm(a,b))/2;
}

//-----------------------------------------------------------------------------
// (a dyad4 b)_ijkl = 0.5(a_ik b_jl + a_il b_jk)
inline tens4dmm dyad4mm(const mat3d& a, const mat3d& b)
{
    tens4dmm c;
    c.d[ 0] = 2*a(0,0)*b(0,0);
    c.d[ 1] = 2*a(1,0)*b(1,0);
    c.d[ 2] = 2*a(2,0)*b(2,0);
    c.d[ 3] = a(1,0)*b(0,0) + a(0,0)*b(1,0);
    c.d[ 4] = a(2,0)*b(1,0) + a(1,0)*b(2,0);
    c.d[ 5] = a(2,0)*b(0,0) + a(0,0)*b(2,0);

    c.d[ 6] = 2*a(0,1)*b(0,1);
    c.d[ 7] = 2*a(1,1)*b(1,1);
    c.d[ 8] = 2*a(2,1)*b(2,1);
    c.d[ 9] = a(1,1)*b(0,1) + a(0,1)*b(1,1);
    c.d[10] = a(2,1)*b(1,1) + a(1,1)*b(2,1);
    c.d[11] = a(2,1)*b(0,1) + a(0,1)*b(2,1);

    c.d[12] = 2*a(0,2)*b(0,2);
    c.d[13] = 2*a(1,2)*b(1,2);
    c.d[14] = 2*a(2,2)*b(2,2);
    c.d[15] = a(1,2)*b(0,2) + a(0,2)*b(1,2);
    c.d[16] = a(2,2)*b(1,2) + a(1,2)*b(2,2);
    c.d[17] = a(2,2)*b(0,2) + a(0,2)*b(2,2);

    c.d[18] = a(0,1)*b(0,0) + a(0,0)*b(0,1);
    c.d[19] = a(1,1)*b(1,0) + a(1,0)*b(1,1);
    c.d[20] = a(2,1)*b(2,0) + a(2,0)*b(2,1);
    c.d[21] = (a(1,1)*b(0,0) + a(1,0)*b(0,1) + a(0,1)*b(1,0) + a(0,0)*b(1,1))/2.;
    c.d[22] = (a(2,1)*b(1,0) + a(2,0)*b(1,1) + a(1,1)*b(2,0) + a(1,0)*b(2,1))/2.;
    c.d[23] = (a(2,1)*b(0,0) + a(2,0)*b(0,1) + a(0,1)*b(2,0) + a(0,0)*b(2,1))/2.;

    c.d[24] = a(0,2)*b(0,1) + a(0,1)*b(0,2);
    c.d[25] = a(1,2)*b(1,1) + a(1,1)*b(1,2);
    c.d[26] = a(2,2)*b(2,1) + a(2,1)*b(2,2);
    c.d[27] = (a(1,2)*b(0,1) + a(1,1)*b(0,2) + a(0,2)*b(1,1) + a(0,1)*b(1,2))/2.;
    c.d[28] = (a(2,2)*b(1,1) + a(2,1)*b(1,2) + a(1,2)*b(2,1) + a(1,1)*b(2,2))/2.;
    c.d[29] = (a(2,2)*b(0,1) + a(2,1)*b(0,2) + a(0,2)*b(2,1) + a(0,1)*b(2,2))/2.;

    c.d[30] = a(0,2)*b(0,0) + a(0,0)*b(0,2);
    c.d[31] = a(1,2)*b(1,0) + a(1,0)*b(1,2);
    c.d[32] = a(2,2)*b(2,0) + a(2,0)*b(2,2);
    c.d[33] = (a(1,2)*b(0,0) + a(1,0)*b(0,2) + a(0,2)*b(1,0) + a(0,0)*b(1,2))/2.;
    c.d[34] = (a(2,2)*b(1,0) + a(2,0)*b(1,2) + a(1,2)*b(2,0) + a(1,0)*b(2,2))/2.;
    c.d[35] = (a(2,2)*b(0,0) + a(2,0)*b(0,2) + a(0,2)*b(2,0) + a(0,0)*b(2,2))/2.;
    
    return c;
}

//-----------------------------------------------------------------------------
// (a ddot b)_ijkl = a_ijmn b_mnkl
inline tens4dmm ddot(const tens4dmm& a, const tens4dmm& b)
{
    tens4dmm c;
    
    c.d[ 0] = a.d[0]*b.d[0] + a.d[6]*b.d[1] + a.d[12]*b.d[2] + a.d[18]*b.d[3] + a.d[24]*b.d[4] + a.d[30]*b.d[5];
    c.d[ 1] = a.d[1]*b.d[0] + a.d[7]*b.d[1] + a.d[13]*b.d[2] + a.d[19]*b.d[3] + a.d[25]*b.d[4] + a.d[31]*b.d[5];
    c.d[ 2] = a.d[2]*b.d[0] + a.d[8]*b.d[1] + a.d[14]*b.d[2] + a.d[20]*b.d[3] + a.d[26]*b.d[4] + a.d[32]*b.d[5];
    c.d[ 3] = a.d[3]*b.d[0] + a.d[9]*b.d[1] + a.d[15]*b.d[2] + a.d[21]*b.d[3] + a.d[27]*b.d[4] + a.d[33]*b.d[5];
    c.d[ 4] = a.d[4]*b.d[0] + a.d[10]*b.d[1] + a.d[16]*b.d[2] + a.d[22]*b.d[3] + a.d[28]*b.d[4] + a.d[34]*b.d[5];
    c.d[ 5] = a.d[5]*b.d[0] + a.d[11]*b.d[1] + a.d[17]*b.d[2] + a.d[23]*b.d[3] + a.d[29]*b.d[4] + a.d[35]*b.d[5];
    c.d[ 6] = a.d[24]*b.d[10] + a.d[30]*b.d[11] + a.d[0]*b.d[6] + a.d[6]*b.d[7] + a.d[12]*b.d[8] + a.d[18]*b.d[9];
    c.d[ 7] = a.d[25]*b.d[10] + a.d[31]*b.d[11] + a.d[1]*b.d[6] + a.d[7]*b.d[7] + a.d[13]*b.d[8] + a.d[19]*b.d[9];
    c.d[ 8] = a.d[26]*b.d[10] + a.d[32]*b.d[11] + a.d[2]*b.d[6] + a.d[8]*b.d[7] + a.d[14]*b.d[8] + a.d[20]*b.d[9];
    c.d[ 9] = a.d[27]*b.d[10] + a.d[33]*b.d[11] + a.d[3]*b.d[6] + a.d[9]*b.d[7] + a.d[15]*b.d[8] + a.d[21]*b.d[9];
    c.d[10] = a.d[28]*b.d[10] + a.d[34]*b.d[11] + a.d[4]*b.d[6] + a.d[10]*b.d[7] + a.d[16]*b.d[8] + a.d[22]*b.d[9];
    c.d[11] = a.d[29]*b.d[10] + a.d[35]*b.d[11] + a.d[5]*b.d[6] + a.d[11]*b.d[7] + a.d[17]*b.d[8] + a.d[23]*b.d[9];
    c.d[12] = a.d[0]*b.d[12] + a.d[6]*b.d[13] + a.d[12]*b.d[14] + a.d[18]*b.d[15] + a.d[24]*b.d[16] + a.d[30]*b.d[17];
    c.d[13] = a.d[1]*b.d[12] + a.d[7]*b.d[13] + a.d[13]*b.d[14] + a.d[19]*b.d[15] + a.d[25]*b.d[16] + a.d[31]*b.d[17];
    c.d[14] = a.d[2]*b.d[12] + a.d[8]*b.d[13] + a.d[14]*b.d[14] + a.d[20]*b.d[15] + a.d[26]*b.d[16] + a.d[32]*b.d[17];
    c.d[15] = a.d[3]*b.d[12] + a.d[9]*b.d[13] + a.d[15]*b.d[14] + a.d[21]*b.d[15] + a.d[27]*b.d[16] + a.d[33]*b.d[17];
    c.d[16] = a.d[4]*b.d[12] + a.d[10]*b.d[13] + a.d[16]*b.d[14] + a.d[22]*b.d[15] + a.d[28]*b.d[16] + a.d[34]*b.d[17];
    c.d[17] = a.d[5]*b.d[12] + a.d[11]*b.d[13] + a.d[17]*b.d[14] + a.d[23]*b.d[15] + a.d[29]*b.d[16] + a.d[35]*b.d[17];
    c.d[18] = a.d[0]*b.d[18] + a.d[6]*b.d[19] + a.d[12]*b.d[20] + a.d[18]*b.d[21] + a.d[24]*b.d[22] + a.d[30]*b.d[23];
    c.d[19] = a.d[1]*b.d[18] + a.d[7]*b.d[19] + a.d[13]*b.d[20] + a.d[19]*b.d[21] + a.d[25]*b.d[22] + a.d[31]*b.d[23];
    c.d[20] = a.d[2]*b.d[18] + a.d[8]*b.d[19] + a.d[14]*b.d[20] + a.d[20]*b.d[21] + a.d[26]*b.d[22] + a.d[32]*b.d[23];
    c.d[21] = a.d[3]*b.d[18] + a.d[9]*b.d[19] + a.d[15]*b.d[20] + a.d[21]*b.d[21] + a.d[27]*b.d[22] + a.d[33]*b.d[23];
    c.d[22] = a.d[4]*b.d[18] + a.d[10]*b.d[19] + a.d[16]*b.d[20] + a.d[22]*b.d[21] + a.d[28]*b.d[22] + a.d[34]*b.d[23];
    c.d[23] = a.d[5]*b.d[18] + a.d[11]*b.d[19] + a.d[17]*b.d[20] + a.d[23]*b.d[21] + a.d[29]*b.d[22] + a.d[35]*b.d[23];
    c.d[24] = a.d[0]*b.d[24] + a.d[6]*b.d[25] + a.d[12]*b.d[26] + a.d[18]*b.d[27] + a.d[24]*b.d[28] + a.d[30]*b.d[29];
    c.d[25] = a.d[1]*b.d[24] + a.d[7]*b.d[25] + a.d[13]*b.d[26] + a.d[19]*b.d[27] + a.d[25]*b.d[28] + a.d[31]*b.d[29];
    c.d[26] = a.d[2]*b.d[24] + a.d[8]*b.d[25] + a.d[14]*b.d[26] + a.d[20]*b.d[27] + a.d[26]*b.d[28] + a.d[32]*b.d[29];
    c.d[27] = a.d[3]*b.d[24] + a.d[9]*b.d[25] + a.d[15]*b.d[26] + a.d[21]*b.d[27] + a.d[27]*b.d[28] + a.d[33]*b.d[29];
    c.d[28] = a.d[4]*b.d[24] + a.d[10]*b.d[25] + a.d[16]*b.d[26] + a.d[22]*b.d[27] + a.d[28]*b.d[28] + a.d[34]*b.d[29];
    c.d[29] = a.d[5]*b.d[24] + a.d[11]*b.d[25] + a.d[17]*b.d[26] + a.d[23]*b.d[27] + a.d[29]*b.d[28] + a.d[35]*b.d[29];
    c.d[30] = a.d[0]*b.d[30] + a.d[6]*b.d[31] + a.d[12]*b.d[32] + a.d[18]*b.d[33] + a.d[24]*b.d[34] + a.d[30]*b.d[35];
    c.d[31] = a.d[1]*b.d[30] + a.d[7]*b.d[31] + a.d[13]*b.d[32] + a.d[19]*b.d[33] + a.d[25]*b.d[34] + a.d[31]*b.d[35];
    c.d[32] = a.d[2]*b.d[30] + a.d[8]*b.d[31] + a.d[14]*b.d[32] + a.d[20]*b.d[33] + a.d[26]*b.d[34] + a.d[32]*b.d[35];
    c.d[33] = a.d[3]*b.d[30] + a.d[9]*b.d[31] + a.d[15]*b.d[32] + a.d[21]*b.d[33] + a.d[27]*b.d[34] + a.d[33]*b.d[35];
    c.d[34] = a.d[4]*b.d[30] + a.d[10]*b.d[31] + a.d[16]*b.d[32] + a.d[22]*b.d[33] + a.d[28]*b.d[34] + a.d[34]*b.d[35];
    c.d[35] = a.d[5]*b.d[30] + a.d[11]*b.d[31] + a.d[17]*b.d[32] + a.d[23]*b.d[33] + a.d[29]*b.d[34] + a.d[35]*b.d[35];
    
    return c;
    
}

//-----------------------------------------------------------------------------
// (a ddot b)_ijkl = a_ijmn b_mnkl where b is super-symmetric
inline tens4dmm ddot(const tens4dmm& a, const tens4ds& b)
{
    tens4dmm c;
    
    c.d[ 0] = a.d[0]*b.d[0] + a.d[6]*b.d[1] + a.d[12]*b.d[3] + a.d[18]*b.d[6] + a.d[24]*b.d[10] + a.d[30]*b.d[15];
    c.d[ 1] = a.d[1]*b.d[0] + a.d[7]*b.d[1] + a.d[13]*b.d[3] + a.d[19]*b.d[6] + a.d[25]*b.d[10] + a.d[31]*b.d[15];
    c.d[ 2] = a.d[2]*b.d[0] + a.d[8]*b.d[1] + a.d[14]*b.d[3] + a.d[20]*b.d[6] + a.d[26]*b.d[10] + a.d[32]*b.d[15];
    c.d[ 3] = a.d[3]*b.d[0] + a.d[9]*b.d[1] + a.d[15]*b.d[3] + a.d[21]*b.d[6] + a.d[27]*b.d[10] + a.d[33]*b.d[15];
    c.d[ 4] = a.d[4]*b.d[0] + a.d[10]*b.d[1] + a.d[16]*b.d[3] + a.d[22]*b.d[6] + a.d[28]*b.d[10] + a.d[34]*b.d[15];
    c.d[ 5] = a.d[5]*b.d[0] + a.d[11]*b.d[1] + a.d[17]*b.d[3] + a.d[23]*b.d[6] + a.d[29]*b.d[10] + a.d[35]*b.d[15];
    c.d[ 6] = a.d[0]*b.d[1] + a.d[6]*b.d[2] + a.d[12]*b.d[4] + a.d[18]*b.d[7] + a.d[24]*b.d[11] + a.d[30]*b.d[16];
    c.d[ 7] = a.d[1]*b.d[1] + a.d[7]*b.d[2] + a.d[13]*b.d[4] + a.d[19]*b.d[7] + a.d[25]*b.d[11] + a.d[31]*b.d[16];
    c.d[ 8] = a.d[2]*b.d[1] + a.d[8]*b.d[2] + a.d[14]*b.d[4] + a.d[20]*b.d[7] + a.d[26]*b.d[11] + a.d[32]*b.d[16];
    c.d[ 9] = a.d[3]*b.d[1] + a.d[9]*b.d[2] + a.d[15]*b.d[4] + a.d[21]*b.d[7] + a.d[27]*b.d[11] + a.d[33]*b.d[16];
    c.d[10] = a.d[4]*b.d[1] + a.d[10]*b.d[2] + a.d[16]*b.d[4] + a.d[22]*b.d[7] + a.d[28]*b.d[11] + a.d[34]*b.d[16];
    c.d[11] = a.d[5]*b.d[1] + a.d[11]*b.d[2] + a.d[17]*b.d[4] + a.d[23]*b.d[7] + a.d[29]*b.d[11] + a.d[35]*b.d[16];
    c.d[12] = a.d[0]*b.d[3] + a.d[6]*b.d[4] + a.d[12]*b.d[5] + a.d[18]*b.d[8] + a.d[24]*b.d[12] + a.d[30]*b.d[17];
    c.d[13] = a.d[1]*b.d[3] + a.d[7]*b.d[4] + a.d[13]*b.d[5] + a.d[19]*b.d[8] + a.d[25]*b.d[12] + a.d[31]*b.d[17];
    c.d[14] = a.d[2]*b.d[3] + a.d[8]*b.d[4] + a.d[14]*b.d[5] + a.d[20]*b.d[8] + a.d[26]*b.d[12] + a.d[32]*b.d[17];
    c.d[15] = a.d[3]*b.d[3] + a.d[9]*b.d[4] + a.d[15]*b.d[5] + a.d[21]*b.d[8] + a.d[27]*b.d[12] + a.d[33]*b.d[17];
    c.d[16] = a.d[4]*b.d[3] + a.d[10]*b.d[4] + a.d[16]*b.d[5] + a.d[22]*b.d[8] + a.d[28]*b.d[12] + a.d[34]*b.d[17];
    c.d[17] = a.d[5]*b.d[3] + a.d[11]*b.d[4] + a.d[17]*b.d[5] + a.d[23]*b.d[8] + a.d[29]*b.d[12] + a.d[35]*b.d[17];
    c.d[18] = a.d[0]*b.d[6] + a.d[6]*b.d[7] + a.d[12]*b.d[8] + a.d[18]*b.d[9] + a.d[24]*b.d[13] + a.d[30]*b.d[18];
    c.d[19] = a.d[1]*b.d[6] + a.d[7]*b.d[7] + a.d[13]*b.d[8] + a.d[19]*b.d[9] + a.d[25]*b.d[13] + a.d[31]*b.d[18];
    c.d[20] = a.d[2]*b.d[6] + a.d[8]*b.d[7] + a.d[14]*b.d[8] + a.d[20]*b.d[9] + a.d[26]*b.d[13] + a.d[32]*b.d[18];
    c.d[21] = a.d[3]*b.d[6] + a.d[9]*b.d[7] + a.d[15]*b.d[8] + a.d[21]*b.d[9] + a.d[27]*b.d[13] + a.d[33]*b.d[18];
    c.d[22] = a.d[4]*b.d[6] + a.d[10]*b.d[7] + a.d[16]*b.d[8] + a.d[22]*b.d[9] + a.d[28]*b.d[13] + a.d[34]*b.d[18];
    c.d[23] = a.d[5]*b.d[6] + a.d[11]*b.d[7] + a.d[17]*b.d[8] + a.d[23]*b.d[9] + a.d[29]*b.d[13] + a.d[35]*b.d[18];
    c.d[24] = a.d[0]*b.d[10] + a.d[6]*b.d[11] + a.d[12]*b.d[12] + a.d[18]*b.d[13] + a.d[24]*b.d[14] + a.d[30]*b.d[19];
    c.d[25] = a.d[1]*b.d[10] + a.d[7]*b.d[11] + a.d[13]*b.d[12] + a.d[19]*b.d[13] + a.d[25]*b.d[14] + a.d[31]*b.d[19];
    c.d[26] = a.d[2]*b.d[10] + a.d[8]*b.d[11] + a.d[14]*b.d[12] + a.d[20]*b.d[13] + a.d[26]*b.d[14] + a.d[32]*b.d[19];
    c.d[27] = a.d[3]*b.d[10] + a.d[9]*b.d[11] + a.d[15]*b.d[12] + a.d[21]*b.d[13] + a.d[27]*b.d[14] + a.d[33]*b.d[19];
    c.d[28] = a.d[4]*b.d[10] + a.d[10]*b.d[11] + a.d[16]*b.d[12] + a.d[22]*b.d[13] + a.d[28]*b.d[14] + a.d[34]*b.d[19];
    c.d[29] = a.d[5]*b.d[10] + a.d[11]*b.d[11] + a.d[17]*b.d[12] + a.d[23]*b.d[13] + a.d[29]*b.d[14] + a.d[35]*b.d[19];
    c.d[30] = a.d[0]*b.d[15] + a.d[6]*b.d[16] + a.d[12]*b.d[17] + a.d[18]*b.d[18] + a.d[24]*b.d[19] + a.d[30]*b.d[20];
    c.d[31] = a.d[1]*b.d[15] + a.d[7]*b.d[16] + a.d[13]*b.d[17] + a.d[19]*b.d[18] + a.d[25]*b.d[19] + a.d[31]*b.d[20];
    c.d[32] = a.d[2]*b.d[15] + a.d[8]*b.d[16] + a.d[14]*b.d[17] + a.d[20]*b.d[18] + a.d[26]*b.d[19] + a.d[32]*b.d[20];
    c.d[33] = a.d[3]*b.d[15] + a.d[9]*b.d[16] + a.d[15]*b.d[17] + a.d[21]*b.d[18] + a.d[27]*b.d[19] + a.d[33]*b.d[20];
    c.d[34] = a.d[4]*b.d[15] + a.d[10]*b.d[16] + a.d[16]*b.d[17] + a.d[22]*b.d[18] + a.d[28]*b.d[19] + a.d[34]*b.d[20];
    c.d[35] = a.d[5]*b.d[15] + a.d[11]*b.d[16] + a.d[17]*b.d[17] + a.d[23]*b.d[18] + a.d[29]*b.d[19] + a.d[35]*b.d[20];
    
    return c;
    
}

//-----------------------------------------------------------------------------
// inverse
inline tens4dmm tens4dmm::inverse() const
{
    matrix c(6,6);
    
    // populate c
    c(0,0) = d[ 0]; c(0,1) = d[ 6]; c(0,2) = d[12]; c(0,3) = d[18]; c(0,4) = d[24]; c(0,5) = d[30];
    c(1,0) = d[ 1]; c(1,1) = d[ 7]; c(1,2) = d[13]; c(1,3) = d[19]; c(1,4) = d[25]; c(1,5) = d[31];
    c(2,0) = d[ 2]; c(2,1) = d[ 8]; c(2,2) = d[14]; c(2,3) = d[20]; c(2,4) = d[26]; c(2,5) = d[32];
    c(3,0) = d[ 3]; c(3,1) = d[ 9]; c(3,2) = d[15]; c(3,3) = d[21]; c(3,4) = d[27]; c(3,5) = d[33];
    c(4,0) = d[ 4]; c(4,1) = d[10]; c(4,2) = d[16]; c(4,3) = d[22]; c(4,4) = d[28]; c(4,5) = d[34];
    c(5,0) = d[ 5]; c(5,1) = d[11]; c(5,2) = d[17]; c(5,3) = d[23]; c(5,4) = d[29]; c(5,5) = d[35];
    
    // invert c
    matrix s = c.inverse();
    
    // return inverse
    tens4dmm S;
    S.d[ 0] = s(0,0); S.d[ 6] = s(0,1); S.d[12] = s(0,2); S.d[18] = s(0,3); S.d[24] = s(0,4); S.d[30] = s(0,5);
    S.d[ 1] = s(1,0); S.d[ 7] = s(1,1); S.d[13] = s(1,2); S.d[19] = s(1,3); S.d[25] = s(1,4); S.d[31] = s(1,5);
    S.d[ 2] = s(2,0); S.d[ 8] = s(2,1); S.d[14] = s(2,2); S.d[20] = s(2,3); S.d[26] = s(2,4); S.d[32] = s(2,5);
    S.d[ 3] = s(3,0); S.d[ 9] = s(3,1); S.d[15] = s(3,2); S.d[21] = s(3,3); S.d[27] = s(3,4); S.d[33] = s(3,5);
    S.d[ 4] = s(4,0); S.d[10] = s(4,1); S.d[16] = s(4,2); S.d[22] = s(4,3); S.d[28] = s(4,4); S.d[34] = s(4,5);
    S.d[ 5] = s(5,0); S.d[11] = s(5,1); S.d[17] = s(5,2); S.d[23] = s(5,3); S.d[29] = s(5,4); S.d[35] = s(5,5);

    return S;
}

//-----------------------------------------------------------------------------
// evaluate push/pull operation
// c_ijpq = F_ik F_jl C_klmn F_pm F_qn
inline tens4dmm tens4dmm::pp(const mat3d& F)
{
    tens4dmm c;

    c.d[0] = d[0]*pow(F(0,0),4) + 2*d[3]*pow(F(0,0),3)*F(0,1) + 2*d[18]*pow(F(0,0),3)*F(0,1) + d[1]*pow(F(0,0),2)*pow(F(0,1),2) + d[6]*pow(F(0,0),2)*pow(F(0,1),2) + 4*d[21]*pow(F(0,0),2)*pow(F(0,1),2) + 2*d[9]*F(0,0)*pow(F(0,1),3) + 2*d[19]*F(0,0)*pow(F(0,1),3) + d[7]*pow(F(0,1),4) + 2*d[5]*pow(F(0,0),3)*F(0,2) + 2*d[30]*pow(F(0,0),3)*F(0,2) + 2*d[4]*pow(F(0,0),2)*F(0,1)*F(0,2) + 4*d[23]*pow(F(0,0),2)*F(0,1)*F(0,2) + 2*d[24]*pow(F(0,0),2)*F(0,1)*F(0,2) + 4*d[33]*pow(F(0,0),2)*F(0,1)*F(0,2) + 2*d[11]*F(0,0)*pow(F(0,1),2)*F(0,2) + 4*d[22]*F(0,0)*pow(F(0,1),2)*F(0,2) + 4*d[27]*F(0,0)*pow(F(0,1),2)*F(0,2) + 2*d[31]*F(0,0)*pow(F(0,1),2)*F(0,2) + 2*d[10]*pow(F(0,1),3)*F(0,2) + 2*d[25]*pow(F(0,1),3)*F(0,2) + d[2]*pow(F(0,0),2)*pow(F(0,2),2) + d[12]*pow(F(0,0),2)*pow(F(0,2),2) + 4*d[35]*pow(F(0,0),2)*pow(F(0,2),2) + 2*d[15]*F(0,0)*F(0,1)*pow(F(0,2),2) + 2*d[20]*F(0,0)*F(0,1)*pow(F(0,2),2) + 4*d[29]*F(0,0)*F(0,1)*pow(F(0,2),2) + 4*d[34]*F(0,0)*F(0,1)*pow(F(0,2),2) + d[8]*pow(F(0,1),2)*pow(F(0,2),2) + d[13]*pow(F(0,1),2)*pow(F(0,2),2) + 4*d[28]*pow(F(0,1),2)*pow(F(0,2),2) + 2*d[17]*F(0,0)*pow(F(0,2),3) + 2*d[32]*F(0,0)*pow(F(0,2),3) + 2*d[16]*F(0,1)*pow(F(0,2),3) + 2*d[26]*F(0,1)*pow(F(0,2),3) + d[14]*pow(F(0,2),4);

    c.d[1] = d[0]*pow(F(0,0),2)*pow(F(1,0),2) + 2*d[18]*F(0,0)*F(0,1)*pow(F(1,0),2) + d[6]*pow(F(0,1),2)*pow(F(1,0),2) + 2*d[30]*F(0,0)*F(0,2)*pow(F(1,0),2) + 2*d[24]*F(0,1)*F(0,2)*pow(F(1,0),2) + d[12]*pow(F(0,2),2)*pow(F(1,0),2) + 2*d[3]*pow(F(0,0),2)*F(1,0)*F(1,1) + 4*d[21]*F(0,0)*F(0,1)*F(1,0)*F(1,1) + 2*d[9]*pow(F(0,1),2)*F(1,0)*F(1,1) + 4*d[33]*F(0,0)*F(0,2)*F(1,0)*F(1,1) + 4*d[27]*F(0,1)*F(0,2)*F(1,0)*F(1,1) + 2*d[15]*pow(F(0,2),2)*F(1,0)*F(1,1) + d[1]*pow(F(0,0),2)*pow(F(1,1),2) + 2*d[19]*F(0,0)*F(0,1)*pow(F(1,1),2) + d[7]*pow(F(0,1),2)*pow(F(1,1),2) + 2*d[31]*F(0,0)*F(0,2)*pow(F(1,1),2) + 2*d[25]*F(0,1)*F(0,2)*pow(F(1,1),2) + d[13]*pow(F(0,2),2)*pow(F(1,1),2) + 2*d[5]*pow(F(0,0),2)*F(1,0)*F(1,2) + 4*d[23]*F(0,0)*F(0,1)*F(1,0)*F(1,2) + 2*d[11]*pow(F(0,1),2)*F(1,0)*F(1,2) + 4*d[35]*F(0,0)*F(0,2)*F(1,0)*F(1,2) + 4*d[29]*F(0,1)*F(0,2)*F(1,0)*F(1,2) + 2*d[17]*pow(F(0,2),2)*F(1,0)*F(1,2) + 2*d[4]*pow(F(0,0),2)*F(1,1)*F(1,2) + 4*d[22]*F(0,0)*F(0,1)*F(1,1)*F(1,2) + 2*d[10]*pow(F(0,1),2)*F(1,1)*F(1,2) + 4*d[34]*F(0,0)*F(0,2)*F(1,1)*F(1,2) + 4*d[28]*F(0,1)*F(0,2)*F(1,1)*F(1,2) + 2*d[16]*pow(F(0,2),2)*F(1,1)*F(1,2) + d[2]*pow(F(0,0),2)*pow(F(1,2),2) + 2*d[20]*F(0,0)*F(0,1)*pow(F(1,2),2) + d[8]*pow(F(0,1),2)*pow(F(1,2),2) + 2*d[32]*F(0,0)*F(0,2)*pow(F(1,2),2) + 2*d[26]*F(0,1)*F(0,2)*pow(F(1,2),2) + d[14]*pow(F(0,2),2)*pow(F(1,2),2);

    c.d[2] = d[0]*pow(F(0,0),2)*pow(F(2,0),2) + 2*d[18]*F(0,0)*F(0,1)*pow(F(2,0),2) + d[6]*pow(F(0,1),2)*pow(F(2,0),2) + 2*d[30]*F(0,0)*F(0,2)*pow(F(2,0),2) + 2*d[24]*F(0,1)*F(0,2)*pow(F(2,0),2) + d[12]*pow(F(0,2),2)*pow(F(2,0),2) + 2*d[3]*pow(F(0,0),2)*F(2,0)*F(2,1) + 4*d[21]*F(0,0)*F(0,1)*F(2,0)*F(2,1) + 2*d[9]*pow(F(0,1),2)*F(2,0)*F(2,1) + 4*d[33]*F(0,0)*F(0,2)*F(2,0)*F(2,1) + 4*d[27]*F(0,1)*F(0,2)*F(2,0)*F(2,1) + 2*d[15]*pow(F(0,2),2)*F(2,0)*F(2,1) + d[1]*pow(F(0,0),2)*pow(F(2,1),2) + 2*d[19]*F(0,0)*F(0,1)*pow(F(2,1),2) + d[7]*pow(F(0,1),2)*pow(F(2,1),2) + 2*d[31]*F(0,0)*F(0,2)*pow(F(2,1),2) + 2*d[25]*F(0,1)*F(0,2)*pow(F(2,1),2) + d[13]*pow(F(0,2),2)*pow(F(2,1),2) + 2*d[5]*pow(F(0,0),2)*F(2,0)*F(2,2) + 4*d[23]*F(0,0)*F(0,1)*F(2,0)*F(2,2) + 2*d[11]*pow(F(0,1),2)*F(2,0)*F(2,2) + 4*d[35]*F(0,0)*F(0,2)*F(2,0)*F(2,2) + 4*d[29]*F(0,1)*F(0,2)*F(2,0)*F(2,2) + 2*d[17]*pow(F(0,2),2)*F(2,0)*F(2,2) + 2*d[4]*pow(F(0,0),2)*F(2,1)*F(2,2) + 4*d[22]*F(0,0)*F(0,1)*F(2,1)*F(2,2) + 2*d[10]*pow(F(0,1),2)*F(2,1)*F(2,2) + 4*d[34]*F(0,0)*F(0,2)*F(2,1)*F(2,2) + 4*d[28]*F(0,1)*F(0,2)*F(2,1)*F(2,2) + 2*d[16]*pow(F(0,2),2)*F(2,1)*F(2,2) + d[2]*pow(F(0,0),2)*pow(F(2,2),2) + 2*d[20]*F(0,0)*F(0,1)*pow(F(2,2),2) + d[8]*pow(F(0,1),2)*pow(F(2,2),2) + 2*d[32]*F(0,0)*F(0,2)*pow(F(2,2),2) + 2*d[26]*F(0,1)*F(0,2)*pow(F(2,2),2) + d[14]*pow(F(0,2),2)*pow(F(2,2),2);

    c.d[3] = d[0]*pow(F(0,0),3)*F(1,0) + 2*d[18]*pow(F(0,0),2)*F(0,1)*F(1,0) + d[6]*F(0,0)*pow(F(0,1),2)*F(1,0) + 2*d[21]*F(0,0)*pow(F(0,1),2)*F(1,0) + d[9]*pow(F(0,1),3)*F(1,0) + d[5]*pow(F(0,0),2)*F(0,2)*F(1,0) + 2*d[30]*pow(F(0,0),2)*F(0,2)*F(1,0) + 2*d[23]*F(0,0)*F(0,1)*F(0,2)*F(1,0) + 2*d[24]*F(0,0)*F(0,1)*F(0,2)*F(1,0) + 2*d[33]*F(0,0)*F(0,1)*F(0,2)*F(1,0) + d[11]*pow(F(0,1),2)*F(0,2)*F(1,0) + 2*d[27]*pow(F(0,1),2)*F(0,2)*F(1,0) + d[12]*F(0,0)*pow(F(0,2),2)*F(1,0) + 2*d[35]*F(0,0)*pow(F(0,2),2)*F(1,0) + d[15]*F(0,1)*pow(F(0,2),2)*F(1,0) + 2*d[29]*F(0,1)*pow(F(0,2),2)*F(1,0) + d[17]*pow(F(0,2),3)*F(1,0) + d[1]*pow(F(0,0),2)*F(0,1)*F(1,1) + 2*d[21]*pow(F(0,0),2)*F(0,1)*F(1,1) + d[9]*F(0,0)*pow(F(0,1),2)*F(1,1) + 2*d[19]*F(0,0)*pow(F(0,1),2)*F(1,1) + d[7]*pow(F(0,1),3)*F(1,1) + d[4]*pow(F(0,0),2)*F(0,2)*F(1,1) + 2*d[33]*pow(F(0,0),2)*F(0,2)*F(1,1) + 2*d[22]*F(0,0)*F(0,1)*F(0,2)*F(1,1) + 2*d[27]*F(0,0)*F(0,1)*F(0,2)*F(1,1) + 2*d[31]*F(0,0)*F(0,1)*F(0,2)*F(1,1) + d[10]*pow(F(0,1),2)*F(0,2)*F(1,1) + 2*d[25]*pow(F(0,1),2)*F(0,2)*F(1,1) + d[15]*F(0,0)*pow(F(0,2),2)*F(1,1) + 2*d[34]*F(0,0)*pow(F(0,2),2)*F(1,1) + d[13]*F(0,1)*pow(F(0,2),2)*F(1,1) + 2*d[28]*F(0,1)*pow(F(0,2),2)*F(1,1) + d[16]*pow(F(0,2),3)*F(1,1) + d[3]*pow(F(0,0),2)*(F(0,1)*F(1,0) + F(0,0)*F(1,1)) + d[5]*pow(F(0,0),3)*F(1,2) + d[4]*pow(F(0,0),2)*F(0,1)*F(1,2) + 2*d[23]*pow(F(0,0),2)*F(0,1)*F(1,2) + d[11]*F(0,0)*pow(F(0,1),2)*F(1,2) + 2*d[22]*F(0,0)*pow(F(0,1),2)*F(1,2) + d[10]*pow(F(0,1),3)*F(1,2) + d[2]*pow(F(0,0),2)*F(0,2)*F(1,2) + 2*d[35]*pow(F(0,0),2)*F(0,2)*F(1,2) + 2*d[20]*F(0,0)*F(0,1)*F(0,2)*F(1,2) + 2*d[29]*F(0,0)*F(0,1)*F(0,2)*F(1,2) + 2*d[34]*F(0,0)*F(0,1)*F(0,2)*F(1,2) + d[8]*pow(F(0,1),2)*F(0,2)*F(1,2) + 2*d[28]*pow(F(0,1),2)*F(0,2)*F(1,2) + d[17]*F(0,0)*pow(F(0,2),2)*F(1,2) + 2*d[32]*F(0,0)*pow(F(0,2),2)*F(1,2) + d[16]*F(0,1)*pow(F(0,2),2)*F(1,2) + 2*d[26]*F(0,1)*pow(F(0,2),2)*F(1,2) + d[14]*pow(F(0,2),3)*F(1,2);

    c.d[4] = d[0]*pow(F(0,0),2)*F(1,0)*F(2,0) + 2*d[18]*F(0,0)*F(0,1)*F(1,0)*F(2,0) + d[6]*pow(F(0,1),2)*F(1,0)*F(2,0) + 2*d[30]*F(0,0)*F(0,2)*F(1,0)*F(2,0) + 2*d[24]*F(0,1)*F(0,2)*F(1,0)*F(2,0) + d[12]*pow(F(0,2),2)*F(1,0)*F(2,0) + d[3]*pow(F(0,0),2)*F(1,1)*F(2,0) + 2*d[21]*F(0,0)*F(0,1)*F(1,1)*F(2,0) + d[9]*pow(F(0,1),2)*F(1,1)*F(2,0) + 2*d[33]*F(0,0)*F(0,2)*F(1,1)*F(2,0) + 2*d[27]*F(0,1)*F(0,2)*F(1,1)*F(2,0) + d[15]*pow(F(0,2),2)*F(1,1)*F(2,0) + d[5]*pow(F(0,0),2)*F(1,2)*F(2,0) + 2*d[23]*F(0,0)*F(0,1)*F(1,2)*F(2,0) + d[11]*pow(F(0,1),2)*F(1,2)*F(2,0) + 2*d[35]*F(0,0)*F(0,2)*F(1,2)*F(2,0) + 2*d[29]*F(0,1)*F(0,2)*F(1,2)*F(2,0) + d[17]*pow(F(0,2),2)*F(1,2)*F(2,0) + d[3]*pow(F(0,0),2)*F(1,0)*F(2,1) + 2*d[21]*F(0,0)*F(0,1)*F(1,0)*F(2,1) + d[9]*pow(F(0,1),2)*F(1,0)*F(2,1) + 2*d[33]*F(0,0)*F(0,2)*F(1,0)*F(2,1) + 2*d[27]*F(0,1)*F(0,2)*F(1,0)*F(2,1) + d[15]*pow(F(0,2),2)*F(1,0)*F(2,1) + d[1]*pow(F(0,0),2)*F(1,1)*F(2,1) + 2*d[19]*F(0,0)*F(0,1)*F(1,1)*F(2,1) + d[7]*pow(F(0,1),2)*F(1,1)*F(2,1) + 2*d[31]*F(0,0)*F(0,2)*F(1,1)*F(2,1) + 2*d[25]*F(0,1)*F(0,2)*F(1,1)*F(2,1) + d[13]*pow(F(0,2),2)*F(1,1)*F(2,1) + d[4]*pow(F(0,0),2)*F(1,2)*F(2,1) + 2*d[22]*F(0,0)*F(0,1)*F(1,2)*F(2,1) + d[10]*pow(F(0,1),2)*F(1,2)*F(2,1) + 2*d[34]*F(0,0)*F(0,2)*F(1,2)*F(2,1) + 2*d[28]*F(0,1)*F(0,2)*F(1,2)*F(2,1) + d[16]*pow(F(0,2),2)*F(1,2)*F(2,1) + d[5]*pow(F(0,0),2)*F(1,0)*F(2,2) + 2*d[23]*F(0,0)*F(0,1)*F(1,0)*F(2,2) + d[11]*pow(F(0,1),2)*F(1,0)*F(2,2) + 2*d[35]*F(0,0)*F(0,2)*F(1,0)*F(2,2) + 2*d[29]*F(0,1)*F(0,2)*F(1,0)*F(2,2) + d[17]*pow(F(0,2),2)*F(1,0)*F(2,2) + d[4]*pow(F(0,0),2)*F(1,1)*F(2,2) + 2*d[22]*F(0,0)*F(0,1)*F(1,1)*F(2,2) + d[10]*pow(F(0,1),2)*F(1,1)*F(2,2) + 2*d[34]*F(0,0)*F(0,2)*F(1,1)*F(2,2) + 2*d[28]*F(0,1)*F(0,2)*F(1,1)*F(2,2) + d[16]*pow(F(0,2),2)*F(1,1)*F(2,2) + d[2]*pow(F(0,0),2)*F(1,2)*F(2,2) + 2*d[20]*F(0,0)*F(0,1)*F(1,2)*F(2,2) + d[8]*pow(F(0,1),2)*F(1,2)*F(2,2) + 2*d[32]*F(0,0)*F(0,2)*F(1,2)*F(2,2) + 2*d[26]*F(0,1)*F(0,2)*F(1,2)*F(2,2) + d[14]*pow(F(0,2),2)*F(1,2)*F(2,2);

    c.d[5] = d[0]*pow(F(0,0),3)*F(2,0) + 2*d[18]*pow(F(0,0),2)*F(0,1)*F(2,0) + d[6]*F(0,0)*pow(F(0,1),2)*F(2,0) + 2*d[21]*F(0,0)*pow(F(0,1),2)*F(2,0) + d[9]*pow(F(0,1),3)*F(2,0) + d[5]*pow(F(0,0),2)*F(0,2)*F(2,0) + 2*d[30]*pow(F(0,0),2)*F(0,2)*F(2,0) + 2*d[23]*F(0,0)*F(0,1)*F(0,2)*F(2,0) + 2*d[24]*F(0,0)*F(0,1)*F(0,2)*F(2,0) + 2*d[33]*F(0,0)*F(0,1)*F(0,2)*F(2,0) + d[11]*pow(F(0,1),2)*F(0,2)*F(2,0) + 2*d[27]*pow(F(0,1),2)*F(0,2)*F(2,0) + d[12]*F(0,0)*pow(F(0,2),2)*F(2,0) + 2*d[35]*F(0,0)*pow(F(0,2),2)*F(2,0) + d[15]*F(0,1)*pow(F(0,2),2)*F(2,0) + 2*d[29]*F(0,1)*pow(F(0,2),2)*F(2,0) + d[17]*pow(F(0,2),3)*F(2,0) + d[1]*pow(F(0,0),2)*F(0,1)*F(2,1) + 2*d[21]*pow(F(0,0),2)*F(0,1)*F(2,1) + d[9]*F(0,0)*pow(F(0,1),2)*F(2,1) + 2*d[19]*F(0,0)*pow(F(0,1),2)*F(2,1) + d[7]*pow(F(0,1),3)*F(2,1) + d[4]*pow(F(0,0),2)*F(0,2)*F(2,1) + 2*d[33]*pow(F(0,0),2)*F(0,2)*F(2,1) + 2*d[22]*F(0,0)*F(0,1)*F(0,2)*F(2,1) + 2*d[27]*F(0,0)*F(0,1)*F(0,2)*F(2,1) + 2*d[31]*F(0,0)*F(0,1)*F(0,2)*F(2,1) + d[10]*pow(F(0,1),2)*F(0,2)*F(2,1) + 2*d[25]*pow(F(0,1),2)*F(0,2)*F(2,1) + d[15]*F(0,0)*pow(F(0,2),2)*F(2,1) + 2*d[34]*F(0,0)*pow(F(0,2),2)*F(2,1) + d[13]*F(0,1)*pow(F(0,2),2)*F(2,1) + 2*d[28]*F(0,1)*pow(F(0,2),2)*F(2,1) + d[16]*pow(F(0,2),3)*F(2,1) + d[3]*pow(F(0,0),2)*(F(0,1)*F(2,0) + F(0,0)*F(2,1)) + d[5]*pow(F(0,0),3)*F(2,2) + d[4]*pow(F(0,0),2)*F(0,1)*F(2,2) + 2*d[23]*pow(F(0,0),2)*F(0,1)*F(2,2) + d[11]*F(0,0)*pow(F(0,1),2)*F(2,2) + 2*d[22]*F(0,0)*pow(F(0,1),2)*F(2,2) + d[10]*pow(F(0,1),3)*F(2,2) + d[2]*pow(F(0,0),2)*F(0,2)*F(2,2) + 2*d[35]*pow(F(0,0),2)*F(0,2)*F(2,2) + 2*d[20]*F(0,0)*F(0,1)*F(0,2)*F(2,2) + 2*d[29]*F(0,0)*F(0,1)*F(0,2)*F(2,2) + 2*d[34]*F(0,0)*F(0,1)*F(0,2)*F(2,2) + d[8]*pow(F(0,1),2)*F(0,2)*F(2,2) + 2*d[28]*pow(F(0,1),2)*F(0,2)*F(2,2) + d[17]*F(0,0)*pow(F(0,2),2)*F(2,2) + 2*d[32]*F(0,0)*pow(F(0,2),2)*F(2,2) + d[16]*F(0,1)*pow(F(0,2),2)*F(2,2) + 2*d[26]*F(0,1)*pow(F(0,2),2)*F(2,2) + d[14]*pow(F(0,2),3)*F(2,2);

    c.d[6] = d[0]*pow(F(0,0),2)*pow(F(1,0),2) + 2*d[3]*F(0,0)*F(0,1)*pow(F(1,0),2) + d[1]*pow(F(0,1),2)*pow(F(1,0),2) + 2*d[5]*F(0,0)*F(0,2)*pow(F(1,0),2) + 2*d[4]*F(0,1)*F(0,2)*pow(F(1,0),2) + d[2]*pow(F(0,2),2)*pow(F(1,0),2) + 2*d[18]*pow(F(0,0),2)*F(1,0)*F(1,1) + 4*d[21]*F(0,0)*F(0,1)*F(1,0)*F(1,1) + 2*d[19]*pow(F(0,1),2)*F(1,0)*F(1,1) + 4*d[23]*F(0,0)*F(0,2)*F(1,0)*F(1,1) + 4*d[22]*F(0,1)*F(0,2)*F(1,0)*F(1,1) + 2*d[20]*pow(F(0,2),2)*F(1,0)*F(1,1) + d[6]*pow(F(0,0),2)*pow(F(1,1),2) + 2*d[9]*F(0,0)*F(0,1)*pow(F(1,1),2) + d[7]*pow(F(0,1),2)*pow(F(1,1),2) + 2*d[11]*F(0,0)*F(0,2)*pow(F(1,1),2) + 2*d[10]*F(0,1)*F(0,2)*pow(F(1,1),2) + d[8]*pow(F(0,2),2)*pow(F(1,1),2) + 2*d[30]*pow(F(0,0),2)*F(1,0)*F(1,2) + 4*d[33]*F(0,0)*F(0,1)*F(1,0)*F(1,2) + 2*d[31]*pow(F(0,1),2)*F(1,0)*F(1,2) + 4*d[35]*F(0,0)*F(0,2)*F(1,0)*F(1,2) + 4*d[34]*F(0,1)*F(0,2)*F(1,0)*F(1,2) + 2*d[32]*pow(F(0,2),2)*F(1,0)*F(1,2) + 2*d[24]*pow(F(0,0),2)*F(1,1)*F(1,2) + 4*d[27]*F(0,0)*F(0,1)*F(1,1)*F(1,2) + 2*d[25]*pow(F(0,1),2)*F(1,1)*F(1,2) + 4*d[29]*F(0,0)*F(0,2)*F(1,1)*F(1,2) + 4*d[28]*F(0,1)*F(0,2)*F(1,1)*F(1,2) + 2*d[26]*pow(F(0,2),2)*F(1,1)*F(1,2) + d[12]*pow(F(0,0),2)*pow(F(1,2),2) + 2*d[15]*F(0,0)*F(0,1)*pow(F(1,2),2) + d[13]*pow(F(0,1),2)*pow(F(1,2),2) + 2*d[17]*F(0,0)*F(0,2)*pow(F(1,2),2) + 2*d[16]*F(0,1)*F(0,2)*pow(F(1,2),2) + d[14]*pow(F(0,2),2)*pow(F(1,2),2);

    c.d[7] = d[0]*pow(F(1,0),4) + 2*d[3]*pow(F(1,0),3)*F(1,1) + 2*d[18]*pow(F(1,0),3)*F(1,1) + d[1]*pow(F(1,0),2)*pow(F(1,1),2) + d[6]*pow(F(1,0),2)*pow(F(1,1),2) + 4*d[21]*pow(F(1,0),2)*pow(F(1,1),2) + 2*d[9]*F(1,0)*pow(F(1,1),3) + 2*d[19]*F(1,0)*pow(F(1,1),3) + d[7]*pow(F(1,1),4) + 2*d[5]*pow(F(1,0),3)*F(1,2) + 2*d[30]*pow(F(1,0),3)*F(1,2) + 2*d[4]*pow(F(1,0),2)*F(1,1)*F(1,2) + 4*d[23]*pow(F(1,0),2)*F(1,1)*F(1,2) + 2*d[24]*pow(F(1,0),2)*F(1,1)*F(1,2) + 4*d[33]*pow(F(1,0),2)*F(1,1)*F(1,2) + 2*d[11]*F(1,0)*pow(F(1,1),2)*F(1,2) + 4*d[22]*F(1,0)*pow(F(1,1),2)*F(1,2) + 4*d[27]*F(1,0)*pow(F(1,1),2)*F(1,2) + 2*d[31]*F(1,0)*pow(F(1,1),2)*F(1,2) + 2*d[10]*pow(F(1,1),3)*F(1,2) + 2*d[25]*pow(F(1,1),3)*F(1,2) + d[2]*pow(F(1,0),2)*pow(F(1,2),2) + d[12]*pow(F(1,0),2)*pow(F(1,2),2) + 4*d[35]*pow(F(1,0),2)*pow(F(1,2),2) + 2*d[15]*F(1,0)*F(1,1)*pow(F(1,2),2) + 2*d[20]*F(1,0)*F(1,1)*pow(F(1,2),2) + 4*d[29]*F(1,0)*F(1,1)*pow(F(1,2),2) + 4*d[34]*F(1,0)*F(1,1)*pow(F(1,2),2) + d[8]*pow(F(1,1),2)*pow(F(1,2),2) + d[13]*pow(F(1,1),2)*pow(F(1,2),2) + 4*d[28]*pow(F(1,1),2)*pow(F(1,2),2) + 2*d[17]*F(1,0)*pow(F(1,2),3) + 2*d[32]*F(1,0)*pow(F(1,2),3) + 2*d[16]*F(1,1)*pow(F(1,2),3) + 2*d[26]*F(1,1)*pow(F(1,2),3) + d[14]*pow(F(1,2),4);

    c.d[8] = d[0]*pow(F(1,0),2)*pow(F(2,0),2) + 2*d[18]*F(1,0)*F(1,1)*pow(F(2,0),2) + d[6]*pow(F(1,1),2)*pow(F(2,0),2) + 2*d[30]*F(1,0)*F(1,2)*pow(F(2,0),2) + 2*d[24]*F(1,1)*F(1,2)*pow(F(2,0),2) + d[12]*pow(F(1,2),2)*pow(F(2,0),2) + 2*d[3]*pow(F(1,0),2)*F(2,0)*F(2,1) + 4*d[21]*F(1,0)*F(1,1)*F(2,0)*F(2,1) + 2*d[9]*pow(F(1,1),2)*F(2,0)*F(2,1) + 4*d[33]*F(1,0)*F(1,2)*F(2,0)*F(2,1) + 4*d[27]*F(1,1)*F(1,2)*F(2,0)*F(2,1) + 2*d[15]*pow(F(1,2),2)*F(2,0)*F(2,1) + d[1]*pow(F(1,0),2)*pow(F(2,1),2) + 2*d[19]*F(1,0)*F(1,1)*pow(F(2,1),2) + d[7]*pow(F(1,1),2)*pow(F(2,1),2) + 2*d[31]*F(1,0)*F(1,2)*pow(F(2,1),2) + 2*d[25]*F(1,1)*F(1,2)*pow(F(2,1),2) + d[13]*pow(F(1,2),2)*pow(F(2,1),2) + 2*d[5]*pow(F(1,0),2)*F(2,0)*F(2,2) + 4*d[23]*F(1,0)*F(1,1)*F(2,0)*F(2,2) + 2*d[11]*pow(F(1,1),2)*F(2,0)*F(2,2) + 4*d[35]*F(1,0)*F(1,2)*F(2,0)*F(2,2) + 4*d[29]*F(1,1)*F(1,2)*F(2,0)*F(2,2) + 2*d[17]*pow(F(1,2),2)*F(2,0)*F(2,2) + 2*d[4]*pow(F(1,0),2)*F(2,1)*F(2,2) + 4*d[22]*F(1,0)*F(1,1)*F(2,1)*F(2,2) + 2*d[10]*pow(F(1,1),2)*F(2,1)*F(2,2) + 4*d[34]*F(1,0)*F(1,2)*F(2,1)*F(2,2) + 4*d[28]*F(1,1)*F(1,2)*F(2,1)*F(2,2) + 2*d[16]*pow(F(1,2),2)*F(2,1)*F(2,2) + d[2]*pow(F(1,0),2)*pow(F(2,2),2) + 2*d[20]*F(1,0)*F(1,1)*pow(F(2,2),2) + d[8]*pow(F(1,1),2)*pow(F(2,2),2) + 2*d[32]*F(1,0)*F(1,2)*pow(F(2,2),2) + 2*d[26]*F(1,1)*F(1,2)*pow(F(2,2),2) + d[14]*pow(F(1,2),2)*pow(F(2,2),2);

    c.d[9] = d[0]*F(0,0)*pow(F(1,0),3) + d[5]*F(0,2)*pow(F(1,0),3) + 2*d[18]*F(0,0)*pow(F(1,0),2)*F(1,1) + d[1]*F(0,1)*pow(F(1,0),2)*F(1,1) + 2*d[21]*F(0,1)*pow(F(1,0),2)*F(1,1) + d[4]*F(0,2)*pow(F(1,0),2)*F(1,1) + 2*d[23]*F(0,2)*pow(F(1,0),2)*F(1,1) + d[6]*F(0,0)*F(1,0)*pow(F(1,1),2) + 2*d[21]*F(0,0)*F(1,0)*pow(F(1,1),2) + d[9]*F(0,1)*F(1,0)*pow(F(1,1),2) + 2*d[19]*F(0,1)*F(1,0)*pow(F(1,1),2) + d[11]*F(0,2)*F(1,0)*pow(F(1,1),2) + 2*d[22]*F(0,2)*F(1,0)*pow(F(1,1),2) + d[9]*F(0,0)*pow(F(1,1),3) + d[7]*F(0,1)*pow(F(1,1),3) + d[10]*F(0,2)*pow(F(1,1),3) + d[3]*pow(F(1,0),2)*(F(0,1)*F(1,0) + F(0,0)*F(1,1)) + d[5]*F(0,0)*pow(F(1,0),2)*F(1,2) + 2*d[30]*F(0,0)*pow(F(1,0),2)*F(1,2) + d[4]*F(0,1)*pow(F(1,0),2)*F(1,2) + 2*d[33]*F(0,1)*pow(F(1,0),2)*F(1,2) + d[2]*F(0,2)*pow(F(1,0),2)*F(1,2) + 2*d[35]*F(0,2)*pow(F(1,0),2)*F(1,2) + 2*d[23]*F(0,0)*F(1,0)*F(1,1)*F(1,2) + 2*d[24]*F(0,0)*F(1,0)*F(1,1)*F(1,2) + 2*d[33]*F(0,0)*F(1,0)*F(1,1)*F(1,2) + 2*d[22]*F(0,1)*F(1,0)*F(1,1)*F(1,2) + 2*d[27]*F(0,1)*F(1,0)*F(1,1)*F(1,2) + 2*d[31]*F(0,1)*F(1,0)*F(1,1)*F(1,2) + 2*d[20]*F(0,2)*F(1,0)*F(1,1)*F(1,2) + 2*d[29]*F(0,2)*F(1,0)*F(1,1)*F(1,2) + 2*d[34]*F(0,2)*F(1,0)*F(1,1)*F(1,2) + d[11]*F(0,0)*pow(F(1,1),2)*F(1,2) + 2*d[27]*F(0,0)*pow(F(1,1),2)*F(1,2) + d[10]*F(0,1)*pow(F(1,1),2)*F(1,2) + 2*d[25]*F(0,1)*pow(F(1,1),2)*F(1,2) + d[8]*F(0,2)*pow(F(1,1),2)*F(1,2) + 2*d[28]*F(0,2)*pow(F(1,1),2)*F(1,2) + d[12]*F(0,0)*F(1,0)*pow(F(1,2),2) + 2*d[35]*F(0,0)*F(1,0)*pow(F(1,2),2) + d[15]*F(0,1)*F(1,0)*pow(F(1,2),2) + 2*d[34]*F(0,1)*F(1,0)*pow(F(1,2),2) + d[17]*F(0,2)*F(1,0)*pow(F(1,2),2) + 2*d[32]*F(0,2)*F(1,0)*pow(F(1,2),2) + d[15]*F(0,0)*F(1,1)*pow(F(1,2),2) + 2*d[29]*F(0,0)*F(1,1)*pow(F(1,2),2) + d[13]*F(0,1)*F(1,1)*pow(F(1,2),2) + 2*d[28]*F(0,1)*F(1,1)*pow(F(1,2),2) + d[16]*F(0,2)*F(1,1)*pow(F(1,2),2) + 2*d[26]*F(0,2)*F(1,1)*pow(F(1,2),2) + d[17]*F(0,0)*pow(F(1,2),3) + d[16]*F(0,1)*pow(F(1,2),3) + d[14]*F(0,2)*pow(F(1,2),3);

    c.d[10] = d[0]*pow(F(1,0),3)*F(2,0) + 2*d[18]*pow(F(1,0),2)*F(1,1)*F(2,0) + d[6]*F(1,0)*pow(F(1,1),2)*F(2,0) + 2*d[21]*F(1,0)*pow(F(1,1),2)*F(2,0) + d[9]*pow(F(1,1),3)*F(2,0) + d[5]*pow(F(1,0),2)*F(1,2)*F(2,0) + 2*d[30]*pow(F(1,0),2)*F(1,2)*F(2,0) + 2*d[23]*F(1,0)*F(1,1)*F(1,2)*F(2,0) + 2*d[24]*F(1,0)*F(1,1)*F(1,2)*F(2,0) + 2*d[33]*F(1,0)*F(1,1)*F(1,2)*F(2,0) + d[11]*pow(F(1,1),2)*F(1,2)*F(2,0) + 2*d[27]*pow(F(1,1),2)*F(1,2)*F(2,0) + d[12]*F(1,0)*pow(F(1,2),2)*F(2,0) + 2*d[35]*F(1,0)*pow(F(1,2),2)*F(2,0) + d[15]*F(1,1)*pow(F(1,2),2)*F(2,0) + 2*d[29]*F(1,1)*pow(F(1,2),2)*F(2,0) + d[17]*pow(F(1,2),3)*F(2,0) + d[1]*pow(F(1,0),2)*F(1,1)*F(2,1) + 2*d[21]*pow(F(1,0),2)*F(1,1)*F(2,1) + d[9]*F(1,0)*pow(F(1,1),2)*F(2,1) + 2*d[19]*F(1,0)*pow(F(1,1),2)*F(2,1) + d[7]*pow(F(1,1),3)*F(2,1) + d[4]*pow(F(1,0),2)*F(1,2)*F(2,1) + 2*d[33]*pow(F(1,0),2)*F(1,2)*F(2,1) + 2*d[22]*F(1,0)*F(1,1)*F(1,2)*F(2,1) + 2*d[27]*F(1,0)*F(1,1)*F(1,2)*F(2,1) + 2*d[31]*F(1,0)*F(1,1)*F(1,2)*F(2,1) + d[10]*pow(F(1,1),2)*F(1,2)*F(2,1) + 2*d[25]*pow(F(1,1),2)*F(1,2)*F(2,1) + d[15]*F(1,0)*pow(F(1,2),2)*F(2,1) + 2*d[34]*F(1,0)*pow(F(1,2),2)*F(2,1) + d[13]*F(1,1)*pow(F(1,2),2)*F(2,1) + 2*d[28]*F(1,1)*pow(F(1,2),2)*F(2,1) + d[16]*pow(F(1,2),3)*F(2,1) + d[3]*pow(F(1,0),2)*(F(1,1)*F(2,0) + F(1,0)*F(2,1)) + d[5]*pow(F(1,0),3)*F(2,2) + d[4]*pow(F(1,0),2)*F(1,1)*F(2,2) + 2*d[23]*pow(F(1,0),2)*F(1,1)*F(2,2) + d[11]*F(1,0)*pow(F(1,1),2)*F(2,2) + 2*d[22]*F(1,0)*pow(F(1,1),2)*F(2,2) + d[10]*pow(F(1,1),3)*F(2,2) + d[2]*pow(F(1,0),2)*F(1,2)*F(2,2) + 2*d[35]*pow(F(1,0),2)*F(1,2)*F(2,2) + 2*d[20]*F(1,0)*F(1,1)*F(1,2)*F(2,2) + 2*d[29]*F(1,0)*F(1,1)*F(1,2)*F(2,2) + 2*d[34]*F(1,0)*F(1,1)*F(1,2)*F(2,2) + d[8]*pow(F(1,1),2)*F(1,2)*F(2,2) + 2*d[28]*pow(F(1,1),2)*F(1,2)*F(2,2) + d[17]*F(1,0)*pow(F(1,2),2)*F(2,2) + 2*d[32]*F(1,0)*pow(F(1,2),2)*F(2,2) + d[16]*F(1,1)*pow(F(1,2),2)*F(2,2) + 2*d[26]*F(1,1)*pow(F(1,2),2)*F(2,2) + d[14]*pow(F(1,2),3)*F(2,2);

    c.d[11] = d[0]*F(0,0)*pow(F(1,0),2)*F(2,0) + d[5]*F(0,2)*pow(F(1,0),2)*F(2,0) + 2*d[18]*F(0,0)*F(1,0)*F(1,1)*F(2,0) + 2*d[21]*F(0,1)*F(1,0)*F(1,1)*F(2,0) + 2*d[23]*F(0,2)*F(1,0)*F(1,1)*F(2,0) + d[6]*F(0,0)*pow(F(1,1),2)*F(2,0) + d[9]*F(0,1)*pow(F(1,1),2)*F(2,0) + d[11]*F(0,2)*pow(F(1,1),2)*F(2,0) + 2*d[30]*F(0,0)*F(1,0)*F(1,2)*F(2,0) + 2*d[33]*F(0,1)*F(1,0)*F(1,2)*F(2,0) + 2*d[35]*F(0,2)*F(1,0)*F(1,2)*F(2,0) + 2*d[24]*F(0,0)*F(1,1)*F(1,2)*F(2,0) + 2*d[27]*F(0,1)*F(1,1)*F(1,2)*F(2,0) + 2*d[29]*F(0,2)*F(1,1)*F(1,2)*F(2,0) + d[12]*F(0,0)*pow(F(1,2),2)*F(2,0) + d[15]*F(0,1)*pow(F(1,2),2)*F(2,0) + d[17]*F(0,2)*pow(F(1,2),2)*F(2,0) + d[1]*F(0,1)*pow(F(1,0),2)*F(2,1) + d[4]*F(0,2)*pow(F(1,0),2)*F(2,1) + 2*d[21]*F(0,0)*F(1,0)*F(1,1)*F(2,1) + 2*d[19]*F(0,1)*F(1,0)*F(1,1)*F(2,1) + 2*d[22]*F(0,2)*F(1,0)*F(1,1)*F(2,1) + d[9]*F(0,0)*pow(F(1,1),2)*F(2,1) + d[7]*F(0,1)*pow(F(1,1),2)*F(2,1) + d[10]*F(0,2)*pow(F(1,1),2)*F(2,1) + 2*d[33]*F(0,0)*F(1,0)*F(1,2)*F(2,1) + 2*d[31]*F(0,1)*F(1,0)*F(1,2)*F(2,1) + 2*d[34]*F(0,2)*F(1,0)*F(1,2)*F(2,1) + 2*d[27]*F(0,0)*F(1,1)*F(1,2)*F(2,1) + 2*d[25]*F(0,1)*F(1,1)*F(1,2)*F(2,1) + 2*d[28]*F(0,2)*F(1,1)*F(1,2)*F(2,1) + d[15]*F(0,0)*pow(F(1,2),2)*F(2,1) + d[13]*F(0,1)*pow(F(1,2),2)*F(2,1) + d[16]*F(0,2)*pow(F(1,2),2)*F(2,1) + d[3]*pow(F(1,0),2)*(F(0,1)*F(2,0) + F(0,0)*F(2,1)) + d[5]*F(0,0)*pow(F(1,0),2)*F(2,2) + d[4]*F(0,1)*pow(F(1,0),2)*F(2,2) + d[2]*F(0,2)*pow(F(1,0),2)*F(2,2) + 2*d[23]*F(0,0)*F(1,0)*F(1,1)*F(2,2) + 2*d[22]*F(0,1)*F(1,0)*F(1,1)*F(2,2) + 2*d[20]*F(0,2)*F(1,0)*F(1,1)*F(2,2) + d[11]*F(0,0)*pow(F(1,1),2)*F(2,2) + d[10]*F(0,1)*pow(F(1,1),2)*F(2,2) + d[8]*F(0,2)*pow(F(1,1),2)*F(2,2) + 2*d[35]*F(0,0)*F(1,0)*F(1,2)*F(2,2) + 2*d[34]*F(0,1)*F(1,0)*F(1,2)*F(2,2) + 2*d[32]*F(0,2)*F(1,0)*F(1,2)*F(2,2) + 2*d[29]*F(0,0)*F(1,1)*F(1,2)*F(2,2) + 2*d[28]*F(0,1)*F(1,1)*F(1,2)*F(2,2) + 2*d[26]*F(0,2)*F(1,1)*F(1,2)*F(2,2) + d[17]*F(0,0)*pow(F(1,2),2)*F(2,2) + d[16]*F(0,1)*pow(F(1,2),2)*F(2,2) + d[14]*F(0,2)*pow(F(1,2),2)*F(2,2);

    c.d[12] = d[0]*pow(F(0,0),2)*pow(F(2,0),2) + 2*d[3]*F(0,0)*F(0,1)*pow(F(2,0),2) + d[1]*pow(F(0,1),2)*pow(F(2,0),2) + 2*d[5]*F(0,0)*F(0,2)*pow(F(2,0),2) + 2*d[4]*F(0,1)*F(0,2)*pow(F(2,0),2) + d[2]*pow(F(0,2),2)*pow(F(2,0),2) + 2*d[18]*pow(F(0,0),2)*F(2,0)*F(2,1) + 4*d[21]*F(0,0)*F(0,1)*F(2,0)*F(2,1) + 2*d[19]*pow(F(0,1),2)*F(2,0)*F(2,1) + 4*d[23]*F(0,0)*F(0,2)*F(2,0)*F(2,1) + 4*d[22]*F(0,1)*F(0,2)*F(2,0)*F(2,1) + 2*d[20]*pow(F(0,2),2)*F(2,0)*F(2,1) + d[6]*pow(F(0,0),2)*pow(F(2,1),2) + 2*d[9]*F(0,0)*F(0,1)*pow(F(2,1),2) + d[7]*pow(F(0,1),2)*pow(F(2,1),2) + 2*d[11]*F(0,0)*F(0,2)*pow(F(2,1),2) + 2*d[10]*F(0,1)*F(0,2)*pow(F(2,1),2) + d[8]*pow(F(0,2),2)*pow(F(2,1),2) + 2*d[30]*pow(F(0,0),2)*F(2,0)*F(2,2) + 4*d[33]*F(0,0)*F(0,1)*F(2,0)*F(2,2) + 2*d[31]*pow(F(0,1),2)*F(2,0)*F(2,2) + 4*d[35]*F(0,0)*F(0,2)*F(2,0)*F(2,2) + 4*d[34]*F(0,1)*F(0,2)*F(2,0)*F(2,2) + 2*d[32]*pow(F(0,2),2)*F(2,0)*F(2,2) + 2*d[24]*pow(F(0,0),2)*F(2,1)*F(2,2) + 4*d[27]*F(0,0)*F(0,1)*F(2,1)*F(2,2) + 2*d[25]*pow(F(0,1),2)*F(2,1)*F(2,2) + 4*d[29]*F(0,0)*F(0,2)*F(2,1)*F(2,2) + 4*d[28]*F(0,1)*F(0,2)*F(2,1)*F(2,2) + 2*d[26]*pow(F(0,2),2)*F(2,1)*F(2,2) + d[12]*pow(F(0,0),2)*pow(F(2,2),2) + 2*d[15]*F(0,0)*F(0,1)*pow(F(2,2),2) + d[13]*pow(F(0,1),2)*pow(F(2,2),2) + 2*d[17]*F(0,0)*F(0,2)*pow(F(2,2),2) + 2*d[16]*F(0,1)*F(0,2)*pow(F(2,2),2) + d[14]*pow(F(0,2),2)*pow(F(2,2),2);

    c.d[13] = d[0]*pow(F(1,0),2)*pow(F(2,0),2) + 2*d[3]*F(1,0)*F(1,1)*pow(F(2,0),2) + d[1]*pow(F(1,1),2)*pow(F(2,0),2) + 2*d[5]*F(1,0)*F(1,2)*pow(F(2,0),2) + 2*d[4]*F(1,1)*F(1,2)*pow(F(2,0),2) + d[2]*pow(F(1,2),2)*pow(F(2,0),2) + 2*d[18]*pow(F(1,0),2)*F(2,0)*F(2,1) + 4*d[21]*F(1,0)*F(1,1)*F(2,0)*F(2,1) + 2*d[19]*pow(F(1,1),2)*F(2,0)*F(2,1) + 4*d[23]*F(1,0)*F(1,2)*F(2,0)*F(2,1) + 4*d[22]*F(1,1)*F(1,2)*F(2,0)*F(2,1) + 2*d[20]*pow(F(1,2),2)*F(2,0)*F(2,1) + d[6]*pow(F(1,0),2)*pow(F(2,1),2) + 2*d[9]*F(1,0)*F(1,1)*pow(F(2,1),2) + d[7]*pow(F(1,1),2)*pow(F(2,1),2) + 2*d[11]*F(1,0)*F(1,2)*pow(F(2,1),2) + 2*d[10]*F(1,1)*F(1,2)*pow(F(2,1),2) + d[8]*pow(F(1,2),2)*pow(F(2,1),2) + 2*d[30]*pow(F(1,0),2)*F(2,0)*F(2,2) + 4*d[33]*F(1,0)*F(1,1)*F(2,0)*F(2,2) + 2*d[31]*pow(F(1,1),2)*F(2,0)*F(2,2) + 4*d[35]*F(1,0)*F(1,2)*F(2,0)*F(2,2) + 4*d[34]*F(1,1)*F(1,2)*F(2,0)*F(2,2) + 2*d[32]*pow(F(1,2),2)*F(2,0)*F(2,2) + 2*d[24]*pow(F(1,0),2)*F(2,1)*F(2,2) + 4*d[27]*F(1,0)*F(1,1)*F(2,1)*F(2,2) + 2*d[25]*pow(F(1,1),2)*F(2,1)*F(2,2) + 4*d[29]*F(1,0)*F(1,2)*F(2,1)*F(2,2) + 4*d[28]*F(1,1)*F(1,2)*F(2,1)*F(2,2) + 2*d[26]*pow(F(1,2),2)*F(2,1)*F(2,2) + d[12]*pow(F(1,0),2)*pow(F(2,2),2) + 2*d[15]*F(1,0)*F(1,1)*pow(F(2,2),2) + d[13]*pow(F(1,1),2)*pow(F(2,2),2) + 2*d[17]*F(1,0)*F(1,2)*pow(F(2,2),2) + 2*d[16]*F(1,1)*F(1,2)*pow(F(2,2),2) + d[14]*pow(F(1,2),2)*pow(F(2,2),2);

    c.d[14] = d[0]*pow(F(2,0),4) + 2*d[3]*pow(F(2,0),3)*F(2,1) + 2*d[18]*pow(F(2,0),3)*F(2,1) + d[1]*pow(F(2,0),2)*pow(F(2,1),2) + d[6]*pow(F(2,0),2)*pow(F(2,1),2) + 4*d[21]*pow(F(2,0),2)*pow(F(2,1),2) + 2*d[9]*F(2,0)*pow(F(2,1),3) + 2*d[19]*F(2,0)*pow(F(2,1),3) + d[7]*pow(F(2,1),4) + 2*d[5]*pow(F(2,0),3)*F(2,2) + 2*d[30]*pow(F(2,0),3)*F(2,2) + 2*d[4]*pow(F(2,0),2)*F(2,1)*F(2,2) + 4*d[23]*pow(F(2,0),2)*F(2,1)*F(2,2) + 2*d[24]*pow(F(2,0),2)*F(2,1)*F(2,2) + 4*d[33]*pow(F(2,0),2)*F(2,1)*F(2,2) + 2*d[11]*F(2,0)*pow(F(2,1),2)*F(2,2) + 4*d[22]*F(2,0)*pow(F(2,1),2)*F(2,2) + 4*d[27]*F(2,0)*pow(F(2,1),2)*F(2,2) + 2*d[31]*F(2,0)*pow(F(2,1),2)*F(2,2) + 2*d[10]*pow(F(2,1),3)*F(2,2) + 2*d[25]*pow(F(2,1),3)*F(2,2) + d[2]*pow(F(2,0),2)*pow(F(2,2),2) + d[12]*pow(F(2,0),2)*pow(F(2,2),2) + 4*d[35]*pow(F(2,0),2)*pow(F(2,2),2) + 2*d[15]*F(2,0)*F(2,1)*pow(F(2,2),2) + 2*d[20]*F(2,0)*F(2,1)*pow(F(2,2),2) + 4*d[29]*F(2,0)*F(2,1)*pow(F(2,2),2) + 4*d[34]*F(2,0)*F(2,1)*pow(F(2,2),2) + d[8]*pow(F(2,1),2)*pow(F(2,2),2) + d[13]*pow(F(2,1),2)*pow(F(2,2),2) + 4*d[28]*pow(F(2,1),2)*pow(F(2,2),2) + 2*d[17]*F(2,0)*pow(F(2,2),3) + 2*d[32]*F(2,0)*pow(F(2,2),3) + 2*d[16]*F(2,1)*pow(F(2,2),3) + 2*d[26]*F(2,1)*pow(F(2,2),3) + d[14]*pow(F(2,2),4);

    c.d[15] = d[0]*F(0,0)*F(1,0)*pow(F(2,0),2) + d[5]*F(0,2)*F(1,0)*pow(F(2,0),2) + d[1]*F(0,1)*F(1,1)*pow(F(2,0),2) + d[4]*F(0,2)*F(1,1)*pow(F(2,0),2) + d[3]*(F(0,1)*F(1,0) + F(0,0)*F(1,1))*pow(F(2,0),2) + d[5]*F(0,0)*F(1,2)*pow(F(2,0),2) + d[4]*F(0,1)*F(1,2)*pow(F(2,0),2) + d[2]*F(0,2)*F(1,2)*pow(F(2,0),2) + 2*d[18]*F(0,0)*F(1,0)*F(2,0)*F(2,1) + 2*d[21]*F(0,1)*F(1,0)*F(2,0)*F(2,1) + 2*d[23]*F(0,2)*F(1,0)*F(2,0)*F(2,1) + 2*d[21]*F(0,0)*F(1,1)*F(2,0)*F(2,1) + 2*d[19]*F(0,1)*F(1,1)*F(2,0)*F(2,1) + 2*d[22]*F(0,2)*F(1,1)*F(2,0)*F(2,1) + 2*d[23]*F(0,0)*F(1,2)*F(2,0)*F(2,1) + 2*d[22]*F(0,1)*F(1,2)*F(2,0)*F(2,1) + 2*d[20]*F(0,2)*F(1,2)*F(2,0)*F(2,1) + d[6]*F(0,0)*F(1,0)*pow(F(2,1),2) + d[9]*F(0,1)*F(1,0)*pow(F(2,1),2) + d[11]*F(0,2)*F(1,0)*pow(F(2,1),2) + d[9]*F(0,0)*F(1,1)*pow(F(2,1),2) + d[7]*F(0,1)*F(1,1)*pow(F(2,1),2) + d[10]*F(0,2)*F(1,1)*pow(F(2,1),2) + d[11]*F(0,0)*F(1,2)*pow(F(2,1),2) + d[10]*F(0,1)*F(1,2)*pow(F(2,1),2) + d[8]*F(0,2)*F(1,2)*pow(F(2,1),2) + 2*d[30]*F(0,0)*F(1,0)*F(2,0)*F(2,2) + 2*d[33]*F(0,1)*F(1,0)*F(2,0)*F(2,2) + 2*d[35]*F(0,2)*F(1,0)*F(2,0)*F(2,2) + 2*d[33]*F(0,0)*F(1,1)*F(2,0)*F(2,2) + 2*d[31]*F(0,1)*F(1,1)*F(2,0)*F(2,2) + 2*d[34]*F(0,2)*F(1,1)*F(2,0)*F(2,2) + 2*d[35]*F(0,0)*F(1,2)*F(2,0)*F(2,2) + 2*d[34]*F(0,1)*F(1,2)*F(2,0)*F(2,2) + 2*d[32]*F(0,2)*F(1,2)*F(2,0)*F(2,2) + 2*d[24]*F(0,0)*F(1,0)*F(2,1)*F(2,2) + 2*d[27]*F(0,1)*F(1,0)*F(2,1)*F(2,2) + 2*d[29]*F(0,2)*F(1,0)*F(2,1)*F(2,2) + 2*d[27]*F(0,0)*F(1,1)*F(2,1)*F(2,2) + 2*d[25]*F(0,1)*F(1,1)*F(2,1)*F(2,2) + 2*d[28]*F(0,2)*F(1,1)*F(2,1)*F(2,2) + 2*d[29]*F(0,0)*F(1,2)*F(2,1)*F(2,2) + 2*d[28]*F(0,1)*F(1,2)*F(2,1)*F(2,2) + 2*d[26]*F(0,2)*F(1,2)*F(2,1)*F(2,2) + d[12]*F(0,0)*F(1,0)*pow(F(2,2),2) + d[15]*F(0,1)*F(1,0)*pow(F(2,2),2) + d[17]*F(0,2)*F(1,0)*pow(F(2,2),2) + d[15]*F(0,0)*F(1,1)*pow(F(2,2),2) + d[13]*F(0,1)*F(1,1)*pow(F(2,2),2) + d[16]*F(0,2)*F(1,1)*pow(F(2,2),2) + d[17]*F(0,0)*F(1,2)*pow(F(2,2),2) + d[16]*F(0,1)*F(1,2)*pow(F(2,2),2) + d[14]*F(0,2)*F(1,2)*pow(F(2,2),2);

    c.d[16] = d[0]*F(1,0)*pow(F(2,0),3) + d[5]*F(1,2)*pow(F(2,0),3) + 2*d[18]*F(1,0)*pow(F(2,0),2)*F(2,1) + d[1]*F(1,1)*pow(F(2,0),2)*F(2,1) + 2*d[21]*F(1,1)*pow(F(2,0),2)*F(2,1) + d[4]*F(1,2)*pow(F(2,0),2)*F(2,1) + 2*d[23]*F(1,2)*pow(F(2,0),2)*F(2,1) + d[6]*F(1,0)*F(2,0)*pow(F(2,1),2) + 2*d[21]*F(1,0)*F(2,0)*pow(F(2,1),2) + d[9]*F(1,1)*F(2,0)*pow(F(2,1),2) + 2*d[19]*F(1,1)*F(2,0)*pow(F(2,1),2) + d[11]*F(1,2)*F(2,0)*pow(F(2,1),2) + 2*d[22]*F(1,2)*F(2,0)*pow(F(2,1),2) + d[9]*F(1,0)*pow(F(2,1),3) + d[7]*F(1,1)*pow(F(2,1),3) + d[10]*F(1,2)*pow(F(2,1),3) + d[3]*pow(F(2,0),2)*(F(1,1)*F(2,0) + F(1,0)*F(2,1)) + d[5]*F(1,0)*pow(F(2,0),2)*F(2,2) + 2*d[30]*F(1,0)*pow(F(2,0),2)*F(2,2) + d[4]*F(1,1)*pow(F(2,0),2)*F(2,2) + 2*d[33]*F(1,1)*pow(F(2,0),2)*F(2,2) + d[2]*F(1,2)*pow(F(2,0),2)*F(2,2) + 2*d[35]*F(1,2)*pow(F(2,0),2)*F(2,2) + 2*d[23]*F(1,0)*F(2,0)*F(2,1)*F(2,2) + 2*d[24]*F(1,0)*F(2,0)*F(2,1)*F(2,2) + 2*d[33]*F(1,0)*F(2,0)*F(2,1)*F(2,2) + 2*d[22]*F(1,1)*F(2,0)*F(2,1)*F(2,2) + 2*d[27]*F(1,1)*F(2,0)*F(2,1)*F(2,2) + 2*d[31]*F(1,1)*F(2,0)*F(2,1)*F(2,2) + 2*d[20]*F(1,2)*F(2,0)*F(2,1)*F(2,2) + 2*d[29]*F(1,2)*F(2,0)*F(2,1)*F(2,2) + 2*d[34]*F(1,2)*F(2,0)*F(2,1)*F(2,2) + d[11]*F(1,0)*pow(F(2,1),2)*F(2,2) + 2*d[27]*F(1,0)*pow(F(2,1),2)*F(2,2) + d[10]*F(1,1)*pow(F(2,1),2)*F(2,2) + 2*d[25]*F(1,1)*pow(F(2,1),2)*F(2,2) + d[8]*F(1,2)*pow(F(2,1),2)*F(2,2) + 2*d[28]*F(1,2)*pow(F(2,1),2)*F(2,2) + d[12]*F(1,0)*F(2,0)*pow(F(2,2),2) + 2*d[35]*F(1,0)*F(2,0)*pow(F(2,2),2) + d[15]*F(1,1)*F(2,0)*pow(F(2,2),2) + 2*d[34]*F(1,1)*F(2,0)*pow(F(2,2),2) + d[17]*F(1,2)*F(2,0)*pow(F(2,2),2) + 2*d[32]*F(1,2)*F(2,0)*pow(F(2,2),2) + d[15]*F(1,0)*F(2,1)*pow(F(2,2),2) + 2*d[29]*F(1,0)*F(2,1)*pow(F(2,2),2) + d[13]*F(1,1)*F(2,1)*pow(F(2,2),2) + 2*d[28]*F(1,1)*F(2,1)*pow(F(2,2),2) + d[16]*F(1,2)*F(2,1)*pow(F(2,2),2) + 2*d[26]*F(1,2)*F(2,1)*pow(F(2,2),2) + d[17]*F(1,0)*pow(F(2,2),3) + d[16]*F(1,1)*pow(F(2,2),3) + d[14]*F(1,2)*pow(F(2,2),3);

    c.d[17] = d[0]*F(0,0)*pow(F(2,0),3) + d[5]*F(0,2)*pow(F(2,0),3) + 2*d[18]*F(0,0)*pow(F(2,0),2)*F(2,1) + d[1]*F(0,1)*pow(F(2,0),2)*F(2,1) + 2*d[21]*F(0,1)*pow(F(2,0),2)*F(2,1) + d[4]*F(0,2)*pow(F(2,0),2)*F(2,1) + 2*d[23]*F(0,2)*pow(F(2,0),2)*F(2,1) + d[6]*F(0,0)*F(2,0)*pow(F(2,1),2) + 2*d[21]*F(0,0)*F(2,0)*pow(F(2,1),2) + d[9]*F(0,1)*F(2,0)*pow(F(2,1),2) + 2*d[19]*F(0,1)*F(2,0)*pow(F(2,1),2) + d[11]*F(0,2)*F(2,0)*pow(F(2,1),2) + 2*d[22]*F(0,2)*F(2,0)*pow(F(2,1),2) + d[9]*F(0,0)*pow(F(2,1),3) + d[7]*F(0,1)*pow(F(2,1),3) + d[10]*F(0,2)*pow(F(2,1),3) + d[3]*pow(F(2,0),2)*(F(0,1)*F(2,0) + F(0,0)*F(2,1)) + d[5]*F(0,0)*pow(F(2,0),2)*F(2,2) + 2*d[30]*F(0,0)*pow(F(2,0),2)*F(2,2) + d[4]*F(0,1)*pow(F(2,0),2)*F(2,2) + 2*d[33]*F(0,1)*pow(F(2,0),2)*F(2,2) + d[2]*F(0,2)*pow(F(2,0),2)*F(2,2) + 2*d[35]*F(0,2)*pow(F(2,0),2)*F(2,2) + 2*d[23]*F(0,0)*F(2,0)*F(2,1)*F(2,2) + 2*d[24]*F(0,0)*F(2,0)*F(2,1)*F(2,2) + 2*d[33]*F(0,0)*F(2,0)*F(2,1)*F(2,2) + 2*d[22]*F(0,1)*F(2,0)*F(2,1)*F(2,2) + 2*d[27]*F(0,1)*F(2,0)*F(2,1)*F(2,2) + 2*d[31]*F(0,1)*F(2,0)*F(2,1)*F(2,2) + 2*d[20]*F(0,2)*F(2,0)*F(2,1)*F(2,2) + 2*d[29]*F(0,2)*F(2,0)*F(2,1)*F(2,2) + 2*d[34]*F(0,2)*F(2,0)*F(2,1)*F(2,2) + d[11]*F(0,0)*pow(F(2,1),2)*F(2,2) + 2*d[27]*F(0,0)*pow(F(2,1),2)*F(2,2) + d[10]*F(0,1)*pow(F(2,1),2)*F(2,2) + 2*d[25]*F(0,1)*pow(F(2,1),2)*F(2,2) + d[8]*F(0,2)*pow(F(2,1),2)*F(2,2) + 2*d[28]*F(0,2)*pow(F(2,1),2)*F(2,2) + d[12]*F(0,0)*F(2,0)*pow(F(2,2),2) + 2*d[35]*F(0,0)*F(2,0)*pow(F(2,2),2) + d[15]*F(0,1)*F(2,0)*pow(F(2,2),2) + 2*d[34]*F(0,1)*F(2,0)*pow(F(2,2),2) + d[17]*F(0,2)*F(2,0)*pow(F(2,2),2) + 2*d[32]*F(0,2)*F(2,0)*pow(F(2,2),2) + d[15]*F(0,0)*F(2,1)*pow(F(2,2),2) + 2*d[29]*F(0,0)*F(2,1)*pow(F(2,2),2) + d[13]*F(0,1)*F(2,1)*pow(F(2,2),2) + 2*d[28]*F(0,1)*F(2,1)*pow(F(2,2),2) + d[16]*F(0,2)*F(2,1)*pow(F(2,2),2) + 2*d[26]*F(0,2)*F(2,1)*pow(F(2,2),2) + d[17]*F(0,0)*pow(F(2,2),3) + d[16]*F(0,1)*pow(F(2,2),3) + d[14]*F(0,2)*pow(F(2,2),3);

    c.d[18] = d[0]*pow(F(0,0),3)*F(1,0) + 2*d[3]*pow(F(0,0),2)*F(0,1)*F(1,0) + d[18]*pow(F(0,0),2)*F(0,1)*F(1,0) + d[1]*F(0,0)*pow(F(0,1),2)*F(1,0) + 2*d[21]*F(0,0)*pow(F(0,1),2)*F(1,0) + d[19]*pow(F(0,1),3)*F(1,0) + 2*d[5]*pow(F(0,0),2)*F(0,2)*F(1,0) + d[30]*pow(F(0,0),2)*F(0,2)*F(1,0) + 2*d[4]*F(0,0)*F(0,1)*F(0,2)*F(1,0) + 2*d[23]*F(0,0)*F(0,1)*F(0,2)*F(1,0) + 2*d[33]*F(0,0)*F(0,1)*F(0,2)*F(1,0) + 2*d[22]*pow(F(0,1),2)*F(0,2)*F(1,0) + d[31]*pow(F(0,1),2)*F(0,2)*F(1,0) + d[2]*F(0,0)*pow(F(0,2),2)*F(1,0) + 2*d[35]*F(0,0)*pow(F(0,2),2)*F(1,0) + d[20]*F(0,1)*pow(F(0,2),2)*F(1,0) + 2*d[34]*F(0,1)*pow(F(0,2),2)*F(1,0) + d[32]*pow(F(0,2),3)*F(1,0) + d[18]*pow(F(0,0),3)*F(1,1) + d[6]*pow(F(0,0),2)*F(0,1)*F(1,1) + 2*d[21]*pow(F(0,0),2)*F(0,1)*F(1,1) + 2*d[9]*F(0,0)*pow(F(0,1),2)*F(1,1) + d[19]*F(0,0)*pow(F(0,1),2)*F(1,1) + d[7]*pow(F(0,1),3)*F(1,1) + 2*d[23]*pow(F(0,0),2)*F(0,2)*F(1,1) + d[24]*pow(F(0,0),2)*F(0,2)*F(1,1) + 2*d[11]*F(0,0)*F(0,1)*F(0,2)*F(1,1) + 2*d[22]*F(0,0)*F(0,1)*F(0,2)*F(1,1) + 2*d[27]*F(0,0)*F(0,1)*F(0,2)*F(1,1) + 2*d[10]*pow(F(0,1),2)*F(0,2)*F(1,1) + d[25]*pow(F(0,1),2)*F(0,2)*F(1,1) + d[20]*F(0,0)*pow(F(0,2),2)*F(1,1) + 2*d[29]*F(0,0)*pow(F(0,2),2)*F(1,1) + d[8]*F(0,1)*pow(F(0,2),2)*F(1,1) + 2*d[28]*F(0,1)*pow(F(0,2),2)*F(1,1) + d[26]*pow(F(0,2),3)*F(1,1) + d[30]*pow(F(0,0),3)*F(1,2) + d[24]*pow(F(0,0),2)*F(0,1)*F(1,2) + 2*d[33]*pow(F(0,0),2)*F(0,1)*F(1,2) + 2*d[27]*F(0,0)*pow(F(0,1),2)*F(1,2) + d[31]*F(0,0)*pow(F(0,1),2)*F(1,2) + d[25]*pow(F(0,1),3)*F(1,2) + d[12]*pow(F(0,0),2)*F(0,2)*F(1,2) + 2*d[35]*pow(F(0,0),2)*F(0,2)*F(1,2) + 2*d[15]*F(0,0)*F(0,1)*F(0,2)*F(1,2) + 2*d[29]*F(0,0)*F(0,1)*F(0,2)*F(1,2) + 2*d[34]*F(0,0)*F(0,1)*F(0,2)*F(1,2) + d[13]*pow(F(0,1),2)*F(0,2)*F(1,2) + 2*d[28]*pow(F(0,1),2)*F(0,2)*F(1,2) + 2*d[17]*F(0,0)*pow(F(0,2),2)*F(1,2) + d[32]*F(0,0)*pow(F(0,2),2)*F(1,2) + 2*d[16]*F(0,1)*pow(F(0,2),2)*F(1,2) + d[26]*F(0,1)*pow(F(0,2),2)*F(1,2) + d[14]*pow(F(0,2),3)*F(1,2);

    c.d[19] = d[0]*F(0,0)*pow(F(1,0),3) + d[30]*F(0,2)*pow(F(1,0),3) + 2*d[3]*F(0,0)*pow(F(1,0),2)*F(1,1) + d[6]*F(0,1)*pow(F(1,0),2)*F(1,1) + 2*d[21]*F(0,1)*pow(F(1,0),2)*F(1,1) + d[24]*F(0,2)*pow(F(1,0),2)*F(1,1) + 2*d[33]*F(0,2)*pow(F(1,0),2)*F(1,1) + d[1]*F(0,0)*F(1,0)*pow(F(1,1),2) + 2*d[21]*F(0,0)*F(1,0)*pow(F(1,1),2) + 2*d[9]*F(0,1)*F(1,0)*pow(F(1,1),2) + d[19]*F(0,1)*F(1,0)*pow(F(1,1),2) + 2*d[27]*F(0,2)*F(1,0)*pow(F(1,1),2) + d[31]*F(0,2)*F(1,0)*pow(F(1,1),2) + d[19]*F(0,0)*pow(F(1,1),3) + d[7]*F(0,1)*pow(F(1,1),3) + d[25]*F(0,2)*pow(F(1,1),3) + d[18]*pow(F(1,0),2)*(F(0,1)*F(1,0) + F(0,0)*F(1,1)) + 2*d[5]*F(0,0)*pow(F(1,0),2)*F(1,2) + d[30]*F(0,0)*pow(F(1,0),2)*F(1,2) + 2*d[23]*F(0,1)*pow(F(1,0),2)*F(1,2) + d[24]*F(0,1)*pow(F(1,0),2)*F(1,2) + d[12]*F(0,2)*pow(F(1,0),2)*F(1,2) + 2*d[35]*F(0,2)*pow(F(1,0),2)*F(1,2) + 2*d[4]*F(0,0)*F(1,0)*F(1,1)*F(1,2) + 2*d[23]*F(0,0)*F(1,0)*F(1,1)*F(1,2) + 2*d[33]*F(0,0)*F(1,0)*F(1,1)*F(1,2) + 2*d[11]*F(0,1)*F(1,0)*F(1,1)*F(1,2) + 2*d[22]*F(0,1)*F(1,0)*F(1,1)*F(1,2) + 2*d[27]*F(0,1)*F(1,0)*F(1,1)*F(1,2) + 2*d[15]*F(0,2)*F(1,0)*F(1,1)*F(1,2) + 2*d[29]*F(0,2)*F(1,0)*F(1,1)*F(1,2) + 2*d[34]*F(0,2)*F(1,0)*F(1,1)*F(1,2) + 2*d[22]*F(0,0)*pow(F(1,1),2)*F(1,2) + d[31]*F(0,0)*pow(F(1,1),2)*F(1,2) + 2*d[10]*F(0,1)*pow(F(1,1),2)*F(1,2) + d[25]*F(0,1)*pow(F(1,1),2)*F(1,2) + d[13]*F(0,2)*pow(F(1,1),2)*F(1,2) + 2*d[28]*F(0,2)*pow(F(1,1),2)*F(1,2) + d[2]*F(0,0)*F(1,0)*pow(F(1,2),2) + 2*d[35]*F(0,0)*F(1,0)*pow(F(1,2),2) + d[20]*F(0,1)*F(1,0)*pow(F(1,2),2) + 2*d[29]*F(0,1)*F(1,0)*pow(F(1,2),2) + 2*d[17]*F(0,2)*F(1,0)*pow(F(1,2),2) + d[32]*F(0,2)*F(1,0)*pow(F(1,2),2) + d[20]*F(0,0)*F(1,1)*pow(F(1,2),2) + 2*d[34]*F(0,0)*F(1,1)*pow(F(1,2),2) + d[8]*F(0,1)*F(1,1)*pow(F(1,2),2) + 2*d[28]*F(0,1)*F(1,1)*pow(F(1,2),2) + 2*d[16]*F(0,2)*F(1,1)*pow(F(1,2),2) + d[26]*F(0,2)*F(1,1)*pow(F(1,2),2) + d[32]*F(0,0)*pow(F(1,2),3) + d[26]*F(0,1)*pow(F(1,2),3) + d[14]*F(0,2)*pow(F(1,2),3);

    c.d[20] = d[0]*F(0,0)*F(1,0)*pow(F(2,0),2) + d[30]*F(0,2)*F(1,0)*pow(F(2,0),2) + d[6]*F(0,1)*F(1,1)*pow(F(2,0),2) + d[24]*F(0,2)*F(1,1)*pow(F(2,0),2) + d[18]*(F(0,1)*F(1,0) + F(0,0)*F(1,1))*pow(F(2,0),2) + d[30]*F(0,0)*F(1,2)*pow(F(2,0),2) + d[24]*F(0,1)*F(1,2)*pow(F(2,0),2) + d[12]*F(0,2)*F(1,2)*pow(F(2,0),2) + 2*d[3]*F(0,0)*F(1,0)*F(2,0)*F(2,1) + 2*d[21]*F(0,1)*F(1,0)*F(2,0)*F(2,1) + 2*d[33]*F(0,2)*F(1,0)*F(2,0)*F(2,1) + 2*d[21]*F(0,0)*F(1,1)*F(2,0)*F(2,1) + 2*d[9]*F(0,1)*F(1,1)*F(2,0)*F(2,1) + 2*d[27]*F(0,2)*F(1,1)*F(2,0)*F(2,1) + 2*d[33]*F(0,0)*F(1,2)*F(2,0)*F(2,1) + 2*d[27]*F(0,1)*F(1,2)*F(2,0)*F(2,1) + 2*d[15]*F(0,2)*F(1,2)*F(2,0)*F(2,1) + d[1]*F(0,0)*F(1,0)*pow(F(2,1),2) + d[19]*F(0,1)*F(1,0)*pow(F(2,1),2) + d[31]*F(0,2)*F(1,0)*pow(F(2,1),2) + d[19]*F(0,0)*F(1,1)*pow(F(2,1),2) + d[7]*F(0,1)*F(1,1)*pow(F(2,1),2) + d[25]*F(0,2)*F(1,1)*pow(F(2,1),2) + d[31]*F(0,0)*F(1,2)*pow(F(2,1),2) + d[25]*F(0,1)*F(1,2)*pow(F(2,1),2) + d[13]*F(0,2)*F(1,2)*pow(F(2,1),2) + 2*d[5]*F(0,0)*F(1,0)*F(2,0)*F(2,2) + 2*d[23]*F(0,1)*F(1,0)*F(2,0)*F(2,2) + 2*d[35]*F(0,2)*F(1,0)*F(2,0)*F(2,2) + 2*d[23]*F(0,0)*F(1,1)*F(2,0)*F(2,2) + 2*d[11]*F(0,1)*F(1,1)*F(2,0)*F(2,2) + 2*d[29]*F(0,2)*F(1,1)*F(2,0)*F(2,2) + 2*d[35]*F(0,0)*F(1,2)*F(2,0)*F(2,2) + 2*d[29]*F(0,1)*F(1,2)*F(2,0)*F(2,2) + 2*d[17]*F(0,2)*F(1,2)*F(2,0)*F(2,2) + 2*d[4]*F(0,0)*F(1,0)*F(2,1)*F(2,2) + 2*d[22]*F(0,1)*F(1,0)*F(2,1)*F(2,2) + 2*d[34]*F(0,2)*F(1,0)*F(2,1)*F(2,2) + 2*d[22]*F(0,0)*F(1,1)*F(2,1)*F(2,2) + 2*d[10]*F(0,1)*F(1,1)*F(2,1)*F(2,2) + 2*d[28]*F(0,2)*F(1,1)*F(2,1)*F(2,2) + 2*d[34]*F(0,0)*F(1,2)*F(2,1)*F(2,2) + 2*d[28]*F(0,1)*F(1,2)*F(2,1)*F(2,2) + 2*d[16]*F(0,2)*F(1,2)*F(2,1)*F(2,2) + d[2]*F(0,0)*F(1,0)*pow(F(2,2),2) + d[20]*F(0,1)*F(1,0)*pow(F(2,2),2) + d[32]*F(0,2)*F(1,0)*pow(F(2,2),2) + d[20]*F(0,0)*F(1,1)*pow(F(2,2),2) + d[8]*F(0,1)*F(1,1)*pow(F(2,2),2) + d[26]*F(0,2)*F(1,1)*pow(F(2,2),2) + d[32]*F(0,0)*F(1,2)*pow(F(2,2),2) + d[26]*F(0,1)*F(1,2)*pow(F(2,2),2) + d[14]*F(0,2)*F(1,2)*pow(F(2,2),2);

    c.d[21] = d[0]*pow(F(0,0),2)*pow(F(1,0),2) + d[18]*F(0,0)*F(0,1)*pow(F(1,0),2) + d[21]*pow(F(0,1),2)*pow(F(1,0),2) + d[5]*F(0,0)*F(0,2)*pow(F(1,0),2) + d[30]*F(0,0)*F(0,2)*pow(F(1,0),2) + d[23]*F(0,1)*F(0,2)*pow(F(1,0),2) + d[33]*F(0,1)*F(0,2)*pow(F(1,0),2) + d[35]*pow(F(0,2),2)*pow(F(1,0),2) + d[18]*pow(F(0,0),2)*F(1,0)*F(1,1) + d[1]*F(0,0)*F(0,1)*F(1,0)*F(1,1) + d[6]*F(0,0)*F(0,1)*F(1,0)*F(1,1) + 2*d[21]*F(0,0)*F(0,1)*F(1,0)*F(1,1) + d[9]*pow(F(0,1),2)*F(1,0)*F(1,1) + d[19]*pow(F(0,1),2)*F(1,0)*F(1,1) + d[4]*F(0,0)*F(0,2)*F(1,0)*F(1,1) + d[23]*F(0,0)*F(0,2)*F(1,0)*F(1,1) + d[24]*F(0,0)*F(0,2)*F(1,0)*F(1,1) + d[33]*F(0,0)*F(0,2)*F(1,0)*F(1,1) + d[11]*F(0,1)*F(0,2)*F(1,0)*F(1,1) + d[22]*F(0,1)*F(0,2)*F(1,0)*F(1,1) + d[27]*F(0,1)*F(0,2)*F(1,0)*F(1,1) + d[31]*F(0,1)*F(0,2)*F(1,0)*F(1,1) + d[29]*pow(F(0,2),2)*F(1,0)*F(1,1) + d[34]*pow(F(0,2),2)*F(1,0)*F(1,1) + d[21]*pow(F(0,0),2)*pow(F(1,1),2) + d[9]*F(0,0)*F(0,1)*pow(F(1,1),2) + d[19]*F(0,0)*F(0,1)*pow(F(1,1),2) + d[7]*pow(F(0,1),2)*pow(F(1,1),2) + d[22]*F(0,0)*F(0,2)*pow(F(1,1),2) + d[27]*F(0,0)*F(0,2)*pow(F(1,1),2) + d[10]*F(0,1)*F(0,2)*pow(F(1,1),2) + d[25]*F(0,1)*F(0,2)*pow(F(1,1),2) + d[28]*pow(F(0,2),2)*pow(F(1,1),2) + d[3]*F(0,0)*F(1,0)*(F(0,1)*F(1,0) + F(0,0)*F(1,1)) + d[5]*pow(F(0,0),2)*F(1,0)*F(1,2) + d[30]*pow(F(0,0),2)*F(1,0)*F(1,2) + d[4]*F(0,0)*F(0,1)*F(1,0)*F(1,2) + d[23]*F(0,0)*F(0,1)*F(1,0)*F(1,2) + d[24]*F(0,0)*F(0,1)*F(1,0)*F(1,2) + d[33]*F(0,0)*F(0,1)*F(1,0)*F(1,2) + d[22]*pow(F(0,1),2)*F(1,0)*F(1,2) + d[27]*pow(F(0,1),2)*F(1,0)*F(1,2) + d[2]*F(0,0)*F(0,2)*F(1,0)*F(1,2) + d[12]*F(0,0)*F(0,2)*F(1,0)*F(1,2) + 2*d[35]*F(0,0)*F(0,2)*F(1,0)*F(1,2) + d[15]*F(0,1)*F(0,2)*F(1,0)*F(1,2) + d[20]*F(0,1)*F(0,2)*F(1,0)*F(1,2) + d[29]*F(0,1)*F(0,2)*F(1,0)*F(1,2) + d[34]*F(0,1)*F(0,2)*F(1,0)*F(1,2) + d[17]*pow(F(0,2),2)*F(1,0)*F(1,2) + d[32]*pow(F(0,2),2)*F(1,0)*F(1,2) + d[23]*pow(F(0,0),2)*F(1,1)*F(1,2) + d[33]*pow(F(0,0),2)*F(1,1)*F(1,2) + d[11]*F(0,0)*F(0,1)*F(1,1)*F(1,2) + d[22]*F(0,0)*F(0,1)*F(1,1)*F(1,2) + d[27]*F(0,0)*F(0,1)*F(1,1)*F(1,2) + d[31]*F(0,0)*F(0,1)*F(1,1)*F(1,2) + d[10]*pow(F(0,1),2)*F(1,1)*F(1,2) + d[25]*pow(F(0,1),2)*F(1,1)*F(1,2) + d[15]*F(0,0)*F(0,2)*F(1,1)*F(1,2) + d[20]*F(0,0)*F(0,2)*F(1,1)*F(1,2) + d[29]*F(0,0)*F(0,2)*F(1,1)*F(1,2) + d[34]*F(0,0)*F(0,2)*F(1,1)*F(1,2) + d[8]*F(0,1)*F(0,2)*F(1,1)*F(1,2) + d[13]*F(0,1)*F(0,2)*F(1,1)*F(1,2) + 2*d[28]*F(0,1)*F(0,2)*F(1,1)*F(1,2) + d[16]*pow(F(0,2),2)*F(1,1)*F(1,2) + d[26]*pow(F(0,2),2)*F(1,1)*F(1,2) + d[35]*pow(F(0,0),2)*pow(F(1,2),2) + d[29]*F(0,0)*F(0,1)*pow(F(1,2),2) + d[34]*F(0,0)*F(0,1)*pow(F(1,2),2) + d[28]*pow(F(0,1),2)*pow(F(1,2),2) + d[17]*F(0,0)*F(0,2)*pow(F(1,2),2) + d[32]*F(0,0)*F(0,2)*pow(F(1,2),2) + d[16]*F(0,1)*F(0,2)*pow(F(1,2),2) + d[26]*F(0,1)*F(0,2)*pow(F(1,2),2) + d[14]*pow(F(0,2),2)*pow(F(1,2),2);

    c.d[22] = d[0]*F(0,0)*pow(F(1,0),2)*F(2,0) + d[30]*F(0,2)*pow(F(1,0),2)*F(2,0) + d[3]*F(0,0)*F(1,0)*F(1,1)*F(2,0) + d[6]*F(0,1)*F(1,0)*F(1,1)*F(2,0) + d[21]*F(0,1)*F(1,0)*F(1,1)*F(2,0) + d[24]*F(0,2)*F(1,0)*F(1,1)*F(2,0) + d[33]*F(0,2)*F(1,0)*F(1,1)*F(2,0) + d[21]*F(0,0)*pow(F(1,1),2)*F(2,0) + d[9]*F(0,1)*pow(F(1,1),2)*F(2,0) + d[27]*F(0,2)*pow(F(1,1),2)*F(2,0) + d[18]*F(1,0)*(F(0,1)*F(1,0) + F(0,0)*F(1,1))*F(2,0) + d[5]*F(0,0)*F(1,0)*F(1,2)*F(2,0) + d[30]*F(0,0)*F(1,0)*F(1,2)*F(2,0) + d[23]*F(0,1)*F(1,0)*F(1,2)*F(2,0) + d[24]*F(0,1)*F(1,0)*F(1,2)*F(2,0) + d[12]*F(0,2)*F(1,0)*F(1,2)*F(2,0) + d[35]*F(0,2)*F(1,0)*F(1,2)*F(2,0) + d[23]*F(0,0)*F(1,1)*F(1,2)*F(2,0) + d[33]*F(0,0)*F(1,1)*F(1,2)*F(2,0) + d[11]*F(0,1)*F(1,1)*F(1,2)*F(2,0) + d[27]*F(0,1)*F(1,1)*F(1,2)*F(2,0) + d[15]*F(0,2)*F(1,1)*F(1,2)*F(2,0) + d[29]*F(0,2)*F(1,1)*F(1,2)*F(2,0) + d[35]*F(0,0)*pow(F(1,2),2)*F(2,0) + d[29]*F(0,1)*pow(F(1,2),2)*F(2,0) + d[17]*F(0,2)*pow(F(1,2),2)*F(2,0) + d[3]*F(0,0)*pow(F(1,0),2)*F(2,1) + d[21]*F(0,1)*pow(F(1,0),2)*F(2,1) + d[33]*F(0,2)*pow(F(1,0),2)*F(2,1) + d[1]*F(0,0)*F(1,0)*F(1,1)*F(2,1) + d[21]*F(0,0)*F(1,0)*F(1,1)*F(2,1) + d[9]*F(0,1)*F(1,0)*F(1,1)*F(2,1) + d[19]*F(0,1)*F(1,0)*F(1,1)*F(2,1) + d[27]*F(0,2)*F(1,0)*F(1,1)*F(2,1) + d[31]*F(0,2)*F(1,0)*F(1,1)*F(2,1) + d[19]*F(0,0)*pow(F(1,1),2)*F(2,1) + d[7]*F(0,1)*pow(F(1,1),2)*F(2,1) + d[25]*F(0,2)*pow(F(1,1),2)*F(2,1) + d[4]*F(0,0)*F(1,0)*F(1,2)*F(2,1) + d[33]*F(0,0)*F(1,0)*F(1,2)*F(2,1) + d[22]*F(0,1)*F(1,0)*F(1,2)*F(2,1) + d[27]*F(0,1)*F(1,0)*F(1,2)*F(2,1) + d[15]*F(0,2)*F(1,0)*F(1,2)*F(2,1) + d[34]*F(0,2)*F(1,0)*F(1,2)*F(2,1) + d[22]*F(0,0)*F(1,1)*F(1,2)*F(2,1) + d[31]*F(0,0)*F(1,1)*F(1,2)*F(2,1) + d[10]*F(0,1)*F(1,1)*F(1,2)*F(2,1) + d[25]*F(0,1)*F(1,1)*F(1,2)*F(2,1) + d[13]*F(0,2)*F(1,1)*F(1,2)*F(2,1) + d[28]*F(0,2)*F(1,1)*F(1,2)*F(2,1) + d[34]*F(0,0)*pow(F(1,2),2)*F(2,1) + d[28]*F(0,1)*pow(F(1,2),2)*F(2,1) + d[16]*F(0,2)*pow(F(1,2),2)*F(2,1) + d[5]*F(0,0)*pow(F(1,0),2)*F(2,2) + d[23]*F(0,1)*pow(F(1,0),2)*F(2,2) + d[35]*F(0,2)*pow(F(1,0),2)*F(2,2) + d[4]*F(0,0)*F(1,0)*F(1,1)*F(2,2) + d[23]*F(0,0)*F(1,0)*F(1,1)*F(2,2) + d[11]*F(0,1)*F(1,0)*F(1,1)*F(2,2) + d[22]*F(0,1)*F(1,0)*F(1,1)*F(2,2) + d[29]*F(0,2)*F(1,0)*F(1,1)*F(2,2) + d[34]*F(0,2)*F(1,0)*F(1,1)*F(2,2) + d[22]*F(0,0)*pow(F(1,1),2)*F(2,2) + d[10]*F(0,1)*pow(F(1,1),2)*F(2,2) + d[28]*F(0,2)*pow(F(1,1),2)*F(2,2) + d[2]*F(0,0)*F(1,0)*F(1,2)*F(2,2) + d[35]*F(0,0)*F(1,0)*F(1,2)*F(2,2) + d[20]*F(0,1)*F(1,0)*F(1,2)*F(2,2) + d[29]*F(0,1)*F(1,0)*F(1,2)*F(2,2) + d[17]*F(0,2)*F(1,0)*F(1,2)*F(2,2) + d[32]*F(0,2)*F(1,0)*F(1,2)*F(2,2) + d[20]*F(0,0)*F(1,1)*F(1,2)*F(2,2) + d[34]*F(0,0)*F(1,1)*F(1,2)*F(2,2) + d[8]*F(0,1)*F(1,1)*F(1,2)*F(2,2) + d[28]*F(0,1)*F(1,1)*F(1,2)*F(2,2) + d[16]*F(0,2)*F(1,1)*F(1,2)*F(2,2) + d[26]*F(0,2)*F(1,1)*F(1,2)*F(2,2) + d[32]*F(0,0)*pow(F(1,2),2)*F(2,2) + d[26]*F(0,1)*pow(F(1,2),2)*F(2,2) + d[14]*F(0,2)*pow(F(1,2),2)*F(2,2);

    c.d[23] = d[0]*pow(F(0,0),2)*F(1,0)*F(2,0) + d[18]*F(0,0)*F(0,1)*F(1,0)*F(2,0) + d[21]*pow(F(0,1),2)*F(1,0)*F(2,0) + d[5]*F(0,0)*F(0,2)*F(1,0)*F(2,0) + d[30]*F(0,0)*F(0,2)*F(1,0)*F(2,0) + d[23]*F(0,1)*F(0,2)*F(1,0)*F(2,0) + d[33]*F(0,1)*F(0,2)*F(1,0)*F(2,0) + d[35]*pow(F(0,2),2)*F(1,0)*F(2,0) + d[18]*pow(F(0,0),2)*F(1,1)*F(2,0) + d[6]*F(0,0)*F(0,1)*F(1,1)*F(2,0) + d[21]*F(0,0)*F(0,1)*F(1,1)*F(2,0) + d[9]*pow(F(0,1),2)*F(1,1)*F(2,0) + d[23]*F(0,0)*F(0,2)*F(1,1)*F(2,0) + d[24]*F(0,0)*F(0,2)*F(1,1)*F(2,0) + d[11]*F(0,1)*F(0,2)*F(1,1)*F(2,0) + d[27]*F(0,1)*F(0,2)*F(1,1)*F(2,0) + d[29]*pow(F(0,2),2)*F(1,1)*F(2,0) + d[30]*pow(F(0,0),2)*F(1,2)*F(2,0) + d[24]*F(0,0)*F(0,1)*F(1,2)*F(2,0) + d[33]*F(0,0)*F(0,1)*F(1,2)*F(2,0) + d[27]*pow(F(0,1),2)*F(1,2)*F(2,0) + d[12]*F(0,0)*F(0,2)*F(1,2)*F(2,0) + d[35]*F(0,0)*F(0,2)*F(1,2)*F(2,0) + d[15]*F(0,1)*F(0,2)*F(1,2)*F(2,0) + d[29]*F(0,1)*F(0,2)*F(1,2)*F(2,0) + d[17]*pow(F(0,2),2)*F(1,2)*F(2,0) + d[1]*F(0,0)*F(0,1)*F(1,0)*F(2,1) + d[21]*F(0,0)*F(0,1)*F(1,0)*F(2,1) + d[19]*pow(F(0,1),2)*F(1,0)*F(2,1) + d[4]*F(0,0)*F(0,2)*F(1,0)*F(2,1) + d[33]*F(0,0)*F(0,2)*F(1,0)*F(2,1) + d[22]*F(0,1)*F(0,2)*F(1,0)*F(2,1) + d[31]*F(0,1)*F(0,2)*F(1,0)*F(2,1) + d[34]*pow(F(0,2),2)*F(1,0)*F(2,1) + d[21]*pow(F(0,0),2)*F(1,1)*F(2,1) + d[9]*F(0,0)*F(0,1)*F(1,1)*F(2,1) + d[19]*F(0,0)*F(0,1)*F(1,1)*F(2,1) + d[7]*pow(F(0,1),2)*F(1,1)*F(2,1) + d[22]*F(0,0)*F(0,2)*F(1,1)*F(2,1) + d[27]*F(0,0)*F(0,2)*F(1,1)*F(2,1) + d[10]*F(0,1)*F(0,2)*F(1,1)*F(2,1) + d[25]*F(0,1)*F(0,2)*F(1,1)*F(2,1) + d[28]*pow(F(0,2),2)*F(1,1)*F(2,1) + d[33]*pow(F(0,0),2)*F(1,2)*F(2,1) + d[27]*F(0,0)*F(0,1)*F(1,2)*F(2,1) + d[31]*F(0,0)*F(0,1)*F(1,2)*F(2,1) + d[25]*pow(F(0,1),2)*F(1,2)*F(2,1) + d[15]*F(0,0)*F(0,2)*F(1,2)*F(2,1) + d[34]*F(0,0)*F(0,2)*F(1,2)*F(2,1) + d[13]*F(0,1)*F(0,2)*F(1,2)*F(2,1) + d[28]*F(0,1)*F(0,2)*F(1,2)*F(2,1) + d[16]*pow(F(0,2),2)*F(1,2)*F(2,1) + d[3]*F(0,0)*F(1,0)*(F(0,1)*F(2,0) + F(0,0)*F(2,1)) + d[5]*pow(F(0,0),2)*F(1,0)*F(2,2) + d[4]*F(0,0)*F(0,1)*F(1,0)*F(2,2) + d[23]*F(0,0)*F(0,1)*F(1,0)*F(2,2) + d[22]*pow(F(0,1),2)*F(1,0)*F(2,2) + d[2]*F(0,0)*F(0,2)*F(1,0)*F(2,2) + d[35]*F(0,0)*F(0,2)*F(1,0)*F(2,2) + d[20]*F(0,1)*F(0,2)*F(1,0)*F(2,2) + d[34]*F(0,1)*F(0,2)*F(1,0)*F(2,2) + d[32]*pow(F(0,2),2)*F(1,0)*F(2,2) + d[23]*pow(F(0,0),2)*F(1,1)*F(2,2) + d[11]*F(0,0)*F(0,1)*F(1,1)*F(2,2) + d[22]*F(0,0)*F(0,1)*F(1,1)*F(2,2) + d[10]*pow(F(0,1),2)*F(1,1)*F(2,2) + d[20]*F(0,0)*F(0,2)*F(1,1)*F(2,2) + d[29]*F(0,0)*F(0,2)*F(1,1)*F(2,2) + d[8]*F(0,1)*F(0,2)*F(1,1)*F(2,2) + d[28]*F(0,1)*F(0,2)*F(1,1)*F(2,2) + d[26]*pow(F(0,2),2)*F(1,1)*F(2,2) + d[35]*pow(F(0,0),2)*F(1,2)*F(2,2) + d[29]*F(0,0)*F(0,1)*F(1,2)*F(2,2) + d[34]*F(0,0)*F(0,1)*F(1,2)*F(2,2) + d[28]*pow(F(0,1),2)*F(1,2)*F(2,2) + d[17]*F(0,0)*F(0,2)*F(1,2)*F(2,2) + d[32]*F(0,0)*F(0,2)*F(1,2)*F(2,2) + d[16]*F(0,1)*F(0,2)*F(1,2)*F(2,2) + d[26]*F(0,1)*F(0,2)*F(1,2)*F(2,2) + d[14]*pow(F(0,2),2)*F(1,2)*F(2,2);

    c.d[24] = d[0]*pow(F(0,0),2)*F(1,0)*F(2,0) + 2*d[3]*F(0,0)*F(0,1)*F(1,0)*F(2,0) + d[1]*pow(F(0,1),2)*F(1,0)*F(2,0) + 2*d[5]*F(0,0)*F(0,2)*F(1,0)*F(2,0) + 2*d[4]*F(0,1)*F(0,2)*F(1,0)*F(2,0) + d[2]*pow(F(0,2),2)*F(1,0)*F(2,0) + d[18]*pow(F(0,0),2)*F(1,1)*F(2,0) + 2*d[21]*F(0,0)*F(0,1)*F(1,1)*F(2,0) + d[19]*pow(F(0,1),2)*F(1,1)*F(2,0) + 2*d[23]*F(0,0)*F(0,2)*F(1,1)*F(2,0) + 2*d[22]*F(0,1)*F(0,2)*F(1,1)*F(2,0) + d[20]*pow(F(0,2),2)*F(1,1)*F(2,0) + d[30]*pow(F(0,0),2)*F(1,2)*F(2,0) + 2*d[33]*F(0,0)*F(0,1)*F(1,2)*F(2,0) + d[31]*pow(F(0,1),2)*F(1,2)*F(2,0) + 2*d[35]*F(0,0)*F(0,2)*F(1,2)*F(2,0) + 2*d[34]*F(0,1)*F(0,2)*F(1,2)*F(2,0) + d[32]*pow(F(0,2),2)*F(1,2)*F(2,0) + d[18]*pow(F(0,0),2)*F(1,0)*F(2,1) + 2*d[21]*F(0,0)*F(0,1)*F(1,0)*F(2,1) + d[19]*pow(F(0,1),2)*F(1,0)*F(2,1) + 2*d[23]*F(0,0)*F(0,2)*F(1,0)*F(2,1) + 2*d[22]*F(0,1)*F(0,2)*F(1,0)*F(2,1) + d[20]*pow(F(0,2),2)*F(1,0)*F(2,1) + d[6]*pow(F(0,0),2)*F(1,1)*F(2,1) + 2*d[9]*F(0,0)*F(0,1)*F(1,1)*F(2,1) + d[7]*pow(F(0,1),2)*F(1,1)*F(2,1) + 2*d[11]*F(0,0)*F(0,2)*F(1,1)*F(2,1) + 2*d[10]*F(0,1)*F(0,2)*F(1,1)*F(2,1) + d[8]*pow(F(0,2),2)*F(1,1)*F(2,1) + d[24]*pow(F(0,0),2)*F(1,2)*F(2,1) + 2*d[27]*F(0,0)*F(0,1)*F(1,2)*F(2,1) + d[25]*pow(F(0,1),2)*F(1,2)*F(2,1) + 2*d[29]*F(0,0)*F(0,2)*F(1,2)*F(2,1) + 2*d[28]*F(0,1)*F(0,2)*F(1,2)*F(2,1) + d[26]*pow(F(0,2),2)*F(1,2)*F(2,1) + d[30]*pow(F(0,0),2)*F(1,0)*F(2,2) + 2*d[33]*F(0,0)*F(0,1)*F(1,0)*F(2,2) + d[31]*pow(F(0,1),2)*F(1,0)*F(2,2) + 2*d[35]*F(0,0)*F(0,2)*F(1,0)*F(2,2) + 2*d[34]*F(0,1)*F(0,2)*F(1,0)*F(2,2) + d[32]*pow(F(0,2),2)*F(1,0)*F(2,2) + d[24]*pow(F(0,0),2)*F(1,1)*F(2,2) + 2*d[27]*F(0,0)*F(0,1)*F(1,1)*F(2,2) + d[25]*pow(F(0,1),2)*F(1,1)*F(2,2) + 2*d[29]*F(0,0)*F(0,2)*F(1,1)*F(2,2) + 2*d[28]*F(0,1)*F(0,2)*F(1,1)*F(2,2) + d[26]*pow(F(0,2),2)*F(1,1)*F(2,2) + d[12]*pow(F(0,0),2)*F(1,2)*F(2,2) + 2*d[15]*F(0,0)*F(0,1)*F(1,2)*F(2,2) + d[13]*pow(F(0,1),2)*F(1,2)*F(2,2) + 2*d[17]*F(0,0)*F(0,2)*F(1,2)*F(2,2) + 2*d[16]*F(0,1)*F(0,2)*F(1,2)*F(2,2) + d[14]*pow(F(0,2),2)*F(1,2)*F(2,2);

    c.d[25] = d[0]*pow(F(1,0),3)*F(2,0) + 2*d[3]*pow(F(1,0),2)*F(1,1)*F(2,0) + d[18]*pow(F(1,0),2)*F(1,1)*F(2,0) + d[1]*F(1,0)*pow(F(1,1),2)*F(2,0) + 2*d[21]*F(1,0)*pow(F(1,1),2)*F(2,0) + d[19]*pow(F(1,1),3)*F(2,0) + 2*d[5]*pow(F(1,0),2)*F(1,2)*F(2,0) + d[30]*pow(F(1,0),2)*F(1,2)*F(2,0) + 2*d[4]*F(1,0)*F(1,1)*F(1,2)*F(2,0) + 2*d[23]*F(1,0)*F(1,1)*F(1,2)*F(2,0) + 2*d[33]*F(1,0)*F(1,1)*F(1,2)*F(2,0) + 2*d[22]*pow(F(1,1),2)*F(1,2)*F(2,0) + d[31]*pow(F(1,1),2)*F(1,2)*F(2,0) + d[2]*F(1,0)*pow(F(1,2),2)*F(2,0) + 2*d[35]*F(1,0)*pow(F(1,2),2)*F(2,0) + d[20]*F(1,1)*pow(F(1,2),2)*F(2,0) + 2*d[34]*F(1,1)*pow(F(1,2),2)*F(2,0) + d[32]*pow(F(1,2),3)*F(2,0) + d[18]*pow(F(1,0),3)*F(2,1) + d[6]*pow(F(1,0),2)*F(1,1)*F(2,1) + 2*d[21]*pow(F(1,0),2)*F(1,1)*F(2,1) + 2*d[9]*F(1,0)*pow(F(1,1),2)*F(2,1) + d[19]*F(1,0)*pow(F(1,1),2)*F(2,1) + d[7]*pow(F(1,1),3)*F(2,1) + 2*d[23]*pow(F(1,0),2)*F(1,2)*F(2,1) + d[24]*pow(F(1,0),2)*F(1,2)*F(2,1) + 2*d[11]*F(1,0)*F(1,1)*F(1,2)*F(2,1) + 2*d[22]*F(1,0)*F(1,1)*F(1,2)*F(2,1) + 2*d[27]*F(1,0)*F(1,1)*F(1,2)*F(2,1) + 2*d[10]*pow(F(1,1),2)*F(1,2)*F(2,1) + d[25]*pow(F(1,1),2)*F(1,2)*F(2,1) + d[20]*F(1,0)*pow(F(1,2),2)*F(2,1) + 2*d[29]*F(1,0)*pow(F(1,2),2)*F(2,1) + d[8]*F(1,1)*pow(F(1,2),2)*F(2,1) + 2*d[28]*F(1,1)*pow(F(1,2),2)*F(2,1) + d[26]*pow(F(1,2),3)*F(2,1) + d[30]*pow(F(1,0),3)*F(2,2) + d[24]*pow(F(1,0),2)*F(1,1)*F(2,2) + 2*d[33]*pow(F(1,0),2)*F(1,1)*F(2,2) + 2*d[27]*F(1,0)*pow(F(1,1),2)*F(2,2) + d[31]*F(1,0)*pow(F(1,1),2)*F(2,2) + d[25]*pow(F(1,1),3)*F(2,2) + d[12]*pow(F(1,0),2)*F(1,2)*F(2,2) + 2*d[35]*pow(F(1,0),2)*F(1,2)*F(2,2) + 2*d[15]*F(1,0)*F(1,1)*F(1,2)*F(2,2) + 2*d[29]*F(1,0)*F(1,1)*F(1,2)*F(2,2) + 2*d[34]*F(1,0)*F(1,1)*F(1,2)*F(2,2) + d[13]*pow(F(1,1),2)*F(1,2)*F(2,2) + 2*d[28]*pow(F(1,1),2)*F(1,2)*F(2,2) + 2*d[17]*F(1,0)*pow(F(1,2),2)*F(2,2) + d[32]*F(1,0)*pow(F(1,2),2)*F(2,2) + 2*d[16]*F(1,1)*pow(F(1,2),2)*F(2,2) + d[26]*F(1,1)*pow(F(1,2),2)*F(2,2) + d[14]*pow(F(1,2),3)*F(2,2);

    c.d[26] = d[0]*F(1,0)*pow(F(2,0),3) + d[30]*F(1,2)*pow(F(2,0),3) + 2*d[3]*F(1,0)*pow(F(2,0),2)*F(2,1) + d[6]*F(1,1)*pow(F(2,0),2)*F(2,1) + 2*d[21]*F(1,1)*pow(F(2,0),2)*F(2,1) + d[24]*F(1,2)*pow(F(2,0),2)*F(2,1) + 2*d[33]*F(1,2)*pow(F(2,0),2)*F(2,1) + d[1]*F(1,0)*F(2,0)*pow(F(2,1),2) + 2*d[21]*F(1,0)*F(2,0)*pow(F(2,1),2) + 2*d[9]*F(1,1)*F(2,0)*pow(F(2,1),2) + d[19]*F(1,1)*F(2,0)*pow(F(2,1),2) + 2*d[27]*F(1,2)*F(2,0)*pow(F(2,1),2) + d[31]*F(1,2)*F(2,0)*pow(F(2,1),2) + d[19]*F(1,0)*pow(F(2,1),3) + d[7]*F(1,1)*pow(F(2,1),3) + d[25]*F(1,2)*pow(F(2,1),3) + d[18]*pow(F(2,0),2)*(F(1,1)*F(2,0) + F(1,0)*F(2,1)) + 2*d[5]*F(1,0)*pow(F(2,0),2)*F(2,2) + d[30]*F(1,0)*pow(F(2,0),2)*F(2,2) + 2*d[23]*F(1,1)*pow(F(2,0),2)*F(2,2) + d[24]*F(1,1)*pow(F(2,0),2)*F(2,2) + d[12]*F(1,2)*pow(F(2,0),2)*F(2,2) + 2*d[35]*F(1,2)*pow(F(2,0),2)*F(2,2) + 2*d[4]*F(1,0)*F(2,0)*F(2,1)*F(2,2) + 2*d[23]*F(1,0)*F(2,0)*F(2,1)*F(2,2) + 2*d[33]*F(1,0)*F(2,0)*F(2,1)*F(2,2) + 2*d[11]*F(1,1)*F(2,0)*F(2,1)*F(2,2) + 2*d[22]*F(1,1)*F(2,0)*F(2,1)*F(2,2) + 2*d[27]*F(1,1)*F(2,0)*F(2,1)*F(2,2) + 2*d[15]*F(1,2)*F(2,0)*F(2,1)*F(2,2) + 2*d[29]*F(1,2)*F(2,0)*F(2,1)*F(2,2) + 2*d[34]*F(1,2)*F(2,0)*F(2,1)*F(2,2) + 2*d[22]*F(1,0)*pow(F(2,1),2)*F(2,2) + d[31]*F(1,0)*pow(F(2,1),2)*F(2,2) + 2*d[10]*F(1,1)*pow(F(2,1),2)*F(2,2) + d[25]*F(1,1)*pow(F(2,1),2)*F(2,2) + d[13]*F(1,2)*pow(F(2,1),2)*F(2,2) + 2*d[28]*F(1,2)*pow(F(2,1),2)*F(2,2) + d[2]*F(1,0)*F(2,0)*pow(F(2,2),2) + 2*d[35]*F(1,0)*F(2,0)*pow(F(2,2),2) + d[20]*F(1,1)*F(2,0)*pow(F(2,2),2) + 2*d[29]*F(1,1)*F(2,0)*pow(F(2,2),2) + 2*d[17]*F(1,2)*F(2,0)*pow(F(2,2),2) + d[32]*F(1,2)*F(2,0)*pow(F(2,2),2) + d[20]*F(1,0)*F(2,1)*pow(F(2,2),2) + 2*d[34]*F(1,0)*F(2,1)*pow(F(2,2),2) + d[8]*F(1,1)*F(2,1)*pow(F(2,2),2) + 2*d[28]*F(1,1)*F(2,1)*pow(F(2,2),2) + 2*d[16]*F(1,2)*F(2,1)*pow(F(2,2),2) + d[26]*F(1,2)*F(2,1)*pow(F(2,2),2) + d[32]*F(1,0)*pow(F(2,2),3) + d[26]*F(1,1)*pow(F(2,2),3) + d[14]*F(1,2)*pow(F(2,2),3);

    c.d[27] = d[0]*F(0,0)*pow(F(1,0),2)*F(2,0) + d[5]*F(0,2)*pow(F(1,0),2)*F(2,0) + d[18]*F(0,0)*F(1,0)*F(1,1)*F(2,0) + d[1]*F(0,1)*F(1,0)*F(1,1)*F(2,0) + d[21]*F(0,1)*F(1,0)*F(1,1)*F(2,0) + d[4]*F(0,2)*F(1,0)*F(1,1)*F(2,0) + d[23]*F(0,2)*F(1,0)*F(1,1)*F(2,0) + d[21]*F(0,0)*pow(F(1,1),2)*F(2,0) + d[19]*F(0,1)*pow(F(1,1),2)*F(2,0) + d[22]*F(0,2)*pow(F(1,1),2)*F(2,0) + d[3]*F(1,0)*(F(0,1)*F(1,0) + F(0,0)*F(1,1))*F(2,0) + d[5]*F(0,0)*F(1,0)*F(1,2)*F(2,0) + d[30]*F(0,0)*F(1,0)*F(1,2)*F(2,0) + d[4]*F(0,1)*F(1,0)*F(1,2)*F(2,0) + d[33]*F(0,1)*F(1,0)*F(1,2)*F(2,0) + d[2]*F(0,2)*F(1,0)*F(1,2)*F(2,0) + d[35]*F(0,2)*F(1,0)*F(1,2)*F(2,0) + d[23]*F(0,0)*F(1,1)*F(1,2)*F(2,0) + d[33]*F(0,0)*F(1,1)*F(1,2)*F(2,0) + d[22]*F(0,1)*F(1,1)*F(1,2)*F(2,0) + d[31]*F(0,1)*F(1,1)*F(1,2)*F(2,0) + d[20]*F(0,2)*F(1,1)*F(1,2)*F(2,0) + d[34]*F(0,2)*F(1,1)*F(1,2)*F(2,0) + d[35]*F(0,0)*pow(F(1,2),2)*F(2,0) + d[34]*F(0,1)*pow(F(1,2),2)*F(2,0) + d[32]*F(0,2)*pow(F(1,2),2)*F(2,0) + d[18]*F(0,0)*pow(F(1,0),2)*F(2,1) + d[21]*F(0,1)*pow(F(1,0),2)*F(2,1) + d[23]*F(0,2)*pow(F(1,0),2)*F(2,1) + d[6]*F(0,0)*F(1,0)*F(1,1)*F(2,1) + d[21]*F(0,0)*F(1,0)*F(1,1)*F(2,1) + d[9]*F(0,1)*F(1,0)*F(1,1)*F(2,1) + d[19]*F(0,1)*F(1,0)*F(1,1)*F(2,1) + d[11]*F(0,2)*F(1,0)*F(1,1)*F(2,1) + d[22]*F(0,2)*F(1,0)*F(1,1)*F(2,1) + d[9]*F(0,0)*pow(F(1,1),2)*F(2,1) + d[7]*F(0,1)*pow(F(1,1),2)*F(2,1) + d[10]*F(0,2)*pow(F(1,1),2)*F(2,1) + d[23]*F(0,0)*F(1,0)*F(1,2)*F(2,1) + d[24]*F(0,0)*F(1,0)*F(1,2)*F(2,1) + d[22]*F(0,1)*F(1,0)*F(1,2)*F(2,1) + d[27]*F(0,1)*F(1,0)*F(1,2)*F(2,1) + d[20]*F(0,2)*F(1,0)*F(1,2)*F(2,1) + d[29]*F(0,2)*F(1,0)*F(1,2)*F(2,1) + d[11]*F(0,0)*F(1,1)*F(1,2)*F(2,1) + d[27]*F(0,0)*F(1,1)*F(1,2)*F(2,1) + d[10]*F(0,1)*F(1,1)*F(1,2)*F(2,1) + d[25]*F(0,1)*F(1,1)*F(1,2)*F(2,1) + d[8]*F(0,2)*F(1,1)*F(1,2)*F(2,1) + d[28]*F(0,2)*F(1,1)*F(1,2)*F(2,1) + d[29]*F(0,0)*pow(F(1,2),2)*F(2,1) + d[28]*F(0,1)*pow(F(1,2),2)*F(2,1) + d[26]*F(0,2)*pow(F(1,2),2)*F(2,1) + d[30]*F(0,0)*pow(F(1,0),2)*F(2,2) + d[33]*F(0,1)*pow(F(1,0),2)*F(2,2) + d[35]*F(0,2)*pow(F(1,0),2)*F(2,2) + d[24]*F(0,0)*F(1,0)*F(1,1)*F(2,2) + d[33]*F(0,0)*F(1,0)*F(1,1)*F(2,2) + d[27]*F(0,1)*F(1,0)*F(1,1)*F(2,2) + d[31]*F(0,1)*F(1,0)*F(1,1)*F(2,2) + d[29]*F(0,2)*F(1,0)*F(1,1)*F(2,2) + d[34]*F(0,2)*F(1,0)*F(1,1)*F(2,2) + d[27]*F(0,0)*pow(F(1,1),2)*F(2,2) + d[25]*F(0,1)*pow(F(1,1),2)*F(2,2) + d[28]*F(0,2)*pow(F(1,1),2)*F(2,2) + d[12]*F(0,0)*F(1,0)*F(1,2)*F(2,2) + d[35]*F(0,0)*F(1,0)*F(1,2)*F(2,2) + d[15]*F(0,1)*F(1,0)*F(1,2)*F(2,2) + d[34]*F(0,1)*F(1,0)*F(1,2)*F(2,2) + d[17]*F(0,2)*F(1,0)*F(1,2)*F(2,2) + d[32]*F(0,2)*F(1,0)*F(1,2)*F(2,2) + d[15]*F(0,0)*F(1,1)*F(1,2)*F(2,2) + d[29]*F(0,0)*F(1,1)*F(1,2)*F(2,2) + d[13]*F(0,1)*F(1,1)*F(1,2)*F(2,2) + d[28]*F(0,1)*F(1,1)*F(1,2)*F(2,2) + d[16]*F(0,2)*F(1,1)*F(1,2)*F(2,2) + d[26]*F(0,2)*F(1,1)*F(1,2)*F(2,2) + d[17]*F(0,0)*pow(F(1,2),2)*F(2,2) + d[16]*F(0,1)*pow(F(1,2),2)*F(2,2) + d[14]*F(0,2)*pow(F(1,2),2)*F(2,2);

    c.d[28] = d[0]*pow(F(1,0),2)*pow(F(2,0),2) + d[18]*F(1,0)*F(1,1)*pow(F(2,0),2) + d[21]*pow(F(1,1),2)*pow(F(2,0),2) + d[5]*F(1,0)*F(1,2)*pow(F(2,0),2) + d[30]*F(1,0)*F(1,2)*pow(F(2,0),2) + d[23]*F(1,1)*F(1,2)*pow(F(2,0),2) + d[33]*F(1,1)*F(1,2)*pow(F(2,0),2) + d[35]*pow(F(1,2),2)*pow(F(2,0),2) + d[18]*pow(F(1,0),2)*F(2,0)*F(2,1) + d[1]*F(1,0)*F(1,1)*F(2,0)*F(2,1) + d[6]*F(1,0)*F(1,1)*F(2,0)*F(2,1) + 2*d[21]*F(1,0)*F(1,1)*F(2,0)*F(2,1) + d[9]*pow(F(1,1),2)*F(2,0)*F(2,1) + d[19]*pow(F(1,1),2)*F(2,0)*F(2,1) + d[4]*F(1,0)*F(1,2)*F(2,0)*F(2,1) + d[23]*F(1,0)*F(1,2)*F(2,0)*F(2,1) + d[24]*F(1,0)*F(1,2)*F(2,0)*F(2,1) + d[33]*F(1,0)*F(1,2)*F(2,0)*F(2,1) + d[11]*F(1,1)*F(1,2)*F(2,0)*F(2,1) + d[22]*F(1,1)*F(1,2)*F(2,0)*F(2,1) + d[27]*F(1,1)*F(1,2)*F(2,0)*F(2,1) + d[31]*F(1,1)*F(1,2)*F(2,0)*F(2,1) + d[29]*pow(F(1,2),2)*F(2,0)*F(2,1) + d[34]*pow(F(1,2),2)*F(2,0)*F(2,1) + d[21]*pow(F(1,0),2)*pow(F(2,1),2) + d[9]*F(1,0)*F(1,1)*pow(F(2,1),2) + d[19]*F(1,0)*F(1,1)*pow(F(2,1),2) + d[7]*pow(F(1,1),2)*pow(F(2,1),2) + d[22]*F(1,0)*F(1,2)*pow(F(2,1),2) + d[27]*F(1,0)*F(1,2)*pow(F(2,1),2) + d[10]*F(1,1)*F(1,2)*pow(F(2,1),2) + d[25]*F(1,1)*F(1,2)*pow(F(2,1),2) + d[28]*pow(F(1,2),2)*pow(F(2,1),2) + d[3]*F(1,0)*F(2,0)*(F(1,1)*F(2,0) + F(1,0)*F(2,1)) + d[5]*pow(F(1,0),2)*F(2,0)*F(2,2) + d[30]*pow(F(1,0),2)*F(2,0)*F(2,2) + d[4]*F(1,0)*F(1,1)*F(2,0)*F(2,2) + d[23]*F(1,0)*F(1,1)*F(2,0)*F(2,2) + d[24]*F(1,0)*F(1,1)*F(2,0)*F(2,2) + d[33]*F(1,0)*F(1,1)*F(2,0)*F(2,2) + d[22]*pow(F(1,1),2)*F(2,0)*F(2,2) + d[27]*pow(F(1,1),2)*F(2,0)*F(2,2) + d[2]*F(1,0)*F(1,2)*F(2,0)*F(2,2) + d[12]*F(1,0)*F(1,2)*F(2,0)*F(2,2) + 2*d[35]*F(1,0)*F(1,2)*F(2,0)*F(2,2) + d[15]*F(1,1)*F(1,2)*F(2,0)*F(2,2) + d[20]*F(1,1)*F(1,2)*F(2,0)*F(2,2) + d[29]*F(1,1)*F(1,2)*F(2,0)*F(2,2) + d[34]*F(1,1)*F(1,2)*F(2,0)*F(2,2) + d[17]*pow(F(1,2),2)*F(2,0)*F(2,2) + d[32]*pow(F(1,2),2)*F(2,0)*F(2,2) + d[23]*pow(F(1,0),2)*F(2,1)*F(2,2) + d[33]*pow(F(1,0),2)*F(2,1)*F(2,2) + d[11]*F(1,0)*F(1,1)*F(2,1)*F(2,2) + d[22]*F(1,0)*F(1,1)*F(2,1)*F(2,2) + d[27]*F(1,0)*F(1,1)*F(2,1)*F(2,2) + d[31]*F(1,0)*F(1,1)*F(2,1)*F(2,2) + d[10]*pow(F(1,1),2)*F(2,1)*F(2,2) + d[25]*pow(F(1,1),2)*F(2,1)*F(2,2) + d[15]*F(1,0)*F(1,2)*F(2,1)*F(2,2) + d[20]*F(1,0)*F(1,2)*F(2,1)*F(2,2) + d[29]*F(1,0)*F(1,2)*F(2,1)*F(2,2) + d[34]*F(1,0)*F(1,2)*F(2,1)*F(2,2) + d[8]*F(1,1)*F(1,2)*F(2,1)*F(2,2) + d[13]*F(1,1)*F(1,2)*F(2,1)*F(2,2) + 2*d[28]*F(1,1)*F(1,2)*F(2,1)*F(2,2) + d[16]*pow(F(1,2),2)*F(2,1)*F(2,2) + d[26]*pow(F(1,2),2)*F(2,1)*F(2,2) + d[35]*pow(F(1,0),2)*pow(F(2,2),2) + d[29]*F(1,0)*F(1,1)*pow(F(2,2),2) + d[34]*F(1,0)*F(1,1)*pow(F(2,2),2) + d[28]*pow(F(1,1),2)*pow(F(2,2),2) + d[17]*F(1,0)*F(1,2)*pow(F(2,2),2) + d[32]*F(1,0)*F(1,2)*pow(F(2,2),2) + d[16]*F(1,1)*F(1,2)*pow(F(2,2),2) + d[26]*F(1,1)*F(1,2)*pow(F(2,2),2) + d[14]*pow(F(1,2),2)*pow(F(2,2),2);

    c.d[29] = d[0]*F(0,0)*F(1,0)*pow(F(2,0),2) + d[5]*F(0,2)*F(1,0)*pow(F(2,0),2) + d[18]*F(0,0)*F(1,1)*pow(F(2,0),2) + d[21]*F(0,1)*F(1,1)*pow(F(2,0),2) + d[23]*F(0,2)*F(1,1)*pow(F(2,0),2) + d[30]*F(0,0)*F(1,2)*pow(F(2,0),2) + d[33]*F(0,1)*F(1,2)*pow(F(2,0),2) + d[35]*F(0,2)*F(1,2)*pow(F(2,0),2) + d[18]*F(0,0)*F(1,0)*F(2,0)*F(2,1) + d[1]*F(0,1)*F(1,0)*F(2,0)*F(2,1) + d[21]*F(0,1)*F(1,0)*F(2,0)*F(2,1) + d[4]*F(0,2)*F(1,0)*F(2,0)*F(2,1) + d[23]*F(0,2)*F(1,0)*F(2,0)*F(2,1) + d[6]*F(0,0)*F(1,1)*F(2,0)*F(2,1) + d[21]*F(0,0)*F(1,1)*F(2,0)*F(2,1) + d[9]*F(0,1)*F(1,1)*F(2,0)*F(2,1) + d[19]*F(0,1)*F(1,1)*F(2,0)*F(2,1) + d[11]*F(0,2)*F(1,1)*F(2,0)*F(2,1) + d[22]*F(0,2)*F(1,1)*F(2,0)*F(2,1) + d[24]*F(0,0)*F(1,2)*F(2,0)*F(2,1) + d[33]*F(0,0)*F(1,2)*F(2,0)*F(2,1) + d[27]*F(0,1)*F(1,2)*F(2,0)*F(2,1) + d[31]*F(0,1)*F(1,2)*F(2,0)*F(2,1) + d[29]*F(0,2)*F(1,2)*F(2,0)*F(2,1) + d[34]*F(0,2)*F(1,2)*F(2,0)*F(2,1) + d[21]*F(0,0)*F(1,0)*pow(F(2,1),2) + d[19]*F(0,1)*F(1,0)*pow(F(2,1),2) + d[22]*F(0,2)*F(1,0)*pow(F(2,1),2) + d[9]*F(0,0)*F(1,1)*pow(F(2,1),2) + d[7]*F(0,1)*F(1,1)*pow(F(2,1),2) + d[10]*F(0,2)*F(1,1)*pow(F(2,1),2) + d[27]*F(0,0)*F(1,2)*pow(F(2,1),2) + d[25]*F(0,1)*F(1,2)*pow(F(2,1),2) + d[28]*F(0,2)*F(1,2)*pow(F(2,1),2) + d[3]*F(1,0)*F(2,0)*(F(0,1)*F(2,0) + F(0,0)*F(2,1)) + d[5]*F(0,0)*F(1,0)*F(2,0)*F(2,2) + d[30]*F(0,0)*F(1,0)*F(2,0)*F(2,2) + d[4]*F(0,1)*F(1,0)*F(2,0)*F(2,2) + d[33]*F(0,1)*F(1,0)*F(2,0)*F(2,2) + d[2]*F(0,2)*F(1,0)*F(2,0)*F(2,2) + d[35]*F(0,2)*F(1,0)*F(2,0)*F(2,2) + d[23]*F(0,0)*F(1,1)*F(2,0)*F(2,2) + d[24]*F(0,0)*F(1,1)*F(2,0)*F(2,2) + d[22]*F(0,1)*F(1,1)*F(2,0)*F(2,2) + d[27]*F(0,1)*F(1,1)*F(2,0)*F(2,2) + d[20]*F(0,2)*F(1,1)*F(2,0)*F(2,2) + d[29]*F(0,2)*F(1,1)*F(2,0)*F(2,2) + d[12]*F(0,0)*F(1,2)*F(2,0)*F(2,2) + d[35]*F(0,0)*F(1,2)*F(2,0)*F(2,2) + d[15]*F(0,1)*F(1,2)*F(2,0)*F(2,2) + d[34]*F(0,1)*F(1,2)*F(2,0)*F(2,2) + d[17]*F(0,2)*F(1,2)*F(2,0)*F(2,2) + d[32]*F(0,2)*F(1,2)*F(2,0)*F(2,2) + d[23]*F(0,0)*F(1,0)*F(2,1)*F(2,2) + d[33]*F(0,0)*F(1,0)*F(2,1)*F(2,2) + d[22]*F(0,1)*F(1,0)*F(2,1)*F(2,2) + d[31]*F(0,1)*F(1,0)*F(2,1)*F(2,2) + d[20]*F(0,2)*F(1,0)*F(2,1)*F(2,2) + d[34]*F(0,2)*F(1,0)*F(2,1)*F(2,2) + d[11]*F(0,0)*F(1,1)*F(2,1)*F(2,2) + d[27]*F(0,0)*F(1,1)*F(2,1)*F(2,2) + d[10]*F(0,1)*F(1,1)*F(2,1)*F(2,2) + d[25]*F(0,1)*F(1,1)*F(2,1)*F(2,2) + d[8]*F(0,2)*F(1,1)*F(2,1)*F(2,2) + d[28]*F(0,2)*F(1,1)*F(2,1)*F(2,2) + d[15]*F(0,0)*F(1,2)*F(2,1)*F(2,2) + d[29]*F(0,0)*F(1,2)*F(2,1)*F(2,2) + d[13]*F(0,1)*F(1,2)*F(2,1)*F(2,2) + d[28]*F(0,1)*F(1,2)*F(2,1)*F(2,2) + d[16]*F(0,2)*F(1,2)*F(2,1)*F(2,2) + d[26]*F(0,2)*F(1,2)*F(2,1)*F(2,2) + d[35]*F(0,0)*F(1,0)*pow(F(2,2),2) + d[34]*F(0,1)*F(1,0)*pow(F(2,2),2) + d[32]*F(0,2)*F(1,0)*pow(F(2,2),2) + d[29]*F(0,0)*F(1,1)*pow(F(2,2),2) + d[28]*F(0,1)*F(1,1)*pow(F(2,2),2) + d[26]*F(0,2)*F(1,1)*pow(F(2,2),2) + d[17]*F(0,0)*F(1,2)*pow(F(2,2),2) + d[16]*F(0,1)*F(1,2)*pow(F(2,2),2) + d[14]*F(0,2)*F(1,2)*pow(F(2,2),2);

    c.d[30] = d[0]*pow(F(0,0),3)*F(2,0) + 2*d[3]*pow(F(0,0),2)*F(0,1)*F(2,0) + d[18]*pow(F(0,0),2)*F(0,1)*F(2,0) + d[1]*F(0,0)*pow(F(0,1),2)*F(2,0) + 2*d[21]*F(0,0)*pow(F(0,1),2)*F(2,0) + d[19]*pow(F(0,1),3)*F(2,0) + 2*d[5]*pow(F(0,0),2)*F(0,2)*F(2,0) + d[30]*pow(F(0,0),2)*F(0,2)*F(2,0) + 2*d[4]*F(0,0)*F(0,1)*F(0,2)*F(2,0) + 2*d[23]*F(0,0)*F(0,1)*F(0,2)*F(2,0) + 2*d[33]*F(0,0)*F(0,1)*F(0,2)*F(2,0) + 2*d[22]*pow(F(0,1),2)*F(0,2)*F(2,0) + d[31]*pow(F(0,1),2)*F(0,2)*F(2,0) + d[2]*F(0,0)*pow(F(0,2),2)*F(2,0) + 2*d[35]*F(0,0)*pow(F(0,2),2)*F(2,0) + d[20]*F(0,1)*pow(F(0,2),2)*F(2,0) + 2*d[34]*F(0,1)*pow(F(0,2),2)*F(2,0) + d[32]*pow(F(0,2),3)*F(2,0) + d[18]*pow(F(0,0),3)*F(2,1) + d[6]*pow(F(0,0),2)*F(0,1)*F(2,1) + 2*d[21]*pow(F(0,0),2)*F(0,1)*F(2,1) + 2*d[9]*F(0,0)*pow(F(0,1),2)*F(2,1) + d[19]*F(0,0)*pow(F(0,1),2)*F(2,1) + d[7]*pow(F(0,1),3)*F(2,1) + 2*d[23]*pow(F(0,0),2)*F(0,2)*F(2,1) + d[24]*pow(F(0,0),2)*F(0,2)*F(2,1) + 2*d[11]*F(0,0)*F(0,1)*F(0,2)*F(2,1) + 2*d[22]*F(0,0)*F(0,1)*F(0,2)*F(2,1) + 2*d[27]*F(0,0)*F(0,1)*F(0,2)*F(2,1) + 2*d[10]*pow(F(0,1),2)*F(0,2)*F(2,1) + d[25]*pow(F(0,1),2)*F(0,2)*F(2,1) + d[20]*F(0,0)*pow(F(0,2),2)*F(2,1) + 2*d[29]*F(0,0)*pow(F(0,2),2)*F(2,1) + d[8]*F(0,1)*pow(F(0,2),2)*F(2,1) + 2*d[28]*F(0,1)*pow(F(0,2),2)*F(2,1) + d[26]*pow(F(0,2),3)*F(2,1) + d[30]*pow(F(0,0),3)*F(2,2) + d[24]*pow(F(0,0),2)*F(0,1)*F(2,2) + 2*d[33]*pow(F(0,0),2)*F(0,1)*F(2,2) + 2*d[27]*F(0,0)*pow(F(0,1),2)*F(2,2) + d[31]*F(0,0)*pow(F(0,1),2)*F(2,2) + d[25]*pow(F(0,1),3)*F(2,2) + d[12]*pow(F(0,0),2)*F(0,2)*F(2,2) + 2*d[35]*pow(F(0,0),2)*F(0,2)*F(2,2) + 2*d[15]*F(0,0)*F(0,1)*F(0,2)*F(2,2) + 2*d[29]*F(0,0)*F(0,1)*F(0,2)*F(2,2) + 2*d[34]*F(0,0)*F(0,1)*F(0,2)*F(2,2) + d[13]*pow(F(0,1),2)*F(0,2)*F(2,2) + 2*d[28]*pow(F(0,1),2)*F(0,2)*F(2,2) + 2*d[17]*F(0,0)*pow(F(0,2),2)*F(2,2) + d[32]*F(0,0)*pow(F(0,2),2)*F(2,2) + 2*d[16]*F(0,1)*pow(F(0,2),2)*F(2,2) + d[26]*F(0,1)*pow(F(0,2),2)*F(2,2) + d[14]*pow(F(0,2),3)*F(2,2);

    c.d[31] = d[0]*F(0,0)*pow(F(1,0),2)*F(2,0) + d[30]*F(0,2)*pow(F(1,0),2)*F(2,0) + 2*d[3]*F(0,0)*F(1,0)*F(1,1)*F(2,0) + 2*d[21]*F(0,1)*F(1,0)*F(1,1)*F(2,0) + 2*d[33]*F(0,2)*F(1,0)*F(1,1)*F(2,0) + d[1]*F(0,0)*pow(F(1,1),2)*F(2,0) + d[19]*F(0,1)*pow(F(1,1),2)*F(2,0) + d[31]*F(0,2)*pow(F(1,1),2)*F(2,0) + 2*d[5]*F(0,0)*F(1,0)*F(1,2)*F(2,0) + 2*d[23]*F(0,1)*F(1,0)*F(1,2)*F(2,0) + 2*d[35]*F(0,2)*F(1,0)*F(1,2)*F(2,0) + 2*d[4]*F(0,0)*F(1,1)*F(1,2)*F(2,0) + 2*d[22]*F(0,1)*F(1,1)*F(1,2)*F(2,0) + 2*d[34]*F(0,2)*F(1,1)*F(1,2)*F(2,0) + d[2]*F(0,0)*pow(F(1,2),2)*F(2,0) + d[20]*F(0,1)*pow(F(1,2),2)*F(2,0) + d[32]*F(0,2)*pow(F(1,2),2)*F(2,0) + d[6]*F(0,1)*pow(F(1,0),2)*F(2,1) + d[24]*F(0,2)*pow(F(1,0),2)*F(2,1) + 2*d[21]*F(0,0)*F(1,0)*F(1,1)*F(2,1) + 2*d[9]*F(0,1)*F(1,0)*F(1,1)*F(2,1) + 2*d[27]*F(0,2)*F(1,0)*F(1,1)*F(2,1) + d[19]*F(0,0)*pow(F(1,1),2)*F(2,1) + d[7]*F(0,1)*pow(F(1,1),2)*F(2,1) + d[25]*F(0,2)*pow(F(1,1),2)*F(2,1) + 2*d[23]*F(0,0)*F(1,0)*F(1,2)*F(2,1) + 2*d[11]*F(0,1)*F(1,0)*F(1,2)*F(2,1) + 2*d[29]*F(0,2)*F(1,0)*F(1,2)*F(2,1) + 2*d[22]*F(0,0)*F(1,1)*F(1,2)*F(2,1) + 2*d[10]*F(0,1)*F(1,1)*F(1,2)*F(2,1) + 2*d[28]*F(0,2)*F(1,1)*F(1,2)*F(2,1) + d[20]*F(0,0)*pow(F(1,2),2)*F(2,1) + d[8]*F(0,1)*pow(F(1,2),2)*F(2,1) + d[26]*F(0,2)*pow(F(1,2),2)*F(2,1) + d[18]*pow(F(1,0),2)*(F(0,1)*F(2,0) + F(0,0)*F(2,1)) + d[30]*F(0,0)*pow(F(1,0),2)*F(2,2) + d[24]*F(0,1)*pow(F(1,0),2)*F(2,2) + d[12]*F(0,2)*pow(F(1,0),2)*F(2,2) + 2*d[33]*F(0,0)*F(1,0)*F(1,1)*F(2,2) + 2*d[27]*F(0,1)*F(1,0)*F(1,1)*F(2,2) + 2*d[15]*F(0,2)*F(1,0)*F(1,1)*F(2,2) + d[31]*F(0,0)*pow(F(1,1),2)*F(2,2) + d[25]*F(0,1)*pow(F(1,1),2)*F(2,2) + d[13]*F(0,2)*pow(F(1,1),2)*F(2,2) + 2*d[35]*F(0,0)*F(1,0)*F(1,2)*F(2,2) + 2*d[29]*F(0,1)*F(1,0)*F(1,2)*F(2,2) + 2*d[17]*F(0,2)*F(1,0)*F(1,2)*F(2,2) + 2*d[34]*F(0,0)*F(1,1)*F(1,2)*F(2,2) + 2*d[28]*F(0,1)*F(1,1)*F(1,2)*F(2,2) + 2*d[16]*F(0,2)*F(1,1)*F(1,2)*F(2,2) + d[32]*F(0,0)*pow(F(1,2),2)*F(2,2) + d[26]*F(0,1)*pow(F(1,2),2)*F(2,2) + d[14]*F(0,2)*pow(F(1,2),2)*F(2,2);

    c.d[32] = d[0]*F(0,0)*pow(F(2,0),3) + d[30]*F(0,2)*pow(F(2,0),3) + 2*d[3]*F(0,0)*pow(F(2,0),2)*F(2,1) + d[6]*F(0,1)*pow(F(2,0),2)*F(2,1) + 2*d[21]*F(0,1)*pow(F(2,0),2)*F(2,1) + d[24]*F(0,2)*pow(F(2,0),2)*F(2,1) + 2*d[33]*F(0,2)*pow(F(2,0),2)*F(2,1) + d[1]*F(0,0)*F(2,0)*pow(F(2,1),2) + 2*d[21]*F(0,0)*F(2,0)*pow(F(2,1),2) + 2*d[9]*F(0,1)*F(2,0)*pow(F(2,1),2) + d[19]*F(0,1)*F(2,0)*pow(F(2,1),2) + 2*d[27]*F(0,2)*F(2,0)*pow(F(2,1),2) + d[31]*F(0,2)*F(2,0)*pow(F(2,1),2) + d[19]*F(0,0)*pow(F(2,1),3) + d[7]*F(0,1)*pow(F(2,1),3) + d[25]*F(0,2)*pow(F(2,1),3) + d[18]*pow(F(2,0),2)*(F(0,1)*F(2,0) + F(0,0)*F(2,1)) + 2*d[5]*F(0,0)*pow(F(2,0),2)*F(2,2) + d[30]*F(0,0)*pow(F(2,0),2)*F(2,2) + 2*d[23]*F(0,1)*pow(F(2,0),2)*F(2,2) + d[24]*F(0,1)*pow(F(2,0),2)*F(2,2) + d[12]*F(0,2)*pow(F(2,0),2)*F(2,2) + 2*d[35]*F(0,2)*pow(F(2,0),2)*F(2,2) + 2*d[4]*F(0,0)*F(2,0)*F(2,1)*F(2,2) + 2*d[23]*F(0,0)*F(2,0)*F(2,1)*F(2,2) + 2*d[33]*F(0,0)*F(2,0)*F(2,1)*F(2,2) + 2*d[11]*F(0,1)*F(2,0)*F(2,1)*F(2,2) + 2*d[22]*F(0,1)*F(2,0)*F(2,1)*F(2,2) + 2*d[27]*F(0,1)*F(2,0)*F(2,1)*F(2,2) + 2*d[15]*F(0,2)*F(2,0)*F(2,1)*F(2,2) + 2*d[29]*F(0,2)*F(2,0)*F(2,1)*F(2,2) + 2*d[34]*F(0,2)*F(2,0)*F(2,1)*F(2,2) + 2*d[22]*F(0,0)*pow(F(2,1),2)*F(2,2) + d[31]*F(0,0)*pow(F(2,1),2)*F(2,2) + 2*d[10]*F(0,1)*pow(F(2,1),2)*F(2,2) + d[25]*F(0,1)*pow(F(2,1),2)*F(2,2) + d[13]*F(0,2)*pow(F(2,1),2)*F(2,2) + 2*d[28]*F(0,2)*pow(F(2,1),2)*F(2,2) + d[2]*F(0,0)*F(2,0)*pow(F(2,2),2) + 2*d[35]*F(0,0)*F(2,0)*pow(F(2,2),2) + d[20]*F(0,1)*F(2,0)*pow(F(2,2),2) + 2*d[29]*F(0,1)*F(2,0)*pow(F(2,2),2) + 2*d[17]*F(0,2)*F(2,0)*pow(F(2,2),2) + d[32]*F(0,2)*F(2,0)*pow(F(2,2),2) + d[20]*F(0,0)*F(2,1)*pow(F(2,2),2) + 2*d[34]*F(0,0)*F(2,1)*pow(F(2,2),2) + d[8]*F(0,1)*F(2,1)*pow(F(2,2),2) + 2*d[28]*F(0,1)*F(2,1)*pow(F(2,2),2) + 2*d[16]*F(0,2)*F(2,1)*pow(F(2,2),2) + d[26]*F(0,2)*F(2,1)*pow(F(2,2),2) + d[32]*F(0,0)*pow(F(2,2),3) + d[26]*F(0,1)*pow(F(2,2),3) + d[14]*F(0,2)*pow(F(2,2),3);

    c.d[33] = d[0]*pow(F(0,0),2)*F(1,0)*F(2,0) + d[18]*F(0,0)*F(0,1)*F(1,0)*F(2,0) + d[21]*pow(F(0,1),2)*F(1,0)*F(2,0) + d[5]*F(0,0)*F(0,2)*F(1,0)*F(2,0) + d[30]*F(0,0)*F(0,2)*F(1,0)*F(2,0) + d[23]*F(0,1)*F(0,2)*F(1,0)*F(2,0) + d[33]*F(0,1)*F(0,2)*F(1,0)*F(2,0) + d[35]*pow(F(0,2),2)*F(1,0)*F(2,0) + d[1]*F(0,0)*F(0,1)*F(1,1)*F(2,0) + d[21]*F(0,0)*F(0,1)*F(1,1)*F(2,0) + d[19]*pow(F(0,1),2)*F(1,1)*F(2,0) + d[4]*F(0,0)*F(0,2)*F(1,1)*F(2,0) + d[33]*F(0,0)*F(0,2)*F(1,1)*F(2,0) + d[22]*F(0,1)*F(0,2)*F(1,1)*F(2,0) + d[31]*F(0,1)*F(0,2)*F(1,1)*F(2,0) + d[34]*pow(F(0,2),2)*F(1,1)*F(2,0) + d[3]*F(0,0)*(F(0,1)*F(1,0) + F(0,0)*F(1,1))*F(2,0) + d[5]*pow(F(0,0),2)*F(1,2)*F(2,0) + d[4]*F(0,0)*F(0,1)*F(1,2)*F(2,0) + d[23]*F(0,0)*F(0,1)*F(1,2)*F(2,0) + d[22]*pow(F(0,1),2)*F(1,2)*F(2,0) + d[2]*F(0,0)*F(0,2)*F(1,2)*F(2,0) + d[35]*F(0,0)*F(0,2)*F(1,2)*F(2,0) + d[20]*F(0,1)*F(0,2)*F(1,2)*F(2,0) + d[34]*F(0,1)*F(0,2)*F(1,2)*F(2,0) + d[32]*pow(F(0,2),2)*F(1,2)*F(2,0) + d[18]*pow(F(0,0),2)*F(1,0)*F(2,1) + d[6]*F(0,0)*F(0,1)*F(1,0)*F(2,1) + d[21]*F(0,0)*F(0,1)*F(1,0)*F(2,1) + d[9]*pow(F(0,1),2)*F(1,0)*F(2,1) + d[23]*F(0,0)*F(0,2)*F(1,0)*F(2,1) + d[24]*F(0,0)*F(0,2)*F(1,0)*F(2,1) + d[11]*F(0,1)*F(0,2)*F(1,0)*F(2,1) + d[27]*F(0,1)*F(0,2)*F(1,0)*F(2,1) + d[29]*pow(F(0,2),2)*F(1,0)*F(2,1) + d[21]*pow(F(0,0),2)*F(1,1)*F(2,1) + d[9]*F(0,0)*F(0,1)*F(1,1)*F(2,1) + d[19]*F(0,0)*F(0,1)*F(1,1)*F(2,1) + d[7]*pow(F(0,1),2)*F(1,1)*F(2,1) + d[22]*F(0,0)*F(0,2)*F(1,1)*F(2,1) + d[27]*F(0,0)*F(0,2)*F(1,1)*F(2,1) + d[10]*F(0,1)*F(0,2)*F(1,1)*F(2,1) + d[25]*F(0,1)*F(0,2)*F(1,1)*F(2,1) + d[28]*pow(F(0,2),2)*F(1,1)*F(2,1) + d[23]*pow(F(0,0),2)*F(1,2)*F(2,1) + d[11]*F(0,0)*F(0,1)*F(1,2)*F(2,1) + d[22]*F(0,0)*F(0,1)*F(1,2)*F(2,1) + d[10]*pow(F(0,1),2)*F(1,2)*F(2,1) + d[20]*F(0,0)*F(0,2)*F(1,2)*F(2,1) + d[29]*F(0,0)*F(0,2)*F(1,2)*F(2,1) + d[8]*F(0,1)*F(0,2)*F(1,2)*F(2,1) + d[28]*F(0,1)*F(0,2)*F(1,2)*F(2,1) + d[26]*pow(F(0,2),2)*F(1,2)*F(2,1) + d[30]*pow(F(0,0),2)*F(1,0)*F(2,2) + d[24]*F(0,0)*F(0,1)*F(1,0)*F(2,2) + d[33]*F(0,0)*F(0,1)*F(1,0)*F(2,2) + d[27]*pow(F(0,1),2)*F(1,0)*F(2,2) + d[12]*F(0,0)*F(0,2)*F(1,0)*F(2,2) + d[35]*F(0,0)*F(0,2)*F(1,0)*F(2,2) + d[15]*F(0,1)*F(0,2)*F(1,0)*F(2,2) + d[29]*F(0,1)*F(0,2)*F(1,0)*F(2,2) + d[17]*pow(F(0,2),2)*F(1,0)*F(2,2) + d[33]*pow(F(0,0),2)*F(1,1)*F(2,2) + d[27]*F(0,0)*F(0,1)*F(1,1)*F(2,2) + d[31]*F(0,0)*F(0,1)*F(1,1)*F(2,2) + d[25]*pow(F(0,1),2)*F(1,1)*F(2,2) + d[15]*F(0,0)*F(0,2)*F(1,1)*F(2,2) + d[34]*F(0,0)*F(0,2)*F(1,1)*F(2,2) + d[13]*F(0,1)*F(0,2)*F(1,1)*F(2,2) + d[28]*F(0,1)*F(0,2)*F(1,1)*F(2,2) + d[16]*pow(F(0,2),2)*F(1,1)*F(2,2) + d[35]*pow(F(0,0),2)*F(1,2)*F(2,2) + d[29]*F(0,0)*F(0,1)*F(1,2)*F(2,2) + d[34]*F(0,0)*F(0,1)*F(1,2)*F(2,2) + d[28]*pow(F(0,1),2)*F(1,2)*F(2,2) + d[17]*F(0,0)*F(0,2)*F(1,2)*F(2,2) + d[32]*F(0,0)*F(0,2)*F(1,2)*F(2,2) + d[16]*F(0,1)*F(0,2)*F(1,2)*F(2,2) + d[26]*F(0,1)*F(0,2)*F(1,2)*F(2,2) + d[14]*pow(F(0,2),2)*F(1,2)*F(2,2);

    c.d[34] = d[0]*F(0,0)*F(1,0)*pow(F(2,0),2) + d[30]*F(0,2)*F(1,0)*pow(F(2,0),2) + d[3]*F(0,0)*F(1,1)*pow(F(2,0),2) + d[21]*F(0,1)*F(1,1)*pow(F(2,0),2) + d[33]*F(0,2)*F(1,1)*pow(F(2,0),2) + d[5]*F(0,0)*F(1,2)*pow(F(2,0),2) + d[23]*F(0,1)*F(1,2)*pow(F(2,0),2) + d[35]*F(0,2)*F(1,2)*pow(F(2,0),2) + d[3]*F(0,0)*F(1,0)*F(2,0)*F(2,1) + d[6]*F(0,1)*F(1,0)*F(2,0)*F(2,1) + d[21]*F(0,1)*F(1,0)*F(2,0)*F(2,1) + d[24]*F(0,2)*F(1,0)*F(2,0)*F(2,1) + d[33]*F(0,2)*F(1,0)*F(2,0)*F(2,1) + d[1]*F(0,0)*F(1,1)*F(2,0)*F(2,1) + d[21]*F(0,0)*F(1,1)*F(2,0)*F(2,1) + d[9]*F(0,1)*F(1,1)*F(2,0)*F(2,1) + d[19]*F(0,1)*F(1,1)*F(2,0)*F(2,1) + d[27]*F(0,2)*F(1,1)*F(2,0)*F(2,1) + d[31]*F(0,2)*F(1,1)*F(2,0)*F(2,1) + d[4]*F(0,0)*F(1,2)*F(2,0)*F(2,1) + d[23]*F(0,0)*F(1,2)*F(2,0)*F(2,1) + d[11]*F(0,1)*F(1,2)*F(2,0)*F(2,1) + d[22]*F(0,1)*F(1,2)*F(2,0)*F(2,1) + d[29]*F(0,2)*F(1,2)*F(2,0)*F(2,1) + d[34]*F(0,2)*F(1,2)*F(2,0)*F(2,1) + d[21]*F(0,0)*F(1,0)*pow(F(2,1),2) + d[9]*F(0,1)*F(1,0)*pow(F(2,1),2) + d[27]*F(0,2)*F(1,0)*pow(F(2,1),2) + d[19]*F(0,0)*F(1,1)*pow(F(2,1),2) + d[7]*F(0,1)*F(1,1)*pow(F(2,1),2) + d[25]*F(0,2)*F(1,1)*pow(F(2,1),2) + d[22]*F(0,0)*F(1,2)*pow(F(2,1),2) + d[10]*F(0,1)*F(1,2)*pow(F(2,1),2) + d[28]*F(0,2)*F(1,2)*pow(F(2,1),2) + d[18]*F(1,0)*F(2,0)*(F(0,1)*F(2,0) + F(0,0)*F(2,1)) + d[5]*F(0,0)*F(1,0)*F(2,0)*F(2,2) + d[30]*F(0,0)*F(1,0)*F(2,0)*F(2,2) + d[23]*F(0,1)*F(1,0)*F(2,0)*F(2,2) + d[24]*F(0,1)*F(1,0)*F(2,0)*F(2,2) + d[12]*F(0,2)*F(1,0)*F(2,0)*F(2,2) + d[35]*F(0,2)*F(1,0)*F(2,0)*F(2,2) + d[4]*F(0,0)*F(1,1)*F(2,0)*F(2,2) + d[33]*F(0,0)*F(1,1)*F(2,0)*F(2,2) + d[22]*F(0,1)*F(1,1)*F(2,0)*F(2,2) + d[27]*F(0,1)*F(1,1)*F(2,0)*F(2,2) + d[15]*F(0,2)*F(1,1)*F(2,0)*F(2,2) + d[34]*F(0,2)*F(1,1)*F(2,0)*F(2,2) + d[2]*F(0,0)*F(1,2)*F(2,0)*F(2,2) + d[35]*F(0,0)*F(1,2)*F(2,0)*F(2,2) + d[20]*F(0,1)*F(1,2)*F(2,0)*F(2,2) + d[29]*F(0,1)*F(1,2)*F(2,0)*F(2,2) + d[17]*F(0,2)*F(1,2)*F(2,0)*F(2,2) + d[32]*F(0,2)*F(1,2)*F(2,0)*F(2,2) + d[23]*F(0,0)*F(1,0)*F(2,1)*F(2,2) + d[33]*F(0,0)*F(1,0)*F(2,1)*F(2,2) + d[11]*F(0,1)*F(1,0)*F(2,1)*F(2,2) + d[27]*F(0,1)*F(1,0)*F(2,1)*F(2,2) + d[15]*F(0,2)*F(1,0)*F(2,1)*F(2,2) + d[29]*F(0,2)*F(1,0)*F(2,1)*F(2,2) + d[22]*F(0,0)*F(1,1)*F(2,1)*F(2,2) + d[31]*F(0,0)*F(1,1)*F(2,1)*F(2,2) + d[10]*F(0,1)*F(1,1)*F(2,1)*F(2,2) + d[25]*F(0,1)*F(1,1)*F(2,1)*F(2,2) + d[13]*F(0,2)*F(1,1)*F(2,1)*F(2,2) + d[28]*F(0,2)*F(1,1)*F(2,1)*F(2,2) + d[20]*F(0,0)*F(1,2)*F(2,1)*F(2,2) + d[34]*F(0,0)*F(1,2)*F(2,1)*F(2,2) + d[8]*F(0,1)*F(1,2)*F(2,1)*F(2,2) + d[28]*F(0,1)*F(1,2)*F(2,1)*F(2,2) + d[16]*F(0,2)*F(1,2)*F(2,1)*F(2,2) + d[26]*F(0,2)*F(1,2)*F(2,1)*F(2,2) + d[35]*F(0,0)*F(1,0)*pow(F(2,2),2) + d[29]*F(0,1)*F(1,0)*pow(F(2,2),2) + d[17]*F(0,2)*F(1,0)*pow(F(2,2),2) + d[34]*F(0,0)*F(1,1)*pow(F(2,2),2) + d[28]*F(0,1)*F(1,1)*pow(F(2,2),2) + d[16]*F(0,2)*F(1,1)*pow(F(2,2),2) + d[32]*F(0,0)*F(1,2)*pow(F(2,2),2) + d[26]*F(0,1)*F(1,2)*pow(F(2,2),2) + d[14]*F(0,2)*F(1,2)*pow(F(2,2),2);

    c.d[35] = d[0]*pow(F(0,0),2)*pow(F(2,0),2) + d[18]*F(0,0)*F(0,1)*pow(F(2,0),2) + d[21]*pow(F(0,1),2)*pow(F(2,0),2) + d[5]*F(0,0)*F(0,2)*pow(F(2,0),2) + d[30]*F(0,0)*F(0,2)*pow(F(2,0),2) + d[23]*F(0,1)*F(0,2)*pow(F(2,0),2) + d[33]*F(0,1)*F(0,2)*pow(F(2,0),2) + d[35]*pow(F(0,2),2)*pow(F(2,0),2) + d[18]*pow(F(0,0),2)*F(2,0)*F(2,1) + d[1]*F(0,0)*F(0,1)*F(2,0)*F(2,1) + d[6]*F(0,0)*F(0,1)*F(2,0)*F(2,1) + 2*d[21]*F(0,0)*F(0,1)*F(2,0)*F(2,1) + d[9]*pow(F(0,1),2)*F(2,0)*F(2,1) + d[19]*pow(F(0,1),2)*F(2,0)*F(2,1) + d[4]*F(0,0)*F(0,2)*F(2,0)*F(2,1) + d[23]*F(0,0)*F(0,2)*F(2,0)*F(2,1) + d[24]*F(0,0)*F(0,2)*F(2,0)*F(2,1) + d[33]*F(0,0)*F(0,2)*F(2,0)*F(2,1) + d[11]*F(0,1)*F(0,2)*F(2,0)*F(2,1) + d[22]*F(0,1)*F(0,2)*F(2,0)*F(2,1) + d[27]*F(0,1)*F(0,2)*F(2,0)*F(2,1) + d[31]*F(0,1)*F(0,2)*F(2,0)*F(2,1) + d[29]*pow(F(0,2),2)*F(2,0)*F(2,1) + d[34]*pow(F(0,2),2)*F(2,0)*F(2,1) + d[21]*pow(F(0,0),2)*pow(F(2,1),2) + d[9]*F(0,0)*F(0,1)*pow(F(2,1),2) + d[19]*F(0,0)*F(0,1)*pow(F(2,1),2) + d[7]*pow(F(0,1),2)*pow(F(2,1),2) + d[22]*F(0,0)*F(0,2)*pow(F(2,1),2) + d[27]*F(0,0)*F(0,2)*pow(F(2,1),2) + d[10]*F(0,1)*F(0,2)*pow(F(2,1),2) + d[25]*F(0,1)*F(0,2)*pow(F(2,1),2) + d[28]*pow(F(0,2),2)*pow(F(2,1),2) + d[3]*F(0,0)*F(2,0)*(F(0,1)*F(2,0) + F(0,0)*F(2,1)) + d[5]*pow(F(0,0),2)*F(2,0)*F(2,2) + d[30]*pow(F(0,0),2)*F(2,0)*F(2,2) + d[4]*F(0,0)*F(0,1)*F(2,0)*F(2,2) + d[23]*F(0,0)*F(0,1)*F(2,0)*F(2,2) + d[24]*F(0,0)*F(0,1)*F(2,0)*F(2,2) + d[33]*F(0,0)*F(0,1)*F(2,0)*F(2,2) + d[22]*pow(F(0,1),2)*F(2,0)*F(2,2) + d[27]*pow(F(0,1),2)*F(2,0)*F(2,2) + d[2]*F(0,0)*F(0,2)*F(2,0)*F(2,2) + d[12]*F(0,0)*F(0,2)*F(2,0)*F(2,2) + 2*d[35]*F(0,0)*F(0,2)*F(2,0)*F(2,2) + d[15]*F(0,1)*F(0,2)*F(2,0)*F(2,2) + d[20]*F(0,1)*F(0,2)*F(2,0)*F(2,2) + d[29]*F(0,1)*F(0,2)*F(2,0)*F(2,2) + d[34]*F(0,1)*F(0,2)*F(2,0)*F(2,2) + d[17]*pow(F(0,2),2)*F(2,0)*F(2,2) + d[32]*pow(F(0,2),2)*F(2,0)*F(2,2) + d[23]*pow(F(0,0),2)*F(2,1)*F(2,2) + d[33]*pow(F(0,0),2)*F(2,1)*F(2,2) + d[11]*F(0,0)*F(0,1)*F(2,1)*F(2,2) + d[22]*F(0,0)*F(0,1)*F(2,1)*F(2,2) + d[27]*F(0,0)*F(0,1)*F(2,1)*F(2,2) + d[31]*F(0,0)*F(0,1)*F(2,1)*F(2,2) + d[10]*pow(F(0,1),2)*F(2,1)*F(2,2) + d[25]*pow(F(0,1),2)*F(2,1)*F(2,2) + d[15]*F(0,0)*F(0,2)*F(2,1)*F(2,2) + d[20]*F(0,0)*F(0,2)*F(2,1)*F(2,2) + d[29]*F(0,0)*F(0,2)*F(2,1)*F(2,2) + d[34]*F(0,0)*F(0,2)*F(2,1)*F(2,2) + d[8]*F(0,1)*F(0,2)*F(2,1)*F(2,2) + d[13]*F(0,1)*F(0,2)*F(2,1)*F(2,2) + 2*d[28]*F(0,1)*F(0,2)*F(2,1)*F(2,2) + d[16]*pow(F(0,2),2)*F(2,1)*F(2,2) + d[26]*pow(F(0,2),2)*F(2,1)*F(2,2) + d[35]*pow(F(0,0),2)*pow(F(2,2),2) + d[29]*F(0,0)*F(0,1)*pow(F(2,2),2) + d[34]*F(0,0)*F(0,1)*pow(F(2,2),2) + d[28]*pow(F(0,1),2)*pow(F(2,2),2) + d[17]*F(0,0)*F(0,2)*pow(F(2,2),2) + d[32]*F(0,0)*F(0,2)*pow(F(2,2),2) + d[16]*F(0,1)*F(0,2)*pow(F(2,2),2) + d[26]*F(0,1)*F(0,2)*pow(F(2,2),2) + d[14]*pow(F(0,2),2)*pow(F(2,2),2);

    return c;
}
