// NOTE: This file is automatically included from tens4d.h
// Users should not include this file manually!

#include "matrix.h"

inline tens4d::tens4d(const double g)
{
	for (int i = 0; i < NNZ; i++)
		d[i] = g;
}

inline tens4d::tens4d(double m[9][9])
{
	d[ 0] = m[0][0]; d[ 1] = m[1][0]; d[ 2] = m[2][0]; d[ 3] = m[3][0]; d[ 4] = m[4][0]; d[ 5] = m[5][0]; d[ 6] = m[6][0]; d[ 7] = m[7][0]; d[ 8] = m[8][0];
	d[ 9] = m[0][1]; d[10] = m[1][1]; d[11] = m[2][1]; d[12] = m[3][1]; d[13] = m[4][1]; d[14] = m[5][1]; d[15] = m[6][1]; d[16] = m[7][1]; d[17] = m[8][1];
	d[18] = m[0][2]; d[19] = m[1][2]; d[20] = m[2][2]; d[21] = m[3][2]; d[22] = m[4][2]; d[23] = m[5][2]; d[24] = m[6][2]; d[25] = m[7][2]; d[26] = m[8][2];
	d[27] = m[0][3]; d[28] = m[1][3]; d[29] = m[2][3]; d[30] = m[3][3]; d[31] = m[4][3]; d[32] = m[5][3]; d[33] = m[6][3]; d[34] = m[7][3]; d[35] = m[8][3];
	d[36] = m[0][4]; d[37] = m[1][4]; d[38] = m[2][4]; d[39] = m[3][4]; d[40] = m[4][4]; d[41] = m[5][4]; d[42] = m[6][4]; d[43] = m[7][4]; d[44] = m[8][4];
	d[45] = m[0][5]; d[46] = m[1][5]; d[47] = m[2][5]; d[48] = m[3][5]; d[49] = m[4][5]; d[50] = m[5][5]; d[51] = m[6][5]; d[52] = m[7][5]; d[53] = m[8][5];
	d[54] = m[0][6]; d[55] = m[1][6]; d[56] = m[2][6]; d[57] = m[3][6]; d[58] = m[4][6]; d[59] = m[5][6]; d[60] = m[6][6]; d[61] = m[7][6]; d[62] = m[8][6];
	d[63] = m[0][7]; d[64] = m[1][7]; d[65] = m[2][7]; d[66] = m[3][7]; d[67] = m[4][7]; d[68] = m[5][7]; d[69] = m[6][7]; d[70] = m[7][7]; d[71] = m[8][7];
	d[72] = m[0][8]; d[73] = m[1][8]; d[74] = m[2][8]; d[75] = m[3][8]; d[76] = m[4][8]; d[77] = m[5][8]; d[78] = m[6][8]; d[79] = m[7][8]; d[80] = m[8][8];
}

inline double& tens4d::operator () (int i, int j, int k, int l)
{
	const int m[3][3] = {{0,3,5},{6,1,4},{8,7,2}};
	tens4d& T = (*this);
	return T(m[i][j], m[k][l]);
}

inline double tens4d::operator () (int i, int j, int k, int l) const
{
	const int m[3][3] = {{0,3,5},{6,1,4},{8,7,2}};
	const tens4d& T = (*this);
	return T(m[i][j], m[k][l]);
}

inline double& tens4d::operator () (int i, int j)
{
	const int m[9] = {0, 9, 19, 27, 36, 45, 54, 63, 72};
	return d[m[j]+i];
}

inline double tens4d::operator () (int i, int j) const
{
	const int m[9] = {0, 9, 19, 27, 36, 45, 54, 63, 72};
	return d[m[j]+i];
}

// operator +
inline tens4d tens4d::operator + (const tens4d& t) const
{
	tens4d s;
	for (int i=0; i<NNZ; i++)
		s.d[i] = d[i] + t.d[i];
	
	return s;
}

// operator -
inline tens4d tens4d::operator - (const tens4d& t) const
{
	tens4d s;
	for (int i=0; i<NNZ; i++)
		s.d[i] = d[i] - t.d[i];

	return s;
}

// operator *
inline tens4d tens4d::operator * (double g) const
{
	tens4d s;
	for (int i=0; i<NNZ; i++)
		s.d[i] = g*d[i];
	
	return s;
}

// operator /
inline tens4d tens4d::operator / (double g) const
{
	tens4d s;
	for (int i=0; i<NNZ; i++)
		s.d[i] = d[i]/g;
	
	return s;
}

// assignment operator +=
inline tens4d& tens4d::operator += (const tens4d& t)
{
	for (int i=0; i<NNZ; i++)
		d[i] += t.d[i];
	
	return (*this);
}

// assignment operator -=
inline tens4d& tens4d::operator -= (const tens4d& t)
{
	for (int i=0; i<NNZ; i++)
		d[i] -= t.d[i];
	
	return (*this);
}

// assignment operator *=
inline tens4d& tens4d::operator *= (double g)
{
	for (int i=0; i<NNZ; i++)
		d[i] *= g;
	
	return (*this);
}

// assignment operator /=
inline tens4d& tens4d::operator /= (double g)
{
	for (int i=0; i<NNZ; i++)
		d[i] /= g;
	
	return (*this);
}

// unary operator -
inline tens4d tens4d::operator - () const
{
	tens4d s;
	for (int i = 0; i < NNZ; i++)
		s.d[i] = -d[i];

	return s;
}

// trace
// C.tr() = I:C:I

inline double tens4d::tr() const
{
	return (d[0] + d[9] + d[18] + d[1] + d[10] + d[19] + d[2] + d[11] + d[20]);
}

// intialize to zero
inline void tens4d::zero()
{
	for (int i = 0; i < NNZ; i++)
		d[i] = 0;
}

// extract 9x9 matrix
inline void tens4d::extract(double D[9][9])
{
	D[0][0] = d[0]; D[0][1] = d[9];  D[0][2] = d[18]; D[0][3] = d[27]; D[0][4] = d[36]; D[0][5] = d[45]; D[0][6] = d[54]; D[0][7] = d[63]; D[0][8] = d[72];
	D[1][0] = d[1]; D[1][1] = d[10]; D[1][2] = d[19]; D[1][3] = d[28]; D[1][4] = d[37]; D[1][5] = d[46]; D[1][6] = d[55]; D[1][7] = d[64]; D[1][8] = d[73];
	D[2][0] = d[2]; D[2][1] = d[11]; D[2][2] = d[20]; D[2][3] = d[29]; D[2][4] = d[38]; D[2][5] = d[47]; D[2][6] = d[56]; D[2][7] = d[65]; D[2][8] = d[74];
	D[3][0] = d[3]; D[3][1] = d[12]; D[3][2] = d[21]; D[3][3] = d[30]; D[3][4] = d[39]; D[3][5] = d[48]; D[3][6] = d[57]; D[3][7] = d[66]; D[3][8] = d[75];
	D[4][0] = d[4]; D[4][1] = d[13]; D[4][2] = d[22]; D[4][3] = d[31]; D[4][4] = d[40]; D[4][5] = d[49]; D[4][6] = d[58]; D[4][7] = d[67]; D[4][8] = d[76];
	D[5][0] = d[5]; D[5][1] = d[14]; D[5][2] = d[23]; D[5][3] = d[32]; D[5][4] = d[41]; D[5][5] = d[50]; D[5][6] = d[59]; D[5][7] = d[68]; D[5][8] = d[77];
	D[6][0] = d[6]; D[6][1] = d[15]; D[6][2] = d[24]; D[6][3] = d[33]; D[6][4] = d[42]; D[6][5] = d[51]; D[6][6] = d[60]; D[6][7] = d[69]; D[6][8] = d[78];
	D[7][0] = d[7]; D[7][1] = d[16]; D[7][2] = d[25]; D[7][3] = d[34]; D[7][4] = d[43]; D[7][5] = d[52]; D[7][6] = d[61]; D[7][7] = d[70]; D[7][8] = d[79];
	D[8][0] = d[8]; D[8][1] = d[17]; D[8][2] = d[26]; D[8][3] = d[35]; D[8][4] = d[44]; D[8][5] = d[53]; D[8][6] = d[62]; D[8][7] = d[71]; D[8][8] = d[80];
}

// compute the super symmetric (major and minor symmetric) component of the tensor
// Sijkl = (1/8)*(Cijkl + Cijlk + Cjikl + Cjilk + Cklij + Clkij + Cklji + Clkji)
inline tens4ds tens4d::supersymm() const
{
	tens4ds s;

	s.d[0] = d[0]; s.d[1] = (d[9] + d[1])/2.; s.d[3] = (d[18] + d[2])/2.;  s.d[6] = (d[27] + d[54] + d[3] + d[6])/4.;   s.d[10] = (d[36] + d[63] + d[4] + d[7])/4.;                                   s.d[15] = (d[45] + d[72] + d[5] + d[8])/4.;
	               s.d[2] = d[10];            s.d[4] = (d[19] + d[11])/2.; s.d[7] = (d[28] + d[55] + d[12] + d[15])/4.; s.d[11] = (d[37] + d[64] + d[13] + d[16])/4.;                                 s.d[16] = (d[46] + d[73] + d[14] + d[17])/4.;
				                              s.d[5] = d[20];              s.d[8] = (d[29] + d[56] + d[21] + d[24])/4.; s.d[12] = (d[38] + d[65] + d[22] + d[25])/4.;                                 s.d[17] = (d[47] + d[74] + d[23] + d[26])/4.;
								                                           s.d[9] = (d[30] + d[57] + d[33] + d[60])/4.; s.d[13] = (d[39] + d[66] + d[42] + d[69] + d[31] + d[34] + d[58] + d[61])/8.; s.d[18] = (d[48] + d[75] + d[51] + d[78] + d[32] + d[35] + d[59] + d[62])/8.;
												                                                                        s.d[14] = (d[40] + d[67] + d[43] + d[70])/4.;                                 s.d[19] = (d[49] + d[76] + d[52] + d[79] + d[41] + d[44] + d[68] + d[71])/4.;
																			                                                                                                                          s.d[20] = (d[50] + d[77] + d[53] + d[80])/4.;

	return s;
}
