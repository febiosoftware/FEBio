// NOTE: This file is automatically included from tens4d.h
// Users should not include this file manually!

#include "matrix.h"

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
inline tens4dms dyad1(const mat3d& a)
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
inline tens4dms dyad1(const mat3d& a, const mat3d& b)
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

//-----------------------------------------------------------------------------
// (a dyad2s a)_ijkl = a_ik a_jl
// The resulting tensor only has major symmetry if a is symmetric
//inline tens4dms dyad2s(const mat3ds& a)
//{
//	tens4dms c;
//	
//	c.d[ 0] = a.xx()*a.xx();
//	c.d[ 1] = a.xy()*a.xy();
//	c.d[ 3] = a.xz()*a.xz();
//	c.d[ 6] = a.xx()*a.xy();
//	c.d[10] = a.xy()*a.xz();
//	c.d[15] = a.xx()*a.xz();
//	c.d[21] = a.xy()*a.xx();
//	c.d[29] = a.xz()*a.xy();
//	c.d[36] = a.xz()*a.xx();
//
//	c.d[ 2] = a.yy()*a.yy();
//	c.d[ 4] = a.yz()*a.yz();
//	c.d[ 7] = a.xy()*a.yy();
//	c.d[11] = a.yy()*a.yz();
//	c.d[16] = a.xy()*a.yz();
//	c.d[22] = a.yy()*a.xy();
//	c.d[29] = a.yz()*a.yy();
//	c.d[37] = a.yz()*a.xy();
//
//	c.d[ 5] = a.zz()*a.zz();
//	c.d[ 8] = a.xz()*a.yz();
//	c.d[12] = a.yz()*a.zz();
//	c.d[17] = a.xz()*a.zz();
//	c.d[23] = a.yz()*a.xz();
//	c.d[30] = a.zz()*a.yz();
//	c.d[38] = a.zz()*a.xz();
//
//	c.d[ 9] = a.xx()*a.yy();
//	c.d[13] = a.xy()*a.yz();
//	c.d[18] = a.xx()*a.yz();
//	c.d[24] = a.xy()*a.xy();
//	c.d[31] = a.xz()*a.yy();
//	c.d[39] = a.xz()*a.xy();
//
//	c.d[14] = a.yy()*a.zz();
//	c.d[19] = a.xy()*a.zz();
//	c.d[25] = a.yy()*a.yz();
//	c.d[32] = a.yz()*a.yz();
//	c.d[40] = a.yz()*a.xz();
//
//	c.d[20] = a.xx()*a.zz();
//	c.d[26] = a.xy()*a.xz();
//	c.d[33] = a.xz()*a.yz();
//	c.d[41] = a.xz()*a.xz();
//
//	c.d[27] = a.xy()*a.xy();
//	c.d[34] = a.yz()*a.xy();
//	c.d[42] = a.yz()*a.xx();
//
//	c.d[35] = a.zz()*a.yy();
//	c.d[43] = a.zz()*a.xy();
//
//	c.d[44] = a.zz()*a.xx();
//
//	return c;
//}


//-----------------------------------------------------------------------------
// (a dyad2s b)_ijkl = a_ik b_jl + b_ik a_jl
// The resulting tensor only has major symmetry if a,b are symmetric
//inline tens4dms dyad2s(const mat3ds& a, const mat3ds& b)
//{
//	tens4dms c;
//	
//	c.d[ 0] = a.xx()*b.xx() + b.xx()*a.xx();
//	c.d[ 1] = a.xy()*b.xy() + b.xy()*a.xy();
//	c.d[ 3] = a.xz()*b.xz() + b.xz()*a.xz();
//	c.d[ 6] = a.xx()*b.xy() + b.xx()*a.xy();
//	c.d[10] = a.xy()*b.xz() + b.xy()*a.xz();
//	c.d[15] = a.xx()*b.xz() + b.xx()*a.xz();
//	c.d[21] = a.xy()*b.xx() + b.xy()*a.xx();
//	c.d[29] = a.xz()*b.xy() + b.xz()*a.xy();
//	c.d[36] = a.xz()*b.xx() + b.xz()*a.xx();
//	
//	c.d[ 2] = a.yy()*b.yy() + b.yy()*a.yy();
//	c.d[ 4] = a.yz()*b.yz() + b.yz()*a.yz();
//	c.d[ 7] = a.xy()*b.yy() + b.xy()*a.yy();
//	c.d[11] = a.yy()*b.yz() + b.yy()*a.yz();
//	c.d[16] = a.xy()*b.yz() + b.xy()*a.yz();
//	c.d[22] = a.yy()*b.xy() + b.yy()*a.xy();
//	c.d[29] = a.yz()*b.yy() + b.yz()*a.yy();
//	c.d[37] = a.yz()*b.xy() + b.yz()*a.xy();
//
//	c.d[ 5] = a.zz()*b.zz() + b.zz()*a.zz();
//	c.d[ 8] = a.xz()*b.yz() + b.xz()*a.yz();
//	c.d[12] = a.yz()*b.zz() + b.yz()*a.zz();
//	c.d[17] = a.xz()*b.zz() + b.xz()*a.zz();
//	c.d[23] = a.yz()*b.xz() + b.yz()*a.xz();
//	c.d[30] = a.zz()*b.yz() + b.zz()*a.yz();
//	c.d[38] = a.zz()*b.xz() + b.zz()*a.xz();
//	
//	c.d[ 9] = a.xx()*b.yy() + b.xx()*a.yy();
//	c.d[13] = a.xy()*b.yz() + b.xy()*a.yz();
//	c.d[18] = a.xx()*b.yz() + b.xx()*a.yz();
//	c.d[24] = a.xy()*b.xy() + b.xy()*a.xy();
//	c.d[31] = a.xz()*b.yy() + b.xz()*a.yy();
//	c.d[39] = a.xz()*b.xy() + b.xz()*a.xy();
//
//	c.d[14] = a.yy()*b.zz() + b.yy()*a.zz();
//	c.d[19] = a.xy()*b.zz() + b.xy()*a.zz();
//	c.d[25] = a.yy()*b.yz() + b.yy()*a.yz();
//	c.d[32] = a.yz()*b.yz() + b.yz()*a.yz();
//	c.d[40] = a.yz()*b.xz() + b.yz()*a.xz();
//
//	c.d[20] = a.xx()*b.zz() + b.xx()*a.zz();
//	c.d[26] = a.xy()*b.xz() + b.xy()*a.xz();
//	c.d[33] = a.xz()*b.yz() + b.xz()*a.yz();
//	c.d[41] = a.xz()*b.xz() + b.xz()*a.xz();
//
//	c.d[27] = a.xy()*b.xy() + b.xy()*a.xy();
//	c.d[34] = a.yz()*b.xy() + b.yz()*a.xy();
//	c.d[42] = a.yz()*b.xx() + b.yz()*a.xx();
//
//	c.d[35] = a.zz()*b.yy() + b.zz()*a.yy();
//	c.d[43] = a.zz()*b.xy() + b.zz()*a.xy();
//
//	c.d[44] = a.zz()*b.xx() + b.zz()*a.xx();
//
//	return c;
//}

////-----------------------------------------------------------------------------
//// (a ddots b)_ijkl = a_ijmn b_mnkl + b_ijmn a_mnkl
//inline tens4ds ddots(const tens4ds& a, const tens4ds& b)
//{
//	tens4ds c;
//	
//	c.d[0] = 2*(a.d[0]*b.d[0] + a.d[1]*b.d[1] + a.d[3]*b.d[3] + 2*a.d[6]*b.d[6] 
//				+ 2*a.d[10]*b.d[10] + 2*a.d[15]*b.d[15]);
//	c.d[1] = a.d[0]*b.d[1] + a.d[2]*b.d[1] + a.d[1]*(b.d[0] + b.d[2]) + a.d[4]*b.d[3] 
//	+ a.d[3]*b.d[4] + 2*a.d[7]*b.d[6] + 2*a.d[6]*b.d[7] + 2*a.d[11]*b.d[10] 
//	+ 2*a.d[10]*b.d[11] + 2*a.d[16]*b.d[15] + 2*a.d[15]*b.d[16];
//	c.d[3] = a.d[4]*b.d[1] + a.d[0]*b.d[3] + a.d[5]*b.d[3] + a.d[1]*b.d[4] 
//	+ a.d[3]*(b.d[0] + b.d[5]) + 2*a.d[8]*b.d[6] + 2*a.d[6]*b.d[8] + 2*a.d[12]*b.d[10] 
//	+ 2*a.d[10]*b.d[12] + 2*a.d[17]*b.d[15] + 2*a.d[15]*b.d[17];
//	c.d[6] = a.d[7]*b.d[1] + a.d[8]*b.d[3] + a.d[0]*b.d[6] + 2*a.d[9]*b.d[6] + a.d[1]*b.d[7] 
//	+ a.d[3]*b.d[8] + a.d[6]*(b.d[0] + 2*b.d[9]) + 2*a.d[13]*b.d[10] + 2*a.d[10]*b.d[13] 
//	+ 2*a.d[18]*b.d[15] + 2*a.d[15]*b.d[18];
//	c.d[10] = a.d[11]*b.d[1] + a.d[12]*b.d[3] + 2*a.d[13]*b.d[6] + a.d[0]*b.d[10] 
//	+ 2*a.d[14]*b.d[10] + a.d[1]*b.d[11] + a.d[3]*b.d[12] + 2*a.d[6]*b.d[13] 
//	+ a.d[10]*(b.d[0] + 2*b.d[14]) + 2*a.d[19]*b.d[15] + 2*a.d[15]*b.d[19];
//	c.d[15] = a.d[16]*b.d[1] + a.d[17]*b.d[3] + 2*a.d[18]*b.d[6] + 2*a.d[19]*b.d[10] 
//	+ a.d[0]*b.d[15] + 2*a.d[20]*b.d[15] + a.d[1]*b.d[16] + a.d[3]*b.d[17] 
//	+ 2*a.d[6]*b.d[18] + 2*a.d[10]*b.d[19] + a.d[15]*(b.d[0] + 2*b.d[20]);
//	
//	c.d[2] = 2*(a.d[1]*b.d[1] + a.d[2]*b.d[2] + a.d[4]*b.d[4] + 2*a.d[7]*b.d[7] 
//				+ 2*a.d[11]*b.d[11] + 2*a.d[16]*b.d[16]);
//	c.d[4] = a.d[3]*b.d[1] + a.d[1]*b.d[3] + a.d[2]*b.d[4] + a.d[5]*b.d[4] 
//	+ a.d[4]*(b.d[2] + b.d[5]) + 2*a.d[8]*b.d[7] + 2*a.d[7]*b.d[8] + 2*a.d[12]*b.d[11] 
//	+ 2*a.d[11]*b.d[12] + 2*a.d[17]*b.d[16] + 2*a.d[16]*b.d[17];
//	c.d[7] = a.d[6]*b.d[1] + a.d[8]*b.d[4] + a.d[1]*b.d[6] + a.d[2]*b.d[7] 
//	+ 2*a.d[9]*b.d[7] + a.d[4]*b.d[8] + a.d[7]*(b.d[2] + 2*b.d[9]) + 2*a.d[13]*b.d[11] 
//	+ 2*a.d[11]*b.d[13] + 2*a.d[18]*b.d[16] + 2*a.d[16]*b.d[18];
//	c.d[11] = a.d[10]*b.d[1] + a.d[12]*b.d[4] + 2*a.d[13]*b.d[7] + a.d[1]*b.d[10] 
//	+ a.d[2]*b.d[11] + 2*a.d[14]*b.d[11] + a.d[4]*b.d[12] + 2*a.d[7]*b.d[13] 
//	+ a.d[11]*(b.d[2] + 2*b.d[14]) + 2*a.d[19]*b.d[16] + 2*a.d[16]*b.d[19];
//	c.d[16] = a.d[15]*b.d[1] + a.d[17]*b.d[4] + 2*a.d[18]*b.d[7] + 2*a.d[19]*b.d[11] 
//	+ a.d[1]*b.d[15] + a.d[2]*b.d[16] + 2*a.d[20]*b.d[16] + a.d[4]*b.d[17] 
//	+ 2*a.d[7]*b.d[18] + 2*a.d[11]*b.d[19] + a.d[16]*(b.d[2] + 2*b.d[20]);
//	
//	c.d[5] = 2*(a.d[3]*b.d[3] + a.d[4]*b.d[4] + a.d[5]*b.d[5] + 2*a.d[8]*b.d[8] 
//				+ 2*a.d[12]*b.d[12] + 2*a.d[17]*b.d[17]);
//	c.d[8] = a.d[6]*b.d[3] + a.d[7]*b.d[4] + a.d[8]*b.d[5] + a.d[3]*b.d[6] 
//	+ a.d[4]*b.d[7] + a.d[5]*b.d[8] + 2*a.d[9]*b.d[8] + 2*a.d[8]*b.d[9] + 2*a.d[13]*b.d[12] 
//	+ 2*a.d[12]*b.d[13] + 2*a.d[18]*b.d[17] + 2*a.d[17]*b.d[18];
//	c.d[12] = a.d[10]*b.d[3] + a.d[11]*b.d[4] + a.d[12]*b.d[5] + 2*a.d[13]*b.d[8] 
//	+ a.d[3]*b.d[10] + a.d[4]*b.d[11] + a.d[5]*b.d[12] + 2*a.d[14]*b.d[12] + 2*a.d[8]*b.d[13] 
//	+ 2*a.d[12]*b.d[14] + 2*a.d[19]*b.d[17] + 2*a.d[17]*b.d[19];
//	c.d[17] = a.d[15]*b.d[3] + a.d[16]*b.d[4] + a.d[17]*b.d[5] + 2*a.d[18]*b.d[8] 
//	+ 2*a.d[19]*b.d[12] + a.d[3]*b.d[15] + a.d[4]*b.d[16] + a.d[5]*b.d[17] 
//	+ 2*a.d[20]*b.d[17] + 2*a.d[8]*b.d[18] + 2*a.d[12]*b.d[19] + 2*a.d[17]*b.d[20];
//	
//	c.d[9] = 2*(a.d[6]*b.d[6] + a.d[7]*b.d[7] + a.d[8]*b.d[8] + 2*a.d[9]*b.d[9] 
//				+ 2*a.d[13]*b.d[13] + 2*a.d[18]*b.d[18]);
//	c.d[13] = a.d[10]*b.d[6] + a.d[11]*b.d[7] + a.d[12]*b.d[8] + 2*a.d[13]*b.d[9] 
//	+ a.d[6]*b.d[10] + a.d[7]*b.d[11] + a.d[8]*b.d[12] + 2*a.d[9]*b.d[13] 
//	+ 2*a.d[14]*b.d[13] + 2*a.d[13]*b.d[14] + 2*a.d[19]*b.d[18] + 2*a.d[18]*b.d[19];
//	c.d[18] = a.d[15]*b.d[6] + a.d[16]*b.d[7] + a.d[17]*b.d[8] + 2*a.d[18]*b.d[9] 
//	+ 2*a.d[19]*b.d[13] + a.d[6]*b.d[15] + a.d[7]*b.d[16] + a.d[8]*b.d[17] 
//	+ 2*a.d[9]*b.d[18] + 2*a.d[20]*b.d[18] + 2*a.d[13]*b.d[19] + 2*a.d[18]*b.d[20];
//	
//	c.d[14] = 2*(a.d[10]*b.d[10] + a.d[11]*b.d[11] + a.d[12]*b.d[12] 
//				 + 2*a.d[13]*b.d[13] + 2*a.d[14]*b.d[14] + 2*a.d[19]*b.d[19]);
//	c.d[19] = a.d[15]*b.d[10] + a.d[16]*b.d[11] + a.d[17]*b.d[12] 
//	+ 2*a.d[18]*b.d[13] + 2*a.d[19]*b.d[14] + a.d[10]*b.d[15] + a.d[11]*b.d[16] 
//	+ a.d[12]*b.d[17] + 2*a.d[13]*b.d[18] + 2*a.d[14]*b.d[19] + 2*a.d[20]*b.d[19] + 2*a.d[19]*b.d[20];
//	
//	c.d[20] = 2*(a.d[15]*b.d[15] + a.d[16]*b.d[16] + a.d[17]*b.d[17] 
//				 + 2*a.d[18]*b.d[18] + 2*a.d[19]*b.d[19] + 2*a.d[20]*b.d[20]);
//	
//	return c;
//	
//}
//
////-----------------------------------------------------------------------------
//// double contraction with a symmetric 2nd-order tensor
//inline mat3ds tens4ds::dot(const mat3ds &m) const
//{
//	mat3ds a;
//	a.xx() = d[ 0]*m.xx() + d[ 1]*m.yy() + d[ 3]*m.zz() + 2*d[ 6]*m.xy() + 2*d[10]*m.yz() + 2*d[15]*m.xz();
//	a.yy() = d[ 1]*m.xx() + d[ 2]*m.yy() + d[ 4]*m.zz() + 2*d[ 7]*m.xy() + 2*d[11]*m.yz() + 2*d[16]*m.xz();
//	a.zz() = d[ 3]*m.xx() + d[ 4]*m.yy() + d[ 5]*m.zz() + 2*d[ 8]*m.xy() + 2*d[12]*m.yz() + 2*d[17]*m.xz();
//	a.xy() = d[ 6]*m.xx() + d[ 7]*m.yy() + d[ 8]*m.zz() + 2*d[ 9]*m.xy() + 2*d[13]*m.yz() + 2*d[18]*m.xz();
//	a.yz() = d[10]*m.xx() + d[11]*m.yy() + d[12]*m.zz() + 2*d[13]*m.xy() + 2*d[14]*m.yz() + 2*d[19]*m.xz();
//	a.xz() = d[15]*m.xx() + d[16]*m.yy() + d[17]*m.zz() + 2*d[18]*m.xy() + 2*d[19]*m.yz() + 2*d[20]*m.xz();
//	return a;
//}
//
////-----------------------------------------------------------------------------
//// vdotTdotv_jk = a_i T_ijkl b_l
//inline mat3d vdotTdotv(const vec3d a, const tens4ds T, const vec3d b)
//{
//	return mat3d(a.x*(b.x*T.d[0] + b.y*T.d[6] + b.z*T.d[15]) + a.y*(b.x*T.d[6] + b.y*T.d[9] + b.z*T.d[18]) + a.z*(b.x*T.d[15] + b.y*T.d[18] + b.z*T.d[20]),
//				 a.x*(b.y*T.d[1] + b.x*T.d[6] + b.z*T.d[10]) + a.y*(b.y*T.d[7] + b.x*T.d[9] + b.z*T.d[13]) + a.z*(b.y*T.d[16] + b.x*T.d[18] + b.z*T.d[19]),
//				 a.x*(b.z*T.d[3] + b.y*T.d[10] + b.x*T.d[15]) + a.y*(b.z*T.d[8] + b.y*T.d[13] + b.x*T.d[18]) + a.z*(b.z*T.d[17] + b.y*T.d[19] + b.x*T.d[20]),
//				 a.y*(b.x*T.d[1] + b.y*T.d[7] + b.z*T.d[16]) + a.x*(b.x*T.d[6] + b.y*T.d[9] + b.z*T.d[18]) + a.z*(b.x*T.d[10] + b.y*T.d[13] + b.z*T.d[19]),
//				 a.y*(b.y*T.d[2] + b.x*T.d[7] + b.z*T.d[11]) + a.x*(b.y*T.d[7] + b.x*T.d[9] + b.z*T.d[13]) + a.z*(b.y*T.d[11] + b.x*T.d[13] + b.z*T.d[14]),
//				 a.y*(b.z*T.d[4] + b.y*T.d[11] + b.x*T.d[16]) + a.x*(b.z*T.d[8] + b.y*T.d[13] + b.x*T.d[18]) + a.z*(b.z*T.d[12] + b.y*T.d[14] + b.x*T.d[19]),
//				 a.z*(b.x*T.d[3] + b.y*T.d[8] + b.z*T.d[17]) + a.y*(b.x*T.d[10] + b.y*T.d[13] + b.z*T.d[19]) + a.x*(b.x*T.d[15] + b.y*T.d[18] + b.z*T.d[20]),
//				 a.z*(b.y*T.d[4] + b.x*T.d[8] + b.z*T.d[12]) + a.y*(b.y*T.d[11] + b.x*T.d[13] + b.z*T.d[14]) + a.x*(b.y*T.d[16] + b.x*T.d[18] + b.z*T.d[19]),
//				 a.z*(b.z*T.d[5] + b.y*T.d[12] + b.x*T.d[17]) + a.y*(b.z*T.d[12] + b.y*T.d[14] + b.x*T.d[19]) + a.x*(b.z*T.d[17] + b.y*T.d[19] + b.x*T.d[20]));
//}
//
//
////-----------------------------------------------------------------------------
//// inverse
//inline tens4ds tens4ds::inverse() const
//{
//	matrix c(6,6);
//	
//	// populate c
//	c(0,0) = d[0];
//	c(1,1) = d[2];
//	c(2,2) = d[5];
//	c(3,3) = d[9];
//	c(4,4) = d[14];
//	c(5,5) = d[20];
//	c(0,1) = c(1,0) = d[1];
//	c(0,2) = c(2,0) = d[3];
//	c(0,3) = c(3,0) = d[6];
//	c(0,4) = c(4,0) = d[10];
//	c(0,5) = c(5,0) = d[15];
//	c(1,2) = c(2,1) = d[4];
//	c(1,3) = c(3,1) = d[7];
//	c(1,4) = c(4,1) = d[11];
//	c(1,5) = c(5,1) = d[16];
//	c(2,3) = c(3,2) = d[8];
//	c(2,4) = c(4,2) = d[12];
//	c(2,5) = c(5,2) = d[17];
//	c(3,4) = c(4,3) = d[13];
//	c(3,5) = c(5,3) = d[18];
//	c(4,5) = c(5,4) = d[19];
//	
//	// invert c
//	matrix s = c.inverse();
//	
//	// return inverse
//	tens4ds S;
//	S.d[ 0] = s(0,0);
//	S.d[ 2] = s(1,1);
//	S.d[ 5] = s(2,2);
//	S.d[ 9] = s(3,3)/4.;
//	S.d[14] = s(4,4)/4.;
//	S.d[20] = s(5,5)/4.;
//	S.d[ 1] = s(0,1);
//	S.d[ 3] = s(0,2);
//	S.d[ 6] = s(0,3)/2.;
//	S.d[10] = s(0,4)/2.;
//	S.d[15] = s(0,5)/2.;
//	S.d[ 4] = s(1,2);
//	S.d[ 7] = s(1,3)/2.;
//	S.d[11] = s(1,4)/2.;
//	S.d[16] = s(1,5)/2.;
//	S.d[ 8] = s(2,3)/2.;
//	S.d[12] = s(2,4)/2.;
//	S.d[17] = s(2,5)/2.;
//	S.d[13] = s(3,4)/4.;
//	S.d[18] = s(3,5)/4.;
//	S.d[19] = s(4,5)/4.;
//	
//	return S;
//}
