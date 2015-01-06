// NOTE: This file is automatically included from tens5d.h
// Users should not include this file manually!

#include "matrix.h"

// operator +
inline tens5d tens5d::operator + (const tens5d& t) const
{
	tens5d s;
	for (int i=0; i<NNZ; i++)
		s.d[i] = d[i] + t.d[i];
	
	return s;
}

// operator -
inline tens5d tens5d::operator - (const tens5d& t) const
{
	tens5d s;
	for (int i=0; i<NNZ; i++)
		s.d[i] = d[i] - t.d[i];

	return s;
}

// operator *
inline tens5d tens5d::operator * (double g) const
{
	tens5d s;
	for (int i=0; i<NNZ; i++)
		s.d[i] = g*d[i];
	
	return s;
}

// operator /
inline tens5d tens5d::operator / (double g) const
{
	tens5d s;
	for (int i=0; i<NNZ; i++)
		s.d[i] = d[i]/g;
	
	return s;
}

// assignment operator +=
inline tens5d& tens5d::operator += (const tens5d& t)
{
	for (int i=0; i<NNZ; i++)
		d[i] += t.d[i];
	
	return (*this);
}

// assignment operator -=
inline tens5d& tens5d::operator -= (const tens5d& t)
{
	for (int i=0; i<NNZ; i++)
		d[i] -= t.d[i];
	
	return (*this);
}

// assignment operator *=
inline tens5d& tens5d::operator *= (double g)
{
	for (int i=0; i<NNZ; i++)
		d[i] *= g;
	
	return (*this);
}

// assignment operator /=
inline tens5d& tens5d::operator /= (double g)
{
	for (int i=0; i<NNZ; i++)
		d[i] /= g;
	
	return (*this);
}

// unary operator -
inline tens5d tens5d::operator - () const
{
	tens5d s;
	for (int i = 0; i < NNZ; i++)
		s.d[i] = -d[i];

	return s;
}

// intialize to zero
inline void tens5d::zero()
{
	for (int i = 0; i < NNZ; i++)
		d[i] = 0;
}

inline void tens5d::unit()
{
	for (int i = 0; i < NNZ; i++)
		if ((i== 0) || (i==63) || (i==107))
			d[i] = 1.;
		else
			d[i] = 0.;
}

// contract on the right by a 3o tensor with RC symmetry  Sij = DijklmHklm
inline mat3ds tens5d::contract3rs(tens3drs H)
{
    mat3ds s;

	s(0,0) = d[0]*H.d[0] + 2*d[1]*H.d[1] + 2*d[2]*H.d[2] + d[3]*H.d[3] + 2*d[4]*H.d[4] + d[5]*H.d[5] + d[6]*H.d[6] + 2*d[7]*H.d[7] + 2*d[8]*H.d[8] + d[9]*H.d[9] + 2*d[10]*H.d[10] + d[11]*H.d[11] + d[12]*H.d[12] + 2*d[13]*H.d[13] + 2*d[14]*H.d[14] + d[15]*H.d[15] + 2*d[16]*H.d[16] + d[17]*H.d[17];
	s(0,1) = d[18]*H.d[0] + 2*d[19]*H.d[1] + 2*d[20]*H.d[2] + d[21]*H.d[3] + 2*d[22]*H.d[4] + d[23]*H.d[5] + d[24]*H.d[6] + 2*d[25]*H.d[7] + 2*d[26]*H.d[8] + d[27]*H.d[9] + 2*d[28]*H.d[10] + d[29]*H.d[11] + d[30]*H.d[12] + 2*d[31]*H.d[13] + 2*d[32]*H.d[14] + d[33]*H.d[15] + 2*d[34]*H.d[16] + d[35]*H.d[17];
	s(0,2) = d[36]*H.d[0] + 2*d[37]*H.d[1] + 2*d[38]*H.d[2] + d[39]*H.d[3] + 2*d[40]*H.d[4] + d[41]*H.d[5] + d[42]*H.d[6] + 2*d[43]*H.d[7] + 2*d[44]*H.d[8] + d[45]*H.d[9] + 2*d[46]*H.d[10] + d[47]*H.d[11] + d[48]*H.d[12] + 2*d[49]*H.d[13] + 2*d[50]*H.d[14] + d[51]*H.d[15] + 2*d[52]*H.d[16] + d[53]*H.d[17];
	s(1,1) = d[54]*H.d[0] + 2*d[55]*H.d[1] + 2*d[56]*H.d[2] + d[57]*H.d[3] + 2*d[58]*H.d[4] + d[59]*H.d[5] + d[60]*H.d[6] + 2*d[61]*H.d[7] + 2*d[62]*H.d[8] + d[63]*H.d[9] + 2*d[64]*H.d[10] + d[65]*H.d[11] + d[66]*H.d[12] + 2*d[67]*H.d[13] + 2*d[68]*H.d[14] + d[69]*H.d[15] + 2*d[70]*H.d[16] + d[71]*H.d[17];
	s(1,2) = d[72]*H.d[0] + 2*d[73]*H.d[1] + 2*d[74]*H.d[2] + d[75]*H.d[3] + 2*d[76]*H.d[4] + d[77]*H.d[5] + d[78]*H.d[6] + 2*d[79]*H.d[7] + 2*d[80]*H.d[8] + d[81]*H.d[9] + 2*d[82]*H.d[10] + d[83]*H.d[11] + d[84]*H.d[12] + 2*d[85]*H.d[13] + 2*d[86]*H.d[14] + d[87]*H.d[15] + 2*d[88]*H.d[16] + d[89]*H.d[17];
	s(2,2) = d[90]*H.d[0] + 2*d[91]*H.d[1] + 2*d[92]*H.d[2] + d[93]*H.d[3] + 2*d[94]*H.d[4] + d[95]*H.d[5] + d[96]*H.d[6] + 2*d[97]*H.d[7] + 2*d[98]*H.d[8] + d[99]*H.d[9] + 2*d[100]*H.d[10] + d[101]*H.d[11] + d[102]*H.d[12] + 2*d[103]*H.d[13] + 2*d[104]*H.d[14] + d[105]*H.d[15] + 2*d[106]*H.d[16] + d[107]*H.d[17];

	return s;
}

// contract on the left by a 2o symmetric tensor  Hklm = SijDijklm
inline tens3drs tens5d::contract2s(mat3ds s)
{
	tens3drs H;

	H.d[0] = s.xx()*d[0] + 2*s.xy()*d[18] + 2*s.xz()*d[36] + s.yy()*d[54] + 2*s.yz()*d[72] + s.zz()*d[90];
	H.d[1] = s.xx()*d[1] + 2*s.xy()*d[19] + 2*s.xz()*d[37] + s.yy()*d[55] + 2*s.yz()*d[73] + s.zz()*d[91];
	H.d[2] = s.xx()*d[2] + 2*s.xy()*d[20] + 2*s.xz()*d[38] + s.yy()*d[56] + 2*s.yz()*d[74] + s.zz()*d[92];
	H.d[3] = s.xx()*d[3] + 2*s.xy()*d[21] + 2*s.xz()*d[39] + s.yy()*d[57] + 2*s.yz()*d[75] + s.zz()*d[93];
	H.d[4] = s.xx()*d[4] + 2*s.xy()*d[22] + 2*s.xz()*d[40] + s.yy()*d[58] + 2*s.yz()*d[76] + s.zz()*d[94];
	H.d[5] = s.xx()*d[5] + 2*s.xy()*d[23] + 2*s.xz()*d[41] + s.yy()*d[59] + 2*s.yz()*d[77] + s.zz()*d[95];
	H.d[6] = s.xx()*d[6] + 2*s.xy()*d[24] + 2*s.xz()*d[42] + s.yy()*d[60] + 2*s.yz()*d[78] + s.zz()*d[96];
	H.d[7] = s.xx()*d[7] + 2*s.xy()*d[25] + 2*s.xz()*d[43] + s.yy()*d[61] + 2*s.yz()*d[79] + s.zz()*d[97];
	H.d[8] = s.xx()*d[8] + 2*s.xy()*d[26] + 2*s.xz()*d[44] + s.yy()*d[62] + 2*s.yz()*d[80] + s.zz()*d[98];
	H.d[9] = s.xx()*d[9] + 2*s.xy()*d[27] + 2*s.xz()*d[45] + s.yy()*d[63] + 2*s.yz()*d[81] + s.zz()*d[99];
	H.d[10] = s.xx()*d[10] + 2*s.xy()*d[28] + 2*s.xz()*d[46] + s.yy()*d[64] + 2*s.yz()*d[82] + s.zz()*d[100];
	H.d[11] = s.xx()*d[11] + 2*s.xy()*d[29] + 2*s.xz()*d[47] + s.yy()*d[65] + 2*s.yz()*d[83] + s.zz()*d[101];
	H.d[12] = s.xx()*d[12] + 2*s.xy()*d[30] + 2*s.xz()*d[48] + s.yy()*d[66] + 2*s.yz()*d[84] + s.zz()*d[102];
	H.d[13] = s.xx()*d[13] + 2*s.xy()*d[31] + 2*s.xz()*d[49] + s.yy()*d[67] + 2*s.yz()*d[85] + s.zz()*d[103];
	H.d[14] = s.xx()*d[14] + 2*s.xy()*d[32] + 2*s.xz()*d[50] + s.yy()*d[68] + 2*s.yz()*d[86] + s.zz()*d[104];
	H.d[15] = s.xx()*d[15] + 2*s.xy()*d[33] + 2*s.xz()*d[51] + s.yy()*d[69] + 2*s.yz()*d[87] + s.zz()*d[105];
	H.d[16] = s.xx()*d[16] + 2*s.xy()*d[34] + 2*s.xz()*d[52] + s.yy()*d[70] + 2*s.yz()*d[88] + s.zz()*d[106];
	H.d[17] = s.xx()*d[17] + 2*s.xy()*d[35] + 2*s.xz()*d[53] + s.yy()*d[71] + 2*s.yz()*d[89] + s.zz()*d[107];

	return H;
}
