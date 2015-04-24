// NOTE: This file is automatically included from tens3d.h
// Users should not include this file manually!

#include "matrix.h"

// operator +
inline tens3dls tens3dls::operator + (const tens3dls& t) const
{
	tens3dls s;
	for (int i=0; i<NNZ; i++)
		s.d[i] = d[i] + t.d[i];
	
	return s;
}

// operator -
inline tens3dls tens3dls::operator - (const tens3dls& t) const
{
	tens3dls s;
	for (int i=0; i<NNZ; i++)
		s.d[i] = d[i] - t.d[i];

	return s;
}

// operator *
inline tens3dls tens3dls::operator * (double g) const
{
	tens3dls s;
	for (int i=0; i<NNZ; i++)
		s.d[i] = g*d[i];
	
	return s;
}

// operator /
inline tens3dls tens3dls::operator / (double g) const
{
	tens3dls s;
	for (int i=0; i<NNZ; i++)
		s.d[i] = d[i]/g;
	
	return s;
}

// assignment operator +=
inline tens3dls& tens3dls::operator += (const tens3dls& t)
{
	for (int i=0; i<NNZ; i++)
		d[i] += t.d[i];
	
	return (*this);
}

// assignment operator -=
inline tens3dls& tens3dls::operator -= (const tens3dls& t)
{
	for (int i=0; i<NNZ; i++)
		d[i] -= t.d[i];
	
	return (*this);
}

// assignment operator *=
inline tens3dls& tens3dls::operator *= (double g)
{
	for (int i=0; i<NNZ; i++)
		d[i] *= g;
	
	return (*this);
}

// assignment operator /=
inline tens3dls& tens3dls::operator /= (double g)
{
	for (int i=0; i<NNZ; i++)
		d[i] /= g;
	
	return (*this);
}

// unary operator -
inline tens3dls tens3dls::operator - () const
{
	tens3dls s;
	for (int i = 0; i < NNZ; i++)
		s.d[i] = -d[i];

	return s;
}

// trace
// Note:  Not quite sure what the right trace for a 3rd order tensor is, this probably isn't it
inline double tens3dls::tr() const
{
	return (d[0] + d[10] + d[17]);
}

// intialize to zero
inline void tens3dls::zero()
{
	for (int i = 0; i < NNZ; i++)
		d[i] = 0;
}

// unit tensor???
inline void tens3dls::unit()
{
	for (int i = 0; i < NNZ; i++)
		if ((i == 0) || (i == 9) || (i == 17))
			d[i] = 1.;
		else
			d[i] = 0;
}

// convert a LS tensor to an unsymmetric tensor class
inline tens3d tens3dls::LStoUnsym()
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

// calculate the transpose ((G_KJi)T = G_iJK)
inline tens3drs tens3dls::transpose()
{
	tens3drs GRC;

	GRC.d[0] =  d[0];
	GRC.d[1] =  d[3];
	GRC.d[2] =  d[6];
	GRC.d[3] =  d[9];
	GRC.d[4] =  d[12];
	GRC.d[5] =  d[15];
	GRC.d[6] =  d[1];
	GRC.d[7] =  d[4];
	GRC.d[8] =  d[7];
	GRC.d[9] =  d[10];
	GRC.d[10] = d[13];
	GRC.d[11] = d[16];
	GRC.d[12] = d[2];
	GRC.d[13] = d[5];
	GRC.d[14] = d[8];
	GRC.d[15] = d[11];
	GRC.d[16] = d[14];
	GRC.d[17] = d[17];

	return GRC;
}

// multiply by a 2o tensor on the right (G_KJi * F_iI)
inline tens3dls tens3dls::multiply2right(const mat3d& F)
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

	G.d[9] = d[9]*F(0,0) + d[10]*F(1,0) + d[11]*F(2,0);
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