// NOTE: This file is automatically included from tens3drs.h
// Users should not include this file manually!

#include "matrix.h"

// operator +
inline tens3drs tens3drs::operator + (const tens3drs& t) const
{
	tens3drs s;
	for (int i=0; i<NNZ; i++)
		s.d[i] = d[i] + t.d[i];
	
	return s;
}

// operator -
inline tens3drs tens3drs::operator - (const tens3drs& t) const
{
	tens3drs s;
	for (int i=0; i<NNZ; i++)
		s.d[i] = d[i] - t.d[i];

	return s;
}

// operator *
inline tens3drs tens3drs::operator * (double g) const
{
	tens3drs s;
	for (int i=0; i<NNZ; i++)
		s.d[i] = g*d[i];
	
	return s;
}

// operator /
inline tens3drs tens3drs::operator / (double g) const
{
	tens3drs s;
	for (int i=0; i<NNZ; i++)
		s.d[i] = d[i]/g;
	
	return s;
}

// assignment operator +=
inline tens3drs& tens3drs::operator += (const tens3drs& t)
{
	for (int i=0; i<NNZ; i++)
		d[i] += t.d[i];
	
	return (*this);
}

// assignment operator -=
inline tens3drs& tens3drs::operator -= (const tens3drs& t)
{
	for (int i=0; i<NNZ; i++)
		d[i] -= t.d[i];
	
	return (*this);
}

// assignment operator *=
inline tens3drs& tens3drs::operator *= (double g)
{
	for (int i=0; i<NNZ; i++)
		d[i] *= g;
	
	return (*this);
}

// assignment operator /=
inline tens3drs& tens3drs::operator /= (double g)
{
	for (int i=0; i<NNZ; i++)
		d[i] /= g;
	
	return (*this);
}

// unary operator -
inline tens3drs tens3drs::operator - () const
{
	tens3drs s;
	for (int i = 0; i < NNZ; i++)
		s.d[i] = -d[i];

	return s;
}

// trace
// Note:  Not quite sure what the right trace for a 3rd order tensor is, this probably isn't it
inline double tens3drs::tr() const
{
	return (d[0] + d[9] + d[17]);
}

// intialize to zero
inline void tens3drs::zero()
{
	for (int i = 0; i < NNZ; i++)
		d[i] = 0;
}

inline void tens3drs::unit()
{
	for (int i = 0; i < NNZ; i++)
		if ((i == 0) || (i == 9) || (i == 17))
			d[i] = 1.;
		else
			d[i] = 0;
}

// contract the right two legs by the dyad formed by a vector  xi = GijkXjXk
inline vec3d tens3drs::contractdyad1(vec3d v)
{
    vec3d x;
	x.x = d[0]*v.x*v.x + 2*d[1]*v.x*v.y + 2*d[2]*v.x*v.z + d[3]*v.y*v.y + 2*d[4]*v.y*v.z + d[5]*v.z*v.z;
	x.y = d[6]*v.x*v.x + 2*d[7]*v.x*v.y + 2*d[8]*v.x*v.z + d[9]*v.y*v.y + 2*d[10]*v.y*v.z + d[11]*v.z*v.z;
	x.z = d[12]*v.x*v.x + 2*d[13]*v.x*v.y + 2*d[14]*v.x*v.z + d[15]*v.y*v.y + 2*d[16]*v.y*v.z + d[17]*v.z*v.z;

	return x;
}

// contract the right two legs by a symmetric 2o tensor  xi = GijkSjk
inline vec3d tens3drs::contract2s(mat3ds s)
{
    vec3d x;
	x.x = d[0]*s.xx() + 2*d[1]*s.xy() + 2*d[2]*s.xz() + d[3]*s.yy() + 2*d[4]*s.yz() + d[5]*s.zz();
	x.y = d[6]*s.xx() + 2*d[7]*s.xy() + 2*d[8]*s.xz() + d[9]*s.yy() + 2*d[10]*s.yz() + d[11]*s.zz();
	x.z = d[12]*s.xx() + 2*d[13]*s.xy() + 2*d[14]*s.xz() + d[15]*s.yy() + 2*d[16]*s.yz() + d[17]*s.zz();

	return x;
}

// triple contraction by a similar 3o tensor m = GijkHijk
inline double tens3drs::tripledot3rs(tens3drs H)
{
	double m;
	m = d[0]*H.d[0] + 2*d[1]*H.d[1] + 2*d[2]*H.d[2] + d[3]*H.d[3] + 2*d[4]*H.d[4] + d[5]*H.d[5]  + d[6]*H.d[6]  + 2*d[7]*H.d[7]  + 2*d[8]*H.d[8]  + d[9]*H.d[9]  + 2*d[10]*H.d[10]  + d[11]*H.d[11]  + d[12]*H.d[12]  + 2*d[13]*H.d[13]  + 2*d[14]*H.d[14]  + d[15]*H.d[15]  + 2*d[16]*H.d[16]  + d[17]*H.d[17];  

	return m;
}

// contract the right two legs by the dyad formed by a vector  xi = GijkXjVk
inline vec3d tens3drs::contractdyad2(vec3d v, vec3d w)
{
    vec3d x;
	x.x = d[0]*v.x*w.x + d[1]*(v.x*w.y + v.y*w.x) + d[2]*(v.x*w.z + v.z*w.x) + d[3]*v.y*w.y + d[4]*(v.y*w.z + v.z*w.y) + d[5]*v.z*w.z;
	x.y = d[6]*v.x*w.x + d[7]*(v.x*w.y + v.y*w.x) + d[8]*(v.x*w.z + v.z*w.x) + d[9]*v.y*w.y + d[10]*(v.y*w.z + v.z*w.y) + d[11]*v.z*w.z;
	x.z = d[12]*v.x*w.x + d[13]*(v.x*w.y + v.y*w.x) + d[14]*(v.x*w.z + v.z*w.x) + d[15]*v.y*w.y + d[16]*(v.y*w.z + v.z*w.y) + d[17]*v.z*w.z;

	return x;
}

// convert a RS tensor to an unsymmetric tensor class
inline tens3d tens3drs::RStoUnsym()
{
	tens3d G;

	G.d[0] = d[0];
	G.d[1] = d[1];
	G.d[2] = d[2];
	G.d[3] = d[1];
	G.d[4] = d[3];
	G.d[5] = d[4];
	G.d[6] = d[2];
	G.d[7] = d[4];
	G.d[8] = d[5];
	G.d[9] = d[6];
	G.d[10] = d[7];
	G.d[11] = d[8];
	G.d[12] = d[7];
	G.d[13] = d[9];
	G.d[14] = d[10];
	G.d[15] = d[8];	
	G.d[16] = d[10];
	G.d[17] = d[11];
	G.d[18] = d[12];
	G.d[19] = d[13];
	G.d[20] = d[14];
	G.d[21] = d[13];
	G.d[22] = d[15];
	G.d[23] = d[16];
	G.d[24] = d[14];
	G.d[25] = d[16];
	G.d[26] = d[17];

	return G;
}

// calculate the transpose ((G_iJK)T = G_KJi)
inline tens3dls tens3drs::transpose()
{
	tens3dls GLC;

	GLC.d[0] =  d[0];
	GLC.d[3] =  d[1];
	GLC.d[6] =  d[2];
	GLC.d[9] =  d[3];
	GLC.d[12] =  d[4];
	GLC.d[15] =  d[5];
	GLC.d[1] =  d[6];
	GLC.d[4] =  d[7];
	GLC.d[7] =  d[8];
	GLC.d[10] =  d[9];
	GLC.d[13] = d[10];
	GLC.d[16] = d[11];
	GLC.d[2] = d[12];
	GLC.d[5] = d[13];
	GLC.d[8] = d[14];
	GLC.d[11] = d[15];
	GLC.d[14] = d[16];
	GLC.d[17] = d[17];

	return GLC;
}

// contract each leg by a 2o tensor (intended to calculate the inverse deformation hessian according to Finv_Ii * G_iJK * Finv_Jj * Fin_Kk)
inline void tens3drs::contractleg2(mat3d F, int leg)
{
	tens3drs G = *this;
	
	if (leg == 1)
	{
		d[0] = F[0][0]*G.d[0] + F[0][1]*G.d[6] + F[0][2]*G.d[12];
		d[1] = F[0][0]*G.d[1] + F[0][1]*G.d[7] + F[0][2]*G.d[13];
		d[2] = F[0][0]*G.d[2] + F[0][1]*G.d[8] + F[0][2]*G.d[14];
		d[3] = F[0][0]*G.d[3] + F[0][1]*G.d[9] + F[0][2]*G.d[15];
		d[4] = F[0][0]*G.d[4] + F[0][1]*G.d[10] + F[0][2]*G.d[16];
		d[5] = F[0][0]*G.d[5] + F[0][1]*G.d[11] + F[0][2]*G.d[17];
		d[6] = F[1][0]*G.d[0] + F[1][1]*G.d[6] + F[1][2]*G.d[12];
		d[7] = F[1][0]*G.d[1] + F[1][1]*G.d[7] + F[1][2]*G.d[13];
		d[8] = F[1][0]*G.d[2] + F[1][1]*G.d[8] + F[1][2]*G.d[14];
		d[9] = F[1][0]*G.d[3] + F[1][1]*G.d[9] + F[1][2]*G.d[15];
		d[10] = F[1][0]*G.d[4] + F[1][1]*G.d[10] + F[1][2]*G.d[16];
		d[11] = F[1][0]*G.d[5] + F[1][1]*G.d[11] + F[1][2]*G.d[17];
		d[12] = F[2][0]*G.d[0] + F[2][1]*G.d[6] + F[2][2]*G.d[12];
		d[13] = F[2][0]*G.d[1] + F[2][1]*G.d[7] + F[2][2]*G.d[13];
		d[14] = F[2][0]*G.d[2] + F[2][1]*G.d[8] + F[2][2]*G.d[14];
		d[15] = F[2][0]*G.d[3] + F[2][1]*G.d[9] + F[2][2]*G.d[15];
		d[16] = F[2][0]*G.d[4] + F[2][1]*G.d[10] + F[2][2]*G.d[16];
		d[17] = F[2][0]*G.d[5] + F[2][1]*G.d[11] + F[2][2]*G.d[17];
	}
	else if (leg == 2)
	{
		d[0] = G.d[0]*F[0][0] + G.d[1]*F[1][0] + G.d[2]*F[2][0];
		d[1] = G.d[1]*F[0][0] + G.d[3]*F[1][0] + G.d[4]*F[2][0];
		d[2] = G.d[2]*F[0][0] + G.d[4]*F[1][0] + G.d[5]*F[2][0];
		d[3] = G.d[1]*F[0][1] + G.d[3]*F[1][1] + G.d[4]*F[2][1];
		d[4] = G.d[2]*F[0][1] + G.d[4]*F[1][1] + G.d[5]*F[2][1];
		d[5] = G.d[2]*F[0][2] + G.d[4]*F[1][2] + G.d[5]*F[2][2];
		d[6] = G.d[6]*F[0][0] + G.d[7]*F[1][0] + G.d[8]*F[2][0];
		d[7] = G.d[7]*F[0][0] + G.d[9]*F[1][0] + G.d[10]*F[2][0];
		d[8] = G.d[8]*F[0][0] + G.d[10]*F[1][0] + G.d[11]*F[2][0];
		d[9] = G.d[7]*F[0][1] + G.d[9]*F[1][1] + G.d[10]*F[2][1];
		d[10] = G.d[8]*F[0][1] + G.d[10]*F[1][1] + G.d[11]*F[2][1];
		d[11] = G.d[8]*F[0][2] + G.d[10]*F[1][2] + G.d[11]*F[2][2];
		d[12] = G.d[12]*F[0][0] + G.d[13]*F[1][0] + G.d[4]*F[2][0];
		d[13] = G.d[13]*F[0][0] + G.d[14]*F[1][0] + G.d[16]*F[2][0];
		d[14] = G.d[14]*F[0][0] + G.d[16]*F[1][0] + G.d[17]*F[2][0];
		d[15] = G.d[13]*F[0][1] + G.d[15]*F[1][1] + G.d[16]*F[2][1];
		d[16] = G.d[14]*F[0][1] + G.d[16]*F[1][1] + G.d[17]*F[2][1];
		d[17] = G.d[14]*F[0][2] + G.d[16]*F[1][2] + G.d[17]*F[2][2];
	}
	else if (leg == 3)
	{
		d[0] = G.d[0]*F[0][0] + G.d[1]*F[1][0] + G.d[2]*F[2][0];
		d[1] = G.d[0]*F[0][1] + G.d[1]*F[1][1] + G.d[2]*F[2][1];
		d[2] = G.d[0]*F[0][2] + G.d[1]*F[1][2] + G.d[2]*F[2][2];
		d[3] = G.d[1]*F[0][1] + G.d[3]*F[1][1] + G.d[4]*F[2][1];
		d[4] = G.d[1]*F[0][2] + G.d[3]*F[1][2] + G.d[4]*F[2][2];
		d[5] = G.d[2]*F[0][2] + G.d[4]*F[1][2] + G.d[5]*F[2][2];
		d[6] = G.d[6]*F[0][0] + G.d[7]*F[1][0] + G.d[8]*F[2][0];
		d[7] = G.d[6]*F[0][1] + G.d[7]*F[1][1] + G.d[8]*F[2][1];
		d[8] = G.d[6]*F[0][2] + G.d[7]*F[1][2] + G.d[8]*F[2][2];
		d[9] = G.d[7]*F[0][1] + G.d[9]*F[1][1] + G.d[10]*F[2][1];
		d[10] = G.d[7]*F[0][2] + G.d[9]*F[1][2] + G.d[10]*F[2][2];
		d[11] = G.d[8]*F[0][2] + G.d[10]*F[1][2] + G.d[11]*F[2][2];
		d[12] = G.d[12]*F[0][0] + G.d[13]*F[1][0] + G.d[14]*F[2][0];
		d[13] = G.d[12]*F[0][1] + G.d[13]*F[1][1] + G.d[14]*F[2][1];
		d[14] = G.d[12]*F[0][2] + G.d[13]*F[1][2] + G.d[14]*F[2][2];
		d[15] = G.d[13]*F[0][1] + G.d[15]*F[1][1] + G.d[16]*F[2][1];
		d[16] = G.d[13]*F[0][2] + G.d[15]*F[1][2] + G.d[16]*F[2][2];
		d[17] = G.d[14]*F[0][2] + G.d[16]*F[1][2] + G.d[17]*F[2][2];
	}
}

// multiply by a 2o tensor on the left (F_Ii * G_iJK)
inline tens3drs tens3drs::multiply2left(mat3d F)
{
	tens3drs G;
/*
	G.d[0] = F[0][0]*d.[0] + F[0][1]*d.[6] + F[0][2]*d.[12];
	G.d[1] = F[0][0]*d.[1] + F[0][1]*d.[7] + F[0][2]*d.[13];
	G.d[2] = F[0][0]*d.[2] + F[0][1]*d.[8] + F[0][2]*d.[14];
	G.d[3] = F[0][0]*d.[3] + F[0][1]*d.[9] + F[0][2]*d.[15];
	G.d[4] = F[0][0]*d.[4] + F[0][1]*d.[10] + F[0][2]*d.[16];
	G.d[5] = F[0][0]*d.[5] + F[0][1]*d.[11] + F[0][2]*d.[17];

	G.d[6] = F[1][0]*d.[0] + F[1][1]*d.[6] + F[1][2]*d.[12];
	G.d[7] = F[1][0]*d.[1] + F[1][1]*d.[7] + F[1][2]*d.[13];
	G.d[8] = F[1][0]*d.[2] + F[1][1]*d.[8] + F[1][2]*d.[14];
	G.d[9] = F[1][0]*d.[3] + F[1][1]*d.[9] + F[1][2]*d.[15];
	G.d[10] = F[1][0]*d.[4] + F[1][1]*d.[10] + F[1][2]*d.[16];
	G.d[11] = F[1][0]*d.[5] + F[1][1]*d.[11] + F[1][2]*d.[17];

	G.d[12] = F[2][0]*d.[0] + F[2][1]*d.[6] + F[2][2]*d.[12];
	G.d[13] = F[2][0]*d.[1] + F[2][1]*d.[7] + F[2][2]*d.[13];
	G.d[14] = F[2][0]*d.[2] + F[2][1]*d.[8] + F[2][2]*d.[14];
	G.d[15] = F[2][0]*d.[3] + F[2][1]*d.[9] + F[2][2]*d.[15];
	G.d[16] = F[2][0]*d.[4] + F[2][1]*d.[10] + F[2][2]*d.[16];
	G.d[17] = F[2][0]*d.[5] + F[2][1]*d.[11] + F[2][2]*d.[17];
*/
	return G;
}
