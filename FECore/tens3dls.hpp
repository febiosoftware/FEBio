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

inline void tens3dls::unit()
{
	for (int i = 0; i < NNZ; i++)
		if ((i == 0) || (i == 9) || (i == 17))
			d[i] = 1.;
		else
			d[i] = 0;
}

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

//// contract the right two legs by the dyad formed by a vector  xi = GijkXjXk
//inline vec3d tens3drs::contractdyad1(vec3d v)
//{
//    vec3d x;
//	x.x = d[0]*v.x*v.x + 2*d[1]*v.x*v.y + 2*d[2]*v.x*v.z + d[3]*v.y*v.y + 2*d[4]*v.y*v.z + d[5]*v.z*v.z;
//	x.y = d[6]*v.x*v.x + 2*d[7]*v.x*v.y + 2*d[8]*v.x*v.z + d[9]*v.y*v.y + 2*d[10]*v.y*v.z + d[11]*v.z*v.z;
//	x.z = d[12]*v.x*v.x + 2*d[13]*v.x*v.y + 2*d[14]*v.x*v.z + d[15]*v.y*v.y + 2*d[16]*v.y*v.z + d[17]*v.z*v.z;
//
//	return x;
//}
//
//// contract the right two legs by a symmetric 2o tensor  xi = GijkSjk
//inline vec3d tens3drs::contract2s(mat3ds s)
//{
//    vec3d x;
//	x.x = d[0]*s.xx() + 2*d[1]*s.xy() + 2*d[2]*s.xz() + d[3]*s.yy() + 2*d[4]*s.yz() + d[5]*s.zz();
//	x.y = d[6]*s.xx() + 2*d[7]*s.xy() + 2*d[8]*s.xz() + d[9]*s.yy() + 2*d[10]*s.yz() + d[11]*s.zz();
//	x.z = d[12]*s.xx() + 2*d[13]*s.xy() + 2*d[14]*s.xz() + d[15]*s.yy() + 2*d[16]*s.yz() + d[17]*s.zz();
//
//	return x;
//}
//
//// triple contraction by a similar 3o tensor m = GijkHijk
//inline double tens3drs::tripledot3rs(tens3drs H)
//{
//	double m;
//	m = d[0]*H.d[0] + 2*d[1]*H.d[1] + 2*d[2]*H.d[2] + d[3]*H.d[3] + 2*d[4]*H.d[4] + d[5]*H.d[5]  + d[6]*H.d[6]  + 2*d[7]*H.d[7]  + 2*d[8]*H.d[8]  + d[9]*H.d[9]  + 2*d[10]*H.d[10]  + d[11]*H.d[11]  + d[12]*H.d[12]  + 2*d[13]*H.d[13]  + 2*d[14]*H.d[14]  + d[15]*H.d[15]  + 2*d[16]*H.d[16]  + d[17]*H.d[17];  
//
//	return m;
//}