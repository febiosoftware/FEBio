// NOTE: This file is automatically included from tens5d.h
// Users should not include this file manually!

#include "matrix.h"

// operator +
inline tens6ds tens6ds::operator + (const tens6ds& t) const
{
	tens6ds s;
	for (int i=0; i<NNZ; i++)
		s.d[i] = d[i] + t.d[i];
	
	return s;
}

// operator -
inline tens6ds tens6ds::operator - (const tens6ds& t) const
{
	tens6ds s;
	for (int i=0; i<NNZ; i++)
		s.d[i] = d[i] - t.d[i];

	return s;
}

// operator *
inline tens6ds tens6ds::operator * (double g) const
{
	tens6ds s;
	for (int i=0; i<NNZ; i++)
		s.d[i] = g*d[i];
	
	return s;
}

// operator /
inline tens6ds tens6ds::operator / (double g) const
{
	tens6ds s;
	for (int i=0; i<NNZ; i++)
		s.d[i] = d[i]/g;
	
	return s;
}

// assignment operator +=
inline tens6ds& tens6ds::operator += (const tens6ds& t)
{
	for (int i=0; i<NNZ; i++)
		d[i] += t.d[i];
	
	return (*this);
}

// assignment operator -=
inline tens6ds& tens6ds::operator -= (const tens6ds& t)
{
	for (int i=0; i<NNZ; i++)
		d[i] -= t.d[i];
	
	return (*this);
}

// assignment operator *=
inline tens6ds& tens6ds::operator *= (double g)
{
	for (int i=0; i<NNZ; i++)
		d[i] *= g;
	
	return (*this);
}

// assignment operator /=
inline tens6ds& tens6ds::operator /= (double g)
{
	for (int i=0; i<NNZ; i++)
		d[i] /= g;
	
	return (*this);
}

// unary operator -
inline tens6ds tens6ds::operator - () const
{
	tens6ds s;
	for (int i = 0; i < NNZ; i++)
		s.d[i] = -d[i];

	return s;
}

// intialize to zero
inline void tens6ds::zero()
{
	for (int i = 0; i < NNZ; i++)
		d[i] = 0;
}

//inline void tens6ds::unit()
//{
//	for (int i = 0; i < NNZ; i++)
//		if ((i==0) || (i==171) || (i==323))
//			d[i] = 1.;
//		else
//			d[i] = 0.;
//}

//// contract on the right by a 3o tensor with RC symmetry  Hijk = EijklmnHlmn
//inline tens3drs tens6ds::contract3s(tens3drs H)
//{
//    tens3drs K;
//
//	for (int i = 0; i < 18; i++)
//		K.d[i] = d[0+18*i]*H.d[0] + 2*d[1+18*i]*H.d[1] + 2*d[2+18*i]*H.d[2] + d[3+18*i]*H.d[3] + 2*d[4+18*i]*H.d[4] + d[5+18*i]*H.d[5] + d[6+18*i]*H.d[6] + 2*d[7+18*i]*H.d[7] + 2*d[8+18*i]*H.d[8] + d[9+18*i]*H.d[9] + 2*d[10+18*i]*H.d[10] + d[11+18*i]*H.d[11] + d[12+18*i]*H.d[12] + 2*d[13+18*i]*H.d[13] + 2*d[14+18*i]*H.d[14] + d[15+18*i]*H.d[15] + 2*d[16+18*i]*H.d[16] + d[17+18*i]*H.d[17];
//	
//	return K;
//}