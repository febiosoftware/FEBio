// NOTE: This file is automatically included from tens3drs.h
// Users should not include this file manually!

#include "matrix.h"

// operator +
inline tens3ds tens3ds::operator + (const tens3ds& t) const
{
	tens3ds s;
	for (int i=0; i<NNZ; i++)
		s.d[i] = d[i] + t.d[i];
	
	return s;
}

// operator -
inline tens3ds tens3ds::operator - (const tens3ds& t) const
{
	tens3ds s;
	for (int i=0; i<NNZ; i++)
		s.d[i] = d[i] - t.d[i];

	return s;
}

// operator *
inline tens3ds tens3ds::operator * (double g) const
{
	tens3ds s;
	for (int i=0; i<NNZ; i++)
		s.d[i] = g*d[i];
	
	return s;
}

// operator /
inline tens3ds tens3ds::operator / (double g) const
{
	tens3ds s;
	for (int i=0; i<NNZ; i++)
		s.d[i] = d[i]/g;
	
	return s;
}

// assignment operator +=
inline tens3ds& tens3ds::operator += (const tens3ds& t)
{
	for (int i=0; i<NNZ; i++)
		d[i] += t.d[i];
	
	return (*this);
}

// assignment operator -=
inline tens3ds& tens3ds::operator -= (const tens3ds& t)
{
	for (int i=0; i<NNZ; i++)
		d[i] -= t.d[i];
	
	return (*this);
}

// assignment operator *=
inline tens3ds& tens3ds::operator *= (double g)
{
	for (int i=0; i<NNZ; i++)
		d[i] *= g;
	
	return (*this);
}

// assignment operator /=
inline tens3ds& tens3ds::operator /= (double g)
{
	for (int i=0; i<NNZ; i++)
		d[i] /= g;
	
	return (*this);
}

// unary operator -
inline tens3ds tens3ds::operator - () const
{
	tens3ds s;
	for (int i = 0; i < NNZ; i++)
		s.d[i] = -d[i];

	return s;
}

// trace
// Note:  Not quite sure what the right trace for a 3rd order tensor is, this probably isn't it
inline double tens3ds::tr() const
{
	return (d[0] + d[9] + d[17]);
}

// intialize to zero
inline void tens3ds::zero()
{
	for (int i = 0; i < NNZ; i++)
		d[i] = 0;
}

inline void tens3ds::unit()
{
	for (int i = 0; i < NNZ; i++)
		if ((i == 0) || (i == 9) || (i == 17))
			d[i] = 1.;
		else
			d[i] = 0;
}

// contract the right two legs by the dyad formed by a vector  xj = XiTijkXk
inline vec3d tens3ds::contractdyad1(vec3d v)
{
    vec3d x;
	x.x = d[0]*v.x*v.x + 2*d[1]*v.x*v.y + 2*d[2]*v.x*v.z + d[3]*v.y*v.y + 2*d[4]*v.y*v.z + d[5]*v.z*v.z;
	x.y = d[1]*v.x*v.x + 2*d[3]*v.x*v.y + 2*d[4]*v.x*v.z + d[6]*v.y*v.y + 2*d[7]*v.y*v.z + d[8]*v.z*v.z;
	x.z = d[2]*v.x*v.x + 2*d[4]*v.x*v.y + 2*d[5]*v.x*v.z + d[7]*v.y*v.y + 2*d[8]*v.y*v.z + d[9]*v.z*v.z;

	return x;
}

// triple contraction by a similar 3o tensor m = TijkHijk
inline double tens3ds::tripledot3s(tens3ds H)
{
	double m;
	m = d[0]*H.d[0] + 3*d[1]*H.d[1] + 3*d[2]*H.d[2] + 3*d[3]*H.d[3] + 6*d[4]*H.d[4] + 3*d[5]*H.d[5]  + d[6]*H.d[6]  + 3*d[7]*H.d[7]  + 3*d[8]*H.d[8]  + d[9]*H.d[9];

	return m;
}


