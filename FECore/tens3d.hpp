// NOTE: This file is automatically included from tens3drs.h
// Users should not include this file manually!

#include "matrix.h"

// operator +
inline tens3d tens3d::operator + (const tens3d& t) const
{
	tens3d s;
	for (int i=0; i<NNZ; i++)
		s.d[i] = d[i] + t.d[i];
	
	return s;
}

// operator -
inline tens3d tens3d::operator - (const tens3d& t) const
{
	tens3d s;
	for (int i=0; i<NNZ; i++)
		s.d[i] = d[i] - t.d[i];

	return s;
}

// operator *
inline tens3d tens3d::operator * (double g) const
{
	tens3d s;
	for (int i=0; i<NNZ; i++)
		s.d[i] = g*d[i];
	
	return s;
}

// operator /
inline tens3d tens3d::operator / (double g) const
{
	tens3d s;
	for (int i=0; i<NNZ; i++)
		s.d[i] = d[i]/g;
	
	return s;
}

// assignment operator +=
inline tens3d& tens3d::operator += (const tens3d& t)
{
	for (int i=0; i<NNZ; i++)
		d[i] += t.d[i];
	
	return (*this);
}

// assignment operator -=
inline tens3d& tens3d::operator -= (const tens3d& t)
{
	for (int i=0; i<NNZ; i++)
		d[i] -= t.d[i];
	
	return (*this);
}

// assignment operator *=
inline tens3d& tens3d::operator *= (double g)
{
	for (int i=0; i<NNZ; i++)
		d[i] *= g;
	
	return (*this);
}

// assignment operator /=
inline tens3d& tens3d::operator /= (double g)
{
	for (int i=0; i<NNZ; i++)
		d[i] /= g;
	
	return (*this);
}

// unary operator -
inline tens3d tens3d::operator - () const
{
	tens3d s;
	for (int i = 0; i < NNZ; i++)
		s.d[i] = -d[i];

	return s;
}

// trace
// Note:  Not quite sure what the right trace for a 3rd order tensor is, this probably isn't it
inline double tens3d::tr() const
{
	return (d[0] + d[13] + d[26]);
}

// intialize to zero
inline void tens3d::zero()
{
	for (int i = 0; i < NNZ; i++)
		d[i] = 0;
}

inline void tens3d::unit()
{
	for (int i = 0; i < NNZ; i++)
		if ((i == 0) || (i == 9) || (i == 17))
			d[i] = 1.;
		else
			d[i] = 0;
}

inline tens3ds tens3d::symm()
{
	tens3ds t;

	t.d[0] +=  d[0]; 
	t.d[1] += (d[1] + d[3]  + d[9])/3.; 
	t.d[2] += (d[2] + d[6]  + d[18])/3.;
	t.d[3] += (d[4] + d[10] + d[12])/3.; 
	t.d[4] += (d[5] + d[11] + d[21] + d[7] + d[19] + d[15])/6.; 
	t.d[5] += (d[8] + d[20] + d[24])/3.;
	t.d[6] +=  d[13]; 
	t.d[7] += (d[14] + d[16] + d[22])/3.;
	t.d[8] += (d[17] + d[23] + d[25])/3.;
	t.d[9] +=  d[26]; 

	return t;
}
