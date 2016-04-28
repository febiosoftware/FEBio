// NOTE: This file is automatically included from tens5d.h
// Users should not include this file manually!

#include <assert.h>

// access operator
// TODO: implement this
inline double tens5ds::operator () (int i, int j, int k, int l, int m) const
{
	assert(false);
	return 0.0;
}

// operator +
inline tens5ds tens5ds::operator + (const tens5ds& t) const
{
	tens5ds s;
	for (int i=0; i<NNZ; i++) s.d[i] = d[i] + t.d[i];
	return s;
}

// operator -
inline tens5ds tens5ds::operator - (const tens5ds& t) const
{
	tens5ds s;
	for (int i=0; i<NNZ; i++) s.d[i] = d[i] - t.d[i];
	return s;
}

// operator *
inline tens5ds tens5ds::operator * (double g) const
{
	tens5ds s;
	for (int i=0; i<NNZ; i++) s.d[i] = g*d[i];

	return s;
}

// operator /
inline tens5ds tens5ds::operator / (double g) const
{
	tens5ds s;
	for (int i=0; i<NNZ; i++) s.d[i] = d[i]/g;
	return s;
}

// assignment operator +=
inline tens5ds& tens5ds::operator += (const tens5ds& t)
{
	for (int i=0; i<NNZ; i++) d[i] += t.d[i];
	return (*this);
}

// assignment operator -=
inline tens5ds& tens5ds::operator -= (const tens5ds& t)
{
	for (int i=0; i<NNZ; i++) d[i] -= t.d[i];
	return (*this);
}

// assignment operator *=
inline tens5ds& tens5ds::operator *= (double g)
{
	for (int i=0; i<NNZ; i++) d[i] *= g;
	return (*this);
}

// assignment operator /=
inline tens5ds& tens5ds::operator /= (double g)
{
	for (int i=0; i<NNZ; i++) d[i] /= g;
	return (*this);
}

// unary operator -
inline tens5ds tens5ds::operator - () const
{
	tens5ds s;
	for (int i = 0; i < NNZ; i++) s.d[i] = -d[i];
	return s;
}

// intialize to zero
inline void tens5ds::zero()
{
	for (int i = 0; i < NNZ; i++) d[i] = 0.0;
}
