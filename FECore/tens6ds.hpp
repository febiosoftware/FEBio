// NOTE: This file is automatically included from tens6d.h
// Users should not include this file manually!

#include "matrix.h"

// access operator
inline double tens6ds::operator() (int i, int j, int k, int l, int m, int n)
{
	// lookup table convering triplets to a row/column index
	const int LUT[3][3][3] = {
		{{0,1,2},{1,3,4},{2,4,5}},
		{{1,3,4},{3,6,7},{4,7,8}},
		{{2,4,5},{4,7,8},{5,8,9}}};

	// index to start of columns
	const int M[10] = {0,1,3,6,10,15,21,28,37,46};

	int I = LUT[i][j][k];
	int J = LUT[l][m][n];
	return (I <= J ? d[M[J]+I] : d[M[I]+J]);
}

// operator +
inline tens6ds tens6ds::operator + (const tens6ds& t) const
{
	tens6ds s;
	for (int i=0; i<NNZ; i++) s.d[i] = d[i] + t.d[i];
	return s;
}

// operator -
inline tens6ds tens6ds::operator - (const tens6ds& t) const
{
	tens6ds s;
	for (int i=0; i<NNZ; i++) s.d[i] = d[i] - t.d[i];
	return s;
}

// operator *
inline tens6ds tens6ds::operator * (double g) const
{
	tens6ds s;
	for (int i=0; i<NNZ; i++) s.d[i] = g*d[i];
	return s;
}

// operator /
inline tens6ds tens6ds::operator / (double g) const
{
	tens6ds s;
	for (int i=0; i<NNZ; i++) s.d[i] = d[i]/g;
	return s;
}

// assignment operator +=
inline tens6ds& tens6ds::operator += (const tens6ds& t)
{
	for (int i=0; i<NNZ; i++) d[i] += t.d[i];
	return (*this);
}

// assignment operator -=
inline tens6ds& tens6ds::operator -= (const tens6ds& t)
{
	for (int i=0; i<NNZ; i++) d[i] -= t.d[i];
	return (*this);
}

// assignment operator *=
inline tens6ds& tens6ds::operator *= (double g)
{
	for (int i=0; i<NNZ; i++) d[i] *= g;
	return (*this);
}

// assignment operator /=
inline tens6ds& tens6ds::operator /= (double g)
{
	for (int i=0; i<NNZ; i++) d[i] /= g;
	return (*this);
}

// unary operator -
inline tens6ds tens6ds::operator - () const
{
	tens6ds s;
	for (int i = 0; i < NNZ; i++) s.d[i] = -d[i];
	return s;
}

// intialize to zero
inline void tens6ds::zero()
{
	for (int i = 0; i < NNZ; i++) d[i] = 0.0;
}
