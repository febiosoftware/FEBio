#include "stdafx.h"
#include "tens4d.h"
#include <math.h>

// operator + 
tens4ds tens4ds::operator + (const tens4ds& t) const
{
	tens4ds s;
	for (int i=0; i<NNZ; i++)
		s.d[i] = d[i] + t.d[i];
	return s;
}

// operator -
tens4ds tens4ds::operator - (const tens4ds& t) const
{
	tens4ds s;
	for (int i=0; i<NNZ; i++)
		s.d[i] = d[i] - t.d[i];
	return s;
}

// operator *
tens4ds tens4ds::operator * (double g) const
{
	tens4ds s;
	for (int i=0; i<NNZ; i++)
		s.d[i] = g*d[i];
	return s;
}

// operator /
tens4ds tens4ds::operator / (double g) const
{
	tens4ds s;
	for (int i=0; i<NNZ; i++)
		s.d[i] = d[i]/g;
	return s;
}

// assignment operator +=
tens4ds& tens4ds::operator += (const tens4ds& t)
{
	for (int i=0; i<NNZ; i++)
		d[i] += t.d[i];
	return (*this);
}

// assignment operator -=
tens4ds& tens4ds::operator -= (const tens4ds& t)
{
	for (int i=0; i<NNZ; i++)
		d[i] -= t.d[i];
	return (*this);
}

// assignment operator *=
tens4ds& tens4ds::operator *= (double g)
{
	for (int i=0; i<NNZ; i++)
		d[i] *= g;
	return (*this);
}

// assignment operator /=
tens4ds& tens4ds::operator /= (double g)
{
	for (int i=0; i<NNZ; i++)
		d[i] /= g;
	return (*this);
}

// extract 6x6 matrix
void tens4ds::extract(double D[6][6])
{
	D[0][0] = d[0];  D[0][1] = d[1];  D[0][2] = d[3];  D[0][3] = d[6];  D[0][4] = d[10]; D[0][5] = d[15];
	D[1][0] = d[1];  D[1][1] = d[2];  D[1][2] = d[4];  D[1][3] = d[7];  D[1][4] = d[11]; D[1][5] = d[16];
	D[2][0] = d[3];  D[2][1] = d[4];  D[2][2] = d[5];  D[2][3] = d[8];  D[2][4] = d[12]; D[2][5] = d[17];
	D[3][0] = d[6];  D[3][1] = d[7];  D[3][2] = d[8];  D[3][3] = d[9];  D[3][4] = d[13]; D[3][5] = d[18];
	D[4][0] = d[10]; D[4][1] = d[11]; D[4][2] = d[12]; D[4][3] = d[13]; D[4][4] = d[14]; D[4][5] = d[19];
	D[5][0] = d[15]; D[5][1] = d[16]; D[5][2] = d[17]; D[5][3] = d[18]; D[5][4] = d[19]; D[5][5] = d[20];
}

//-----------------------------------------------------------------------------
//! This function checks the positive definiteness of a 4th order tensor
//! having both major and minor symmetries. The function does not do an
//! exhaustive test, in the sense it can only detect failure. If a tensor passes
//! it is not guaranteed that the tensor is indeed positive-definite.
bool IsPositiveDefinite(const tens4ds& t)
{
	// test 1. all diagonal entries have to be positive
	if (t(0,0) <= 0) return false;
	if (t(1,1) <= 0) return false;
	if (t(2,2) <= 0) return false;
	if (t(3,3) <= 0) return false;
	if (t(4,4) <= 0) return false;
	if (t(5,5) <= 0) return false;

	// test 2. t(i,i)+t(j,j) > 2t(i,j)
	int i, j;
	for (i=0; i<6; ++i)
	{
		for (j=i+1; j<6; ++j)
		{
			if (t(i,i)+t(j,j) <= 2*t(i,j))
			{
				return false;
			}
		}
	}

	// test 3. the element with largest modulus lies on the main diagonal
	double l = -1, v;
	bool d = false;
	for (i=0; i<6; ++i)
	{
		for (j=i; j<6; ++j)
		{
			v = fabs(t(i,j));
			if (v > l)
			{
				l = v;
				d = (i==j);
			}
		}
	}

	if (d == false) 
	{
		return false;
	}

	// if all tests pass, it is not guaranteed that the tensor is indeed positive-definite
	// but we'd have some good reasons to believe so.
	return true;
}

//-----------------------------------------------------------------------------
// (a dyad1s a)_ijkl = a_ij a_kl
tens4ds dyad1s(const mat3ds& a)
{
	tens4ds c;
	c.d[0] = a.xx()*a.xx();
	c.d[1] = a.xx()*a.yy();
	c.d[3] = a.xx()*a.zz();
	c.d[6] = a.xx()*a.xy();
	c.d[10] = a.xx()*a.xz();
	c.d[15] = a.xx()*a.yz();

	c.d[2] = a.yy()*a.yy();
	c.d[4] = a.yy()*a.zz();
	c.d[7] = a.yy()*a.xy();
	c.d[11] = a.yy()*a.xz();
	c.d[16] = a.yy()*a.yz();

	c.d[5] = a.zz()*a.zz();
	c.d[8] = a.zz()*a.xy();
	c.d[12] = a.zz()*a.xz();
	c.d[17] = a.zz()*a.yz();

	c.d[9] = a.xy()*a.xy();
	c.d[13] = a.xy()*a.xz();
	c.d[18] = a.xy()*a.yz();

	c.d[14] = a.xz()*a.xz();
	c.d[19] = a.xz()*a.yz();

	c.d[20] = a.yz()*a.yz();
	return c;
}

//-----------------------------------------------------------------------------
// (a dyad1s b)_ijkl = a_ij b_kl + b_ij a_kl
tens4ds dyad1s(const mat3ds& a, const mat3ds& b)
{
	tens4ds c;
	c.d[0] = 2*a.xx()*b.xx();
	c.d[1] = a.xx()*b.yy() + b.xx()*a.yy();
	c.d[3] = a.xx()*b.zz() + b.xx()*a.zz();
	c.d[6] = a.xx()*b.xy() + b.xx()*a.xy();
	c.d[10] = a.xx()*b.xz() + b.xx()*a.xz();
	c.d[15] = a.xx()*b.yz() + b.xx()*a.yz();
	
	c.d[2] = 2*a.yy()*b.yy();
	c.d[4] = a.yy()*b.zz() + b.yy()*a.zz();
	c.d[7] = a.yy()*b.xy() + b.yy()*a.xy();
	c.d[11] = a.yy()*b.xz() + b.yy()*a.xz();
	c.d[16] = a.yy()*b.yz() + b.yy()*a.yz();
	
	c.d[5] = 2*a.zz()*b.zz();
	c.d[8] = a.zz()*b.xy() + b.zz()*a.xy();
	c.d[12] = a.zz()*b.xz() + b.zz()*a.xz();
	c.d[17] = a.zz()*b.yz() + b.zz()*a.yz();
	
	c.d[9] = 2*a.xy()*b.xy();
	c.d[13] = a.xy()*b.xz() + b.xy()*a.xz();
	c.d[18] = a.xy()*b.yz() + b.xy()*a.yz();
	
	c.d[14] = 2*a.xz()*b.xz();
	c.d[19] = a.xz()*b.yz() + b.xz()*a.yz();
	
	c.d[20] = 2*a.yz()*b.yz();
	return c;
}

//-----------------------------------------------------------------------------
// (a dyad4s a)_ijkl = (a_ik a_jl + a_il a_jk)/2
tens4ds dyad4s(const mat3ds& a)
{
	tens4ds c;
	c.d[0] =  a.xx()*a.xx();
	c.d[1] =  a.xy()*a.xy();
	c.d[3] =  a.xz()*a.xz();
	c.d[6] =  a.xx()*a.xy();
	c.d[10] = a.xx()*a.xz();
	c.d[15] = a.xy()*a.xz();
	
	c.d[2] =  a.yy()*a.yy();
	c.d[4] =  a.yz()*a.yz();
	c.d[7] =  a.xy()*a.yy();
	c.d[11] = a.xy()*a.yz();
	c.d[16] = a.yy()*a.yz();
	
	c.d[5] =  a.zz()*a.zz();
	c.d[8] =  a.xz()*a.yz();
	c.d[12] = a.xz()*a.zz();
	c.d[17] = a.yz()*a.zz();
	
	c.d[9] =  (a.xx()*a.yy() + a.xy()*a.xy())/2.;
	c.d[13] = (a.xx()*a.yz() + a.xz()*a.xy())/2.;
	c.d[18] = (a.xy()*a.yz() + a.xz()*a.yy())/2.;
	
	c.d[14] = (a.xx()*a.zz() + a.xz()*a.xz())/2.;
	c.d[19] = (a.xy()*a.zz() + a.xz()*a.yz())/2.;
	
	c.d[20] = (a.yy()*a.zz() + a.yz()*a.yz())/2.;
	return c;
	
}

//-----------------------------------------------------------------------------
// (a dyad4s b)_ijkl = (a_ik b_jl + a_il b_jk)/2 +  (b_ik a_jl + b_il a_jk)/2
tens4ds dyad4s(const mat3ds& a, const mat3ds& b)
{
	tens4ds c;
	c.d[0] =  2*a.xx()*b.xx();
	c.d[1] =  2*a.xy()*b.xy();
	c.d[3] =  2*a.xz()*b.xz();
	c.d[6] =  a.xx()*b.xy() + a.xy()*b.xx();
	c.d[10] = a.xx()*b.xz() + a.xz()*b.xx();
	c.d[15] = a.xy()*b.xz() + a.xz()*b.xy();
	
	c.d[2] =  2*a.yy()*b.yy();
	c.d[4] =  2*a.yz()*b.yz();
	c.d[7] =  a.xy()*b.yy() + a.yy()*b.xy();
	c.d[11] = a.xy()*b.yz() + a.yz()*b.xy();
	c.d[16] = a.yy()*b.yz() + a.yz()*b.yy();
	
	c.d[5] =  2*a.zz()*b.zz();
	c.d[8] =  a.xz()*b.yz() + a.yz()*b.xz();
	c.d[12] = a.xz()*b.zz() + a.zz()*b.xz();
	c.d[17] = a.yz()*b.zz() + a.zz()*b.yz();
	
	c.d[9] =  (a.xx()*b.yy() + b.xx()*a.yy())/2. + a.xy()*b.xy();
	c.d[13] = (a.xx()*b.yz() + a.xz()*b.xy() + b.xx()*a.yz() + b.xz()*a.xy())/2.;
	c.d[18] = (a.xy()*b.yz() + a.xz()*b.yy() + b.xy()*a.yz() + b.xz()*a.yy())/2.;
	
	c.d[14] = (a.xx()*b.zz() + b.xx()*a.zz())/2. + a.xz()*b.xz();
	c.d[19] = (a.xy()*b.zz() + a.xz()*b.yz() + b.xy()*a.zz() + b.xz()*a.yz())/2.;
	
	c.d[20] = (a.yy()*b.zz() + b.yy()*a.zz())/2. + a.yz()*b.yz();
	return c;
	
}
