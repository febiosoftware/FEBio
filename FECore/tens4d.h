#pragma once

#include "mat3d.h"

//-----------------------------------------------------------------------------
//! Class for 4th order tensors with major and minor symmetries (i.e., super-symmetry)

// Due to the major symmetry we can store this tensor as a 6x6 matrix.
// The tensor is stored in column major order:
//
//     / 0   1   3   6   10   15  \
//     |     2   4   7   11   16  |
//     |         5   8   12   17  |
// A = |             9   13   18  |
//     |                 14   19  |
//     \                      20  /
//
// Note that due to the minor symmetry we only store the upper matrix of this tensor
//

class tens4ds
{
public:
	enum { NNZ = 21 };

	// default constructor
	tens4ds(){zero();}
	tens4ds(const double g)
	{
		d[ 0] = g;
		d[ 1] = g; d[ 2] = g;
		d[ 3] = g; d[ 4] = g; d[ 5] = g;
		d[ 6] = g; d[ 7] = g; d[ 8] = g; d[ 9] = g;
		d[10] = g; d[11] = g; d[12] = g; d[13] = g; d[14] = g;
		d[15] = g; d[16] = g; d[17] = g; d[18] = g; d[19] = g; d[20] = g;
	}

	tens4ds(double m[6][6])
	{
		d[ 0] = m[0][0];
		d[ 1] = m[0][1]; d[ 2] = m[1][1];
		d[ 3] = m[0][2]; d[ 4] = m[1][2]; d[ 5] = m[2][2];
		d[ 6] = m[0][3]; d[ 7] = m[1][3]; d[ 8] = m[2][3]; d[ 9] = m[3][3];
		d[10] = m[0][4]; d[11] = m[1][4]; d[12] = m[2][4]; d[13] = m[3][4]; d[14] = m[4][4];
		d[15] = m[0][5]; d[16] = m[1][5]; d[17] = m[2][5]; d[18] = m[3][5]; d[19] = m[4][5]; d[20] = m[5][5];
	}

	double& operator () (int i, int j, int k, int l)
	{
		const int m[3][3] = {{0,3,5},{3,1,4},{5,4,2}};
		tens4ds& T = (*this);
		return T(m[i][j], m[k][l]);
	}

	double operator () (int i, int j, int k, int l) const
	{
		const int m[3][3] = {{0,3,5},{3,1,4},{5,4,2}};
		const tens4ds& T = (*this);
		return T(m[i][j], m[k][l]);
	}

	double& operator () (int i, int j)
	{
		const int m[6] = {0, 1, 3, 6, 10, 15};
		if (i<=j) return d[m[j]+i]; else return d[m[i]+j];
	}

	double operator () (int i, int j) const
	{
		const int m[6] = {0, 1, 3, 6, 10, 15};
		if (i<=j) return d[m[j]+i]; else return d[m[i]+j];
	}

	// arithmetic operators
	tens4ds operator + (const tens4ds& t) const;
	tens4ds operator - (const tens4ds& t) const;
	tens4ds operator * (double g) const;
	tens4ds operator / (double g) const;

	// arithmetic assignment operators
	tens4ds& operator += (const tens4ds& t);
	tens4ds& operator -= (const tens4ds& t);
	tens4ds& operator *= (double g);
	tens4ds& operator /= (double g);

	// unary operators
	tens4ds operator - () const;
	
	// double dot product with tensor
	mat3ds dot(const mat3ds& m) const;

	// trace
	double tr() const;
	
	// initialize to zero
	void zero();

	// extract 6x6 matrix
	void extract(double d[6][6]);

	// calculates the inverse
	tens4ds inverse() const;

public:
	double d[NNZ];	// stored in column major order
};

//! Check positive definiteness of a 4th-order symmetric tensor
bool IsPositiveDefinite(const tens4ds& t);

// outer (dyadic) products for symmetric matrices
tens4ds dyad1s(const mat3ds& a);
tens4ds dyad1s(const mat3ds& a, const mat3ds& b);
tens4ds dyad2s(const mat3ds& a);
tens4ds dyad2s(const mat3ds& a, const mat3ds& b);
tens4ds dyad4s(const mat3ds& a);
tens4ds dyad4s(const mat3ds& a, const mat3ds& b);
tens4ds ddots(const tens4ds& a, const tens4ds& b);
mat3d vdotTdotv(const vec3d a, const tens4ds T, const vec3d b);

inline tens4ds operator * (const double g, const tens4ds& a) { return a*g; }

// The following file contains the actual definition of the class functions
#include "tens4ds.hpp"



// LTE - New 4th order tensor class with major symmetry only
//-----------------------------------------------------------------------------
//! Class for 4th order tensors with major symmetry only

// Due to the lack of minor symmetry, we have to store additional components of this tensor as a 9x9 matrix.
// Major symmetry ensures that this storage matrix is symmetric about its main diagonal.
// The tensor is stored in column major order:
//
//     / 0   1   3   6   10   15  21   28   36  \
//     |     2   4   7   11   16  22   29   37  |
//     |         5   8   12   17  23   30   38  |
// A = |             9   13   18  24   31   39  |
//     |                 14   19  25   32   40  |
//     |                      20  26   33   41  |
//     |                          27   34   42  |
//     |                               35   43  |
//     \                                    44  /
//


class tens4dms
{
// Declare and define fields:
public:
	enum { NNZ = 45 };

	// Declare and define functions:
	
	// Default constructor
	tens4dms(){zero();}
	
	tens4dms(const double g)
	{
		for (int i = 0; i < NNZ; i++)
			d[i] = g;
	}

	tens4dms(double m[9][9])
	{
		d[ 0] = m[0][0];
		d[ 1] = m[0][1]; d[ 2] = m[1][1];
		d[ 3] = m[0][2]; d[ 4] = m[1][2]; d[ 5] = m[2][2];
		d[ 6] = m[0][3]; d[ 7] = m[1][3]; d[ 8] = m[2][3]; d[ 9] = m[3][3];
		d[10] = m[0][4]; d[11] = m[1][4]; d[12] = m[2][4]; d[13] = m[3][4]; d[14] = m[4][4];
		d[15] = m[0][5]; d[16] = m[1][5]; d[17] = m[2][5]; d[18] = m[3][5]; d[19] = m[4][5]; d[20] = m[5][5];
		d[21] = m[0][6]; d[22] = m[1][6]; d[23] = m[2][6]; d[24] = m[3][6]; d[25] = m[4][6]; d[26] = m[5][6]; d[27] = m[6][6];
		d[28] = m[0][7]; d[29] = m[1][7]; d[30] = m[2][7]; d[31] = m[3][7]; d[32] = m[4][7]; d[33] = m[5][7]; d[34] = m[6][7]; d[35] = m[7][7];
		d[36] = m[0][8]; d[37] = m[1][8]; d[38] = m[2][8]; d[39] = m[3][8]; d[40] = m[4][8]; d[41] = m[5][8]; d[42] = m[6][8]; d[43] = m[7][8]; d[44] = m[8][8];
	}

	double& operator () (int i, int j, int k, int l)
	{
		const int m[3][3] = {{0,3,5},{6,1,4},{8,7,2}};
		tens4dms& T = (*this);
		return T(m[i][j], m[k][l]);
	}

	double operator () (int i, int j, int k, int l) const
	{
		const int m[3][3] = {{0,3,5},{6,1,4},{8,7,2}};
		const tens4dms& T = (*this);
		return T(m[i][j], m[k][l]);
	}

	double& operator () (int i, int j)
	{
		const int m[9] = {0, 1, 3, 6, 10, 15, 21, 28, 36};
		if (i<=j) return d[m[j]+i]; else return d[m[i]+j];
	}

	double operator () (int i, int j) const
	{
		const int m[9] = {0, 1, 3, 6, 10, 15, 21, 28, 36};
		if (i<=j) return d[m[j]+i]; else return d[m[i]+j];
	}

	// arithmetic operators
	tens4dms operator + (const tens4dms& t) const;
	tens4dms operator - (const tens4dms& t) const;
	tens4dms operator * (double g) const;
	tens4dms operator / (double g) const;

	// arithmetic assignment operators
	tens4dms& operator += (const tens4dms& t);
	tens4dms& operator -= (const tens4dms& t);
	tens4dms& operator *= (double g);
	tens4dms& operator /= (double g);

	// unary operators
	tens4dms operator - () const;
	
	//// double dot product with tensor
	//mat3ds dot(const mat3ds& m) const;

	// trace
	double tr() const;
	
	// initialize to zero
	void zero();

	// extract 9x9 matrix
	void extract(double d[9][9]);

	// compute the super-symmetric (major and minor symmetric) component of the tensor
	tens4ds supersymm() const;

	//// calculates the inverse
	//tens4dms inverse() const;

public:
	double d[NNZ];	// stored in column major order
};

////! Check positive definiteness of a 4th-order symmetric tensor
//bool IsPositiveDefinite(const tens4ds& t);
//
// outer (dyadic) products for symmetric and non-symmetric matrices
tens4dms dyad1(const mat3d& a);
tens4dms dyad1(const mat3ds& a, const mat3ds& b);
//tens4dms dyad2s(const mat3ds& a);
//tens4dms dyad2s(const mat3ds& a, const mat3ds& b);
//tens4ds ddots(const tens4ds& a, const tens4ds& b);
//mat3d vdotTdotv(const vec3d a, const tens4ds T, const vec3d b);

//inline tens4ds operator * (const double g, const tens4ds& a) { return a*g; }

// The following file contains the actual definition of the class functions
#include "tens4dms.hpp"


// LTE - Generalized 4th order tensor class without symmetry
//-----------------------------------------------------------------------------
//! Class for 4th order tensors without symmetry

// We store components of this tensor as a 9x9 matrix.
// The tensor is stored in column major order:
//
//     / 0   9   18   27   36   45   54   63   72  \
//     | 1   10  19   28   37   46   55   64   73  |
//     | 2   11  20   29   38   47   56   65   74  |
// A = | 3   12  21   30   39   48   57   66   75  |
//     | 4   13  22   31   40   49   58   67   76  |
//     | 5   14  23   32   41   50   59   68   77  |
//     | 6   15  24   33   42   51   60   69   78  |
//     | 7   16  25   34   43   52   61   70   79  |
//     \ 8   17  26   35   44   53   62   71   80  /
//


class tens4d
{
// Declare and define fields:
public:
	enum { NNZ = 81 };

	// Declare and define functions:
	
	// Default constructor
	tens4d(){zero();}
	
	tens4d(const double g)
	{
		for (int i = 0; i < NNZ; i++)
			d[i] = g;
	}

	tens4d(double m[9][9])
	{
		d[ 0] = m[0][0]; d[ 1] = m[1][0]; d[ 2] = m[2][0]; d[ 3] = m[3][0]; d[ 4] = m[4][0]; d[ 5] = m[5][0]; d[ 6] = m[6][0]; d[ 7] = m[7][0]; d[ 8] = m[8][0];
		d[ 9] = m[0][1]; d[10] = m[1][1]; d[11] = m[2][1]; d[12] = m[3][1]; d[13] = m[4][1]; d[14] = m[5][1]; d[15] = m[6][1]; d[16] = m[7][1]; d[17] = m[8][1];
		d[18] = m[0][2]; d[19] = m[1][2]; d[20] = m[2][2]; d[21] = m[3][2]; d[22] = m[4][2]; d[23] = m[5][2]; d[24] = m[6][2]; d[25] = m[7][2]; d[26] = m[8][2];
		d[27] = m[0][3]; d[28] = m[1][3]; d[29] = m[2][3]; d[30] = m[3][3]; d[31] = m[4][3]; d[32] = m[5][3]; d[33] = m[6][3]; d[34] = m[7][3]; d[35] = m[8][3];
		d[36] = m[0][4]; d[37] = m[1][4]; d[38] = m[2][4]; d[39] = m[3][4]; d[40] = m[4][4]; d[41] = m[5][4]; d[42] = m[6][4]; d[43] = m[7][4]; d[44] = m[8][4];
		d[45] = m[0][5]; d[46] = m[1][5]; d[47] = m[2][5]; d[48] = m[3][5]; d[49] = m[4][5]; d[50] = m[5][5]; d[51] = m[6][5]; d[52] = m[7][5]; d[53] = m[8][5];
		d[54] = m[0][6]; d[55] = m[1][6]; d[56] = m[2][6]; d[57] = m[3][6]; d[58] = m[4][6]; d[59] = m[5][6]; d[60] = m[6][6]; d[61] = m[7][6]; d[62] = m[8][6];
		d[63] = m[0][7]; d[64] = m[1][7]; d[65] = m[2][7]; d[66] = m[3][7]; d[67] = m[4][7]; d[68] = m[5][7]; d[69] = m[6][7]; d[70] = m[7][7]; d[71] = m[8][7];
		d[72] = m[0][8]; d[73] = m[1][8]; d[74] = m[2][8]; d[75] = m[3][8]; d[76] = m[4][8]; d[77] = m[5][8]; d[78] = m[6][8]; d[79] = m[7][8]; d[80] = m[8][8];
	}

	double& operator () (int i, int j, int k, int l)
	{
		const int m[3][3] = {{0,3,5},{6,1,4},{8,7,2}};
		tens4d& T = (*this);
		return T(m[i][j], m[k][l]);
	}

	double operator () (int i, int j, int k, int l) const
	{
		const int m[3][3] = {{0,3,5},{6,1,4},{8,7,2}};
		const tens4d& T = (*this);
		return T(m[i][j], m[k][l]);
	}

	double& operator () (int i, int j)
	{
		const int m[9] = {0, 9, 19, 27, 36, 45, 54, 63, 72};
		return d[m[j]+i];
	}

	double operator () (int i, int j) const
	{
		const int m[9] = {0, 9, 19, 27, 36, 45, 54, 63, 72};
		return d[m[j]+i];
	}

	// arithmetic operators
	tens4d operator + (const tens4d& t) const;
	tens4d operator - (const tens4d& t) const;
	tens4d operator * (double g) const;
	tens4d operator / (double g) const;

	// arithmetic assignment operators
	tens4d& operator += (const tens4d& t);
	tens4d& operator -= (const tens4d& t);
	tens4d& operator *= (double g);
	tens4d& operator /= (double g);

	// unary operators
	tens4d operator - () const;
	
	//// double dot product with tensor
//	//mat3ds dot(const mat3ds& m) const;
//
	// trace
	double tr() const;
	
	// initialize to zero
	void zero();

	// extract 9x9 matrix
	void extract(double d[9][9]);

	// compute the super-symmetric (major and minor symmetric) component of the tensor
	tens4ds supersymm() const;
//
//	//// calculates the inverse
//	//tens4dms inverse() const;
//
public:
	double d[NNZ];	// stored in column major order
};
//
//////! Check positive definiteness of a 4th-order symmetric tensor
////bool IsPositiveDefinite(const tens4ds& t);
////
//// outer (dyadic) products for symmetric and non-symmetric matrices
//tens4dms dyad1(const mat3d& a);
//tens4dms dyad1(const mat3ds& a, const mat3ds& b);
//tens4dms dyad2s(const mat3ds& a);
//tens4dms dyad2s(const mat3ds& a, const mat3ds& b);
////tens4ds ddots(const tens4ds& a, const tens4ds& b);
////mat3d vdotTdotv(const vec3d a, const tens4ds T, const vec3d b);
//
////inline tens4ds operator * (const double g, const tens4ds& a) { return a*g; }
//
// The following file contains the actual definition of the class functions
#include "tens4d.hpp"
