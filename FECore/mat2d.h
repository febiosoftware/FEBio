#pragma once

#include "vec2d.h"

class mat2d
{
public:
	// constructors
	mat2d(){}
	mat2d(double a00, double a01, double a10, double a11)
	{
		d[0][0] = a00; d[0][1] = a01;
		d[1][0] = a10; d[1][1] = a11;
	}

	// access operators
	double& operator () (int i, int j) { return d[i][j]; }
	double operator () (int i, int j) const { return d[i][j]; }
	double* operator [] (int i) { return d[i]; }

public: // arithmetic operations
	mat2d operator + (const mat2d& m) { return mat2d(d[0][0]+m.d[0][0], d[0][1]+m.d[0][1], d[1][0]+m.d[1][0], d[1][1]+m.d[1][1]); }
	mat2d operator - (const mat2d& m) { return mat2d(d[0][0]-m.d[0][0], d[0][1]-m.d[0][1], d[1][0]-m.d[1][0], d[1][1]-m.d[1][1]); }
	mat2d operator * (double g) { return mat2d(d[0][0]*g, d[0][1]*g, d[1][0]*g, d[1][1]*g); }
	mat2d operator / (double g) { return mat2d(d[0][0]/g, d[0][1]/g, d[1][0]/g, d[1][1]/g); }

	mat2d& operator += (const mat2d& m) { d[0][0] += m.d[0][0]; d[0][1] += m.d[0][1]; d[1][0] += m.d[1][0]; d[1][1] += m.d[1][1]; return *this; }
	mat2d& operator -= (const mat2d& m) { d[0][0] -= m.d[0][0]; d[0][1] -= m.d[0][1]; d[1][0] -= m.d[1][0]; d[1][1] -= m.d[1][1]; return *this; }
	mat2d& operator *= (double g) { d[0][0] *= g; d[0][1] *= g; d[1][0] *= g; d[1][1] *= g; return *this; }
	mat2d& operator /= (double g) { d[0][0] /= g; d[0][1] /= g; d[1][0] /= g; d[1][1] /= g; return *this; }

	mat2d operator * (const mat2d& m) { 
		return mat2d(
			d[0][0]*m.d[0][0]+d[0][1]*m.d[1][0],
			d[0][0]*m.d[0][1]+d[0][1]*m.d[1][1],
			d[1][0]*m.d[0][0]+d[1][1]*m.d[1][0],
			d[1][0]*m.d[0][1]+d[1][1]*m.d[1][1]); 
	}

public:	// matrix operations
	mat2d inverse() const
	{
		double Di = 1/(d[0][0]*d[1][1] - d[0][1]*d[1][0]);
		return mat2d(d[1][1]*Di, -d[0][1]*Di, -d[1][0]*Di, d[0][0]*Di);
	}

	mat2d transpose() const
	{
		return mat2d(d[0][0], d[1][0], d[0][1], d[1][1]);
	}

	void zero()
	{
		d[0][0] = d[0][1] = d[1][0] = d[1][1] = 0.0;
	}

	void identity()
	{
		d[0][0] = d[1][1] = 1.0;
		d[0][1] = d[1][0] = 0.0;
	}
	
protected:
	double	d[2][2];
};

// matrix-vector operations
inline vec2d operator * (mat2d& m, vec2d& a) { return vec2d(m[0][0]*a[0]+m[0][1]*a[1], m[1][0]*a[0]+m[1][1]*a[1]); }

// dyadic product
inline mat2d dyad(vec2d& a, vec2d& b) { return mat2d(a.r[0]*b.r[0], a.r[0]*b.r[1], a.r[1]*b.r[0], a.r[1]*b.r[1]); }
