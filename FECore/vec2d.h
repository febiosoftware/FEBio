#pragma once

class vec2d
{
public:
	// constructor
	vec2d() { r[0] = r[1] = 0; }
	vec2d(double x, double y) { r[0] = x; r[1] = y; }

	// access operators
	double operator [] (int i) const { return r[i]; }
	double& operator [] (int i) { return r[i]; }

	double& x() { return r[0]; }
	double& y() { return r[1]; }

public: // arithmetic operators

	vec2d operator + (const vec2d& v) { return vec2d(r[0]+v.r[0], r[1]+v.r[1]); }
	vec2d operator - (const vec2d& v) { return vec2d(r[0]-v.r[0], r[1]-v.r[1]); }
	vec2d operator * (double g) { return vec2d(r[0]*g, r[1]*g); }
	vec2d operator / (double g) { return vec2d(r[0]/g, r[1]/g); }

	vec2d& operator += (const vec2d& v) { r[0] += v.r[0]; r[1] += v.r[1]; return *this; }
	vec2d& operator -= (const vec2d& v) { r[0] -= v.r[0]; r[1] -= v.r[1]; return *this; }
	vec2d& operator *= (double g) { r[0] *= g; r[1] *= g; return *this; }
	vec2d& operator /= (double g) { r[0] /= g; r[1] /= g; return *this; }

    vec2d operator - () { return vec2d(-r[0], -r[1]); }
    
	// dot product
	double operator * (const vec2d& v) const { return r[0]*v[0] + r[1]*v[1]; }

public:
	double r[2];
};
