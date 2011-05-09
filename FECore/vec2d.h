#pragma once

class vec2d
{
public:
	vec2d() { r[0] = r[1] = 0; }
	vec2d(double x, double y) { r[0] = x; r[1] = y; }

	double operator [] (int i) const { return r[i]; }
	double& operator [] (int i) { return r[i]; }

public:
	double r[2];
};
