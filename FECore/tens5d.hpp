#pragma once
// NOTE: This file is automatically included from tens5d.h
// Users should not include this file manually!

// access operator
inline double tens5d::operator () (int i, int j, int k, int l, int m) const
{
	int R = 3*(3*i + j) + k;
	int C = 3*l + m;
	return d[27*C + R];
}

// access operator
inline double& tens5d::operator () (int i, int j, int k, int l, int m)
{
	int R = 3*(3*i + j) + k;
	int C = 3*l + m;
	return d[27*C + R];
}
