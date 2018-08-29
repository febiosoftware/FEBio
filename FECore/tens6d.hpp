// NOTE: This file is automatically included from tens5d.h
// Users should not include this file manually!

inline double tens6d::operator () (int i, int j, int k, int l, int m, int n) const
{
	int R = 3*(3*i + j) + k;
	int C = 3*(3*l + m) + n;
	return d[27*C + R];
}

inline double& tens6d::operator () (int i, int j, int k, int l, int m, int n)
{
	int R = 3*(3*i + j) + k;
	int C = 3*(3*l + m) + n;
	return d[27*C + R];
}
