// NOTE: This file is automatically included from tens4d.h
// Users should not include this file manually!

inline double& tens4d::operator () (int i, int j, int k, int l)
{
	const int m[3][3] = {{0,3,5},{6,1,4},{8,7,2}};
	tens4d& T = (*this);
	return T(m[i][j], m[k][l]);
}

inline double tens4d::operator () (int i, int j, int k, int l) const
{
	const int m[3][3] = {{0,3,5},{6,1,4},{8,7,2}};
	const tens4d& T = (*this);
	return T(m[i][j], m[k][l]);
}

inline double& tens4d::operator () (int i, int j)
{
	const int m[9] = {0, 9, 18, 27, 36, 45, 54, 63, 72};
	return d[m[j]+i];
}

inline double tens4d::operator () (int i, int j) const
{
	const int m[9] = {0, 9, 18, 27, 36, 45, 54, 63, 72};
	return d[m[j]+i];
}
