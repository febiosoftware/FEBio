#pragma once
// NOTE: This file is automatically included from tens6d.h
// Users should not include this file manually!

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
