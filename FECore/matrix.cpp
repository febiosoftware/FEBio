// matrix.cpp: implementation of the matrix class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "matrix.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

void lubksb(double**a, int n, int *indx, double b[]);
void ludcmp(double**a, int n, int* indx);

vector<double> operator / (vector<double>& b, matrix& m)
{
	int n = b.size();

	vector<double> x(b);
	vector<int> indx(n);
	matrix a(m);

	ludcmp(a, n, indx);
	lubksb(a, n, indx, x);

	return x;
}
