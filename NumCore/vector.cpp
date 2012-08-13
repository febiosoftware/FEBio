#include "stdafx.h"
#include <assert.h>
#include "vector.h"
#include <algorithm>

double operator*(const vector<double>& a, const vector<double>& b)
{
	double sum = 0;
	for (size_t i=0; i<a.size(); i++) sum += a[i]*b[i];

	return sum;
}

vector<double> operator - (vector<double>& a, vector<double>& b)
{
	vector<double> c(a);
	int n = (int) c.size();
	for (int i=0; i<n; ++i) c[i] -= b[i];
	return c;
}

vector<double>& operator += (vector<double>& a, const vector<double>& b)
{
	assert(a.size() == b.size());
	for (size_t i = 0; i < a.size(); ++i) a[i] += b[i];
	return a;
}

vector<double>& operator *= (vector<double>& a, double b)
{
	for (size_t i=0; i<a.size(); ++i) a[i] *= b;
	return a;
}

void vcopys(vector<double>& a, const vector<double>& b, double s)
{
	a = b;
	a *= s;
}

vector<double> operator + (const vector<double>& a, const vector<double>& b)
{
	assert(a.size() == b.size());
	vector<double> s(a);
	for (size_t i = 0; i < s.size(); ++i) s[i] += b[i];
	return s;
}
