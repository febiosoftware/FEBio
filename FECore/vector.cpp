#include "stdafx.h"
#include <assert.h>
#include "vector.h"
#include "FEMesh.h"
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

void operator += (vector<double>& a, const vector<double>& b)
{
	assert(a.size() == b.size());
	for (size_t i = 0; i < a.size(); ++i) a[i] += b[i];
}

void operator -= (vector<double>& a, const vector<double>& b)
{
	assert(a.size() == b.size());
	for (size_t i = 0; i < a.size(); ++i) a[i] -= b[i];
}

void operator *= (vector<double>& a, double b)
{
	for (size_t i=0; i<a.size(); ++i) a[i] *= b;
}

void vcopys(vector<double>& a, const vector<double>& b, double s)
{
	assert(a.size() == b.size());
	for (size_t i=0; i<a.size(); ++i) a[i] = b[i]*s;
}

void vadds(vector<double>& a, const vector<double>& b, double s)
{
	assert(a.size() == b.size());
	for (size_t i = 0; i<a.size(); ++i) a[i] += b[i] * s;
}

void vsubs(vector<double>& a, const vector<double>& b, double s)
{
	assert(a.size() == b.size());
	for (size_t i = 0; i<a.size(); ++i) a[i] -= b[i] * s;
}

void vscale(vector<double>& a, const vector<double>& s)
{
	assert(a.size() == s.size());
	for (size_t i = 0; i<a.size(); ++i) a[i] *= s[i];
}

void vsub(vector<double>& a, const vector<double>& l, const vector<double>& r)
{
	assert((a.size()==l.size())&&(a.size()==r.size()));
	for (size_t i=0; i<a.size(); ++i) a[i] = l[i] - r[i];
}

vector<double> operator + (const vector<double>& a, const vector<double>& b)
{
	assert(a.size() == b.size());
	vector<double> s(a);
	for (size_t i = 0; i < s.size(); ++i) s[i] += b[i];
	return s;
}

void gather(vector<double>& v, FEMesh& mesh, int ndof)
{
	const int NN = mesh.Nodes();
	for (int i=0; i<NN; ++i)
	{
		FENode& node = mesh.Node(i);
		int n = node.m_ID[ndof]; if (n >= 0) v[n] = node.get(ndof);
	}
}

void gather(vector<double>& v, FEMesh& mesh, const vector<int>& dof)
{
	const int NN = mesh.Nodes();
	const int NDOF = (const int) dof.size();
	for (int i=0; i<NN; ++i)
	{
		FENode& node = mesh.Node(i);
		for (int j=0; j<NDOF; ++j)
		{
			int n = node.m_ID[dof[j]]; 
			if (n >= 0) v[n] = node.get(dof[j]);
		}
	}
}

void scatter(vector<double>& v, FEMesh& mesh, int ndof)
{
	const int NN = mesh.Nodes();
	for (int i=0; i<NN; ++i)
	{
		FENode& node = mesh.Node(i);
		int n = node.m_ID[ndof];
		if (n >= 0) node.set(ndof, v[n]);
	}
}

double l2_norm(const vector<double>& v)
{
	double s = 0.0;
	for (auto vi : v) s += vi*vi;
	return sqrt(s);
}
