/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/



#include "stdafx.h"
#include <assert.h>
#include "vector.h"
#include "FEMesh.h"
#include "FEDofList.h"
#include <algorithm>

double operator*(const vector<double>& a, const vector<double>& b)
{
	double sum_p = 0, sum_n = 0;
	for (size_t i = 0; i < a.size(); i++)
	{
		double ab = a[i] * b[i];
		if (ab >= 0.0) sum_p += ab; else sum_n += ab;
	}
	return sum_p + sum_n;
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

vector<double> operator*(const vector<double>& a, double g)
{
	vector<double> s(a.size());
	for (size_t i = 0; i < s.size(); ++i) s[i] = a[i]*g;
	return s;
}

vector<double> FECORE_API operator - (const vector<double>& a)
{
	vector<double> s(a.size());
	for (size_t i = 0; i < s.size(); ++i) s[i] = -a[i];
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

void scatter3(vector<double>& v, FEMesh& mesh, int ndof1, int ndof2, int ndof3)
{
	const int NN = mesh.Nodes();
#pragma omp parallel for 
	for (int i = 0; i<NN; ++i)
	{
		FENode& node = mesh.Node(i);
		int n;
		n = node.m_ID[ndof1]; if (n >= 0) node.set(ndof1, v[n]);
		n = node.m_ID[ndof2]; if (n >= 0) node.set(ndof2, v[n]);
		n = node.m_ID[ndof3]; if (n >= 0) node.set(ndof3, v[n]);
	}
}

void scatter(vector<double>& v, FEMesh& mesh, const FEDofList& dofs)
{
	const int NN = mesh.Nodes();
	for (int i = 0; i<NN; ++i)
	{
		FENode& node = mesh.Node(i);
		for (int j = 0; j < dofs.Size(); ++j)
		{
			int n = node.m_ID[dofs[j]]; if (n >= 0) node.set(dofs[j], v[n]);
		}
	}
}

double l2_norm(const vector<double>& v)
{
	double s = 0.0;
	for (auto vi : v) s += vi*vi;
	return sqrt(s);
}

double l2_sqrnorm(const vector<double>& v)
{
	double s = 0.0;
	for (auto vi : v) s += vi*vi;
	return s;
}

double l2_norm(double* x, int n)
{
	double s = 0.0;
	for (int i = 0; i < n; ++i) s += x[i]*x[i];
	return sqrt(s);
}
