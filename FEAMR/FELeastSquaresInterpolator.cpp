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
#include "FELeastSquaresInterpolator.h"
#include <FECore/FENNQuery.h>
#include <algorithm>
using namespace std;

class KDTree
{
public:
	KDTree()
	{
		m_parent = nullptr;
		m_left = nullptr;
		m_right = nullptr;
	}

	void build(vector<vec3d> pts, int depth = 0)
	{
		size_t n = pts.size();
		if (n == 1)
		{
			m_r = pts[0];
			return;
		}

		int axis = depth % 3;
		std::sort(pts.begin(), pts.end(), [=](const vec3d& a, const vec3d& b) {
			if (axis == 0) return (a.x < b.x);
			if (axis == 1) return (a.y < b.y);
			if (axis == 2) return (a.z < b.z);
			return false;
		});

		int med = n / 2;

		m_r = pts[med];

		// build left list
		if (med > 0)
		{
			vector<vec3d> l(pts.begin(), pts.begin() + med);
			m_left = new KDTree(this);
			m_left->build(l, depth + 1);
		}

		// build right list
		if (med < n - 1)
		{
			vector<vec3d> r(pts.begin() + med + 1, pts.end());
			m_right = new KDTree(this);
			m_right->build(r, depth + 1);
		}
	}

public:
	vec3d	m_r;
	KDTree*	m_parent;
	KDTree*	m_left;
	KDTree*	m_right;

public:
	KDTree(KDTree* parent) : m_parent(parent)
	{
		m_left = nullptr;
		m_right = nullptr;
	}
};

class NearestNeighborSearch
{
public:
	NearestNeighborSearch() {}

	void Init(const std::vector<vec3d>& points, int k)
	{
		m_k = k;
		m_points = points;

//		m_kdtree.build(m_points);
	}

	int findNearestNeighbors(const vec3d& x, std::vector<int>& closestNodes)
	{
		return findNeirestNeighbors(m_points, x, m_k, closestNodes);
	}

protected:
	int				m_k;
	vector<vec3d>	m_points;

	KDTree	m_kdtree;
};

FELeastSquaresInterpolator::Data::Data() {}
FELeastSquaresInterpolator::Data::Data(const Data& d)
{
	A = d.A;
	index = d.index;
	W = d.W;
	X = d.X;
	cpl = d.cpl;
}
void FELeastSquaresInterpolator::Data::operator = (const Data& d)
{
	A = d.A;
	index = d.index;
	W = d.W;
	X = d.X;
	cpl = d.cpl;
}

FELeastSquaresInterpolator::FELeastSquaresInterpolator()
{
	m_dim = 3;
	m_nnc = 8;
	m_checkForMatch = false;
}

//! Set dimension (2 or 3)
void FELeastSquaresInterpolator::SetDimension(int d)
{
	assert((d == 2) || (d == 3));
	m_dim = d;
}

void FELeastSquaresInterpolator::SetNearestNeighborCount(int nnc) { m_nnc = nnc; }

void FELeastSquaresInterpolator::SetCheckForMatch(bool b) { m_checkForMatch = b; }

void FELeastSquaresInterpolator::SetSourcePoints(const vector<vec3d>& srcPoints)
{
	m_src = srcPoints;
}

void FELeastSquaresInterpolator::SetTargetPoints(const vector<vec3d>& trgPoints)
{
	m_trg = trgPoints;
}

bool FELeastSquaresInterpolator::SetTargetPoint(const vec3d& trgPoint)
{
	m_trg.clear();
	m_trg.push_back(trgPoint);
	return Init();
}

bool FELeastSquaresInterpolator::Init()
{
	if (m_nnc < 4) return false;
	if (m_src.empty()) return false;
	if (m_trg.empty()) return false;

	int N0 = m_src.size();
	int N1 = m_trg.size();

	m_data.resize(N1);

	// initialize nearest neighbor search
	NearestNeighborSearch NNS;
	NNS.Init(m_src, m_nnc);

	// do nearest-neighbor search
	for (int i = 0; i < N1; ++i)
	{
		vec3d ri = m_trg[i];
		int M = NNS.findNearestNeighbors(ri, m_data[i].cpl);
		assert(M > 4);
		m_data[i].cpl.resize(M);
	}

	for (int i = 0; i < N1; ++i)
	{
		Data& d = m_data[i];
		vec3d x = m_trg[i];

		vector<int>& closestNodes = m_data[i].cpl;
		int M = closestNodes.size();

		// the last node is the farthest and determines the radius
		vec3d& r = m_src[closestNodes[M - 1]];
		double L = sqrt((r - x)*(r - x));

		// add some offset to make sure none of the points will have a weight of zero.
		// Such points would otherwise be ignored, which reduces the net number of interpolation
		// points and make the MLS system ill-conditioned
		L += L * 0.05;

		// evaluate weights and displacements
		d.X.resize(M);
		d.W.resize(M);
		for (int m = 0; m < M; ++m)
		{
			vec3d rm = m_src[closestNodes[m]];
			vec3d rj = x - rm;

			d.X[m] = rj;

			double D = sqrt(rj*rj);
			double wj = 1.0 - D / L;
			d.W[m] = wj;
		}

		// setup least squares problems
		d.A.resize(m_dim + 1, m_dim + 1);
		d.A.zero();
		for (int m = 0; m < M; ++m)
		{
			vec3d ri = d.X[m];
			double P[4] = { 1.0, ri.x, ri.y, ri.z };
			for (int a = 0; a <= m_dim; ++a)
			{
				for (int b = 0; b <= m_dim; ++b)
				{
					d.A(a, b) += d.W[m] * P[a] * P[b];
				}
			}
		}

		// solve the linear system of equations
		d.index.resize(m_dim + 1);
		d.A.lufactor(d.index);
	}

	return true;
}

bool FELeastSquaresInterpolator::Map(std::vector<double>& tval, function<double(int sourceNode)> src)
{
	if (m_data.size() != m_trg.size()) return false;

	for (int i = 0; i < m_trg.size(); ++i)
	{
		Data& d = m_data[i];

		vector<int>& closestNodes = m_data[i].cpl;
		int M = closestNodes.size();

		// evaluate weights and positions
		vector<vec3d>& X = d.X;
		vector<double>& W = d.W;

		// update nodal values
		vector<double> b(m_dim + 1, 0.0);
		for (int m = 0; m < M; ++m)
		{
			vec3d ri = d.X[m];
			double P[4] = { 1.0, ri.x, ri.y, ri.z };

			double vm = src(closestNodes[m]);

			for (int a = 0; a <= m_dim; ++a)
			{
				b[a] += W[m] * P[a] * vm;
			}
		}

		// solve the linear system of equations
		d.A.lusolve(b, d.index);

		tval[i] = b[0];
	}

	return true;
}

double FELeastSquaresInterpolator::Map(int inode, function<double(int sourceNode)> f)
{
	Data& d = m_data[inode];

	vector<int>& closestNodes = m_data[inode].cpl;
	int M = closestNodes.size();

	// evaluate weights and positions
	vector<vec3d>& X = d.X;
	vector<double>& W = d.W;

	if (m_checkForMatch)
	{
		if (W[0] > 0.9999)
		{
			return f(closestNodes[0]);
		}
	}

	// update nodal values
	vector<double> b(m_dim + 1, 0.0);
	for (int m = 0; m < M; ++m)
	{
		vec3d ri = d.X[m];
		double P[4] = { 1.0, ri.x, ri.y, ri.z };

		double vm = f(closestNodes[m]);

		for (int a = 0; a <= m_dim; ++a)
		{
			b[a] += W[m] * P[a] * vm;
		}
	}

	// solve the linear system of equations
	d.A.lusolve(b, d.index);

	return b[0];
}

vec3d FELeastSquaresInterpolator::MapVec3d(int inode, function<vec3d(int sourceNode)> f)
{
	Data& d = m_data[inode];

	vector<int>& closestNodes = m_data[inode].cpl;
	int M = closestNodes.size();

	// evaluate weights and positions
	vector<vec3d>& X = d.X;
	vector<double>& W = d.W;

	if (m_checkForMatch)
	{
		if (W[0] > 0.9999)
		{
			return f(closestNodes[0]);
		}
	}

	// update nodal values
	vector<double> bx(m_dim+1, 0.0), by(m_dim + 1, 0.0), bz(m_dim + 1, 0.0);
	for (int m = 0; m < M; ++m)
	{
		vec3d ri = d.X[m];
		double P[4] = { 1.0, ri.x, ri.y, ri.z };

		vec3d vm = f(closestNodes[m]);

		for (int a = 0; a <= m_dim; ++a)
		{
			bx[a] += W[m] * P[a] * vm.x;
			by[a] += W[m] * P[a] * vm.y;
			bz[a] += W[m] * P[a] * vm.z;
		}
	}

	// solve the linear system of equations
	d.A.lusolve(bx, d.index);
	d.A.lusolve(by, d.index);

	d.A.lusolve(bz, d.index);

	return vec3d(bx[0], by[0], bz[0]);
}
