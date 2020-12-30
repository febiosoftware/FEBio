/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2020 University of Utah, The Trustees of Columbia University in
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
	m_nnc = 8;
	m_checkForMatch = false;
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

	for (int i = 0; i < N1; ++i)
	{
		vec3d ri = m_trg[i];
		int M = findNeirestNeighbors(m_src, ri, m_nnc, m_data[i].cpl);
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
		d.A.resize(4, 4);
		d.A.zero();
		for (int m = 0; m < M; ++m)
		{
			vec3d ri = d.X[m];
			double P[4] = { 1.0, ri.x, ri.y, ri.z };
			for (int a = 0; a < 4; ++a)
			{
				for (int b = 0; b < 4; ++b)
				{
					d.A(a, b) += d.W[m] * P[a] * P[b];
				}
			}
		}

		// solve the linear system of equations
		d.index.resize(4);
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
		vector<double> b(4, 0.0);
		for (int m = 0; m < M; ++m)
		{
			vec3d ri = d.X[m];
			double P[4] = { 1.0, ri.x, ri.y, ri.z };

			double vm = src(closestNodes[m]);

			for (int a = 0; a < 4; ++a)
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
		if (W[0] > 0.999)
		{
			return f(closestNodes[0]);
		}
	}

	// update nodal values
	vector<double> b(4, 0.0);
	for (int m = 0; m < M; ++m)
	{
		vec3d ri = d.X[m];
		double P[4] = { 1.0, ri.x, ri.y, ri.z };

		double vm = f(closestNodes[m]);

		for (int a = 0; a < 4; ++a)
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
		if (W[0] > 0.999)
		{
			return f(closestNodes[0]);
		}
	}

	// update nodal values
	vector<double> bx(4, 0.0), by(4, 0.0), bz(4, 0.0);
	for (int m = 0; m < M; ++m)
	{
		vec3d ri = d.X[m];
		double P[4] = { 1.0, ri.x, ri.y, ri.z };

		vec3d vm = f(closestNodes[m]);

		for (int a = 0; a < 4; ++a)
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
