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
#include "FENNQuery.h"
#include "FESurface.h"
#include <stdlib.h>
#include "FEMesh.h"
using namespace std;

int cmp_node(const void* e1, const void* e2)
{
	FENNQuery::NODE& n1 = *((FENNQuery::NODE*)e1);
	FENNQuery::NODE& n2 = *((FENNQuery::NODE*)e2);

	return (n1.d1 > n2.d1 ? 1 : -1);
}

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

FENNQuery::FENNQuery(FESurface* ps)
{
	m_ps = ps;
}

FENNQuery::~FENNQuery()
{

}

//-----------------------------------------------------------------------------

void FENNQuery::Init()
{
	assert(m_ps);

	int i;
	vec3d r0, r;
	int N = m_ps->Nodes();

	// pick a random point as pivot
	r0 = m_q1 = m_ps->Node(0).m_rt;

	// find the farthest node of this node
	double dmax = 0, d;
	for (i=0; i<N; ++i)
	{
		r = m_ps->Node(i).m_rt;
		d = (r - r0)*(r - r0);
		if (d > dmax)
		{
			m_q1 = r;
			dmax = d;
		}
	}

	// let's find the farthest node of this node
	r0 = m_q2 = m_q1;
	dmax = 0;
	for (i=0; i<N; ++i)
	{
		r = m_ps->Node(i).m_rt;
		d = (r - r0)*(r - r0);
		if (d > dmax)
		{
			m_q2 = r;
			dmax = d;
		}
	}

	// create the BK-"tree"
	m_bk.resize(N);
	for (i=0; i<N; ++i) 
	{
		r = m_ps->Node(i).m_rt;
		m_bk[i].i = i;
		m_bk[i].r = r;
		m_bk[i].d1 = (m_q1 - r)*(m_q1 - r);
		m_bk[i].d2 = (m_q2 - r)*(m_q2 - r);
	}

	// sort the tree
	qsort(&m_bk[0], N, sizeof(NODE), cmp_node);

	// set the initial search item
	m_imin = 0;
}

//-----------------------------------------------------------------------------

void FENNQuery::InitReference()
{
	assert(m_ps);

	int i;
	vec3d r0, r;
	int N = m_ps->Nodes();

	// pick a random point as pivot
	r0 = m_q1 = m_ps->Node(0).m_r0;

	// find the furtest node of this node
	double dmax = 0, d;
	for (i=0; i<N; ++i)
	{
		r = m_ps->Node(i).m_r0;
		d = (r - r0)*(r - r0);
		if (d > dmax)
		{
			m_q1 = r;
			dmax = d;
		}
	}

	// let's find the furthest node of this node
	r0 = m_q2 = m_q1;
	dmax = 0;
	for (i=0; i<N; ++i)
	{
		r = m_ps->Node(i).m_r0;
		d = (r - r0)*(r - r0);
		if (d > dmax)
		{
			m_q2 = r;
			dmax = d;
		}
	}

	// create the BK-"tree"
	m_bk.resize(N);
	for (i=0; i<N; ++i) 
	{
		r = m_ps->Node(i).m_r0;
		m_bk[i].i = i;
		m_bk[i].r = r;
		m_bk[i].d1 = (m_q1 - r)*(m_q1 - r);
		m_bk[i].d2 = (m_q2 - r)*(m_q2 - r);
	}

	// sort the tree
	qsort(&m_bk[0], N, sizeof(NODE), cmp_node);

	// set the initial search item
	m_imin = 0;
}

//-----------------------------------------------------------------------------

int FENNQuery::Find(vec3d x)
{
	int m_i0 = -1;
	double rmin1, rmin2, rmax1, rmax2;
	double rmin1s, rmin2s, rmax1s, rmax2s;
	double d, d1, d2, dmin;
	vec3d r;

	// set the initial search radii
	d1 = sqrt((m_q1 - x)*(m_q1 - x)); 
	rmin1 = 0;
	rmax1 = 2*d1;

	d2 = sqrt((m_q2 - x)*(m_q2 - x));
	rmin2 = 0;
	rmax2 = 2*d2;

	// check the last found item
	int imin = m_imin;
	r = m_ps->Node(imin).m_rt;
	dmin = (r - x)*(r - x);
	d = sqrt(dmin);
	
	// adjust search radii
	if (d1 - d > rmin1) rmin1 = d1 - d;
	if (d1 + d < rmax1) rmax1 = d1 + d;
	rmin1s = rmin1*rmin1;
	rmax1s = rmax1*rmax1;

	if (d2 - d > rmin2) rmin2 = d2 - d;
	if (d2 + d < rmax2) rmax2 = d2 + d;
	rmin2s = rmin2*rmin2;
	rmax2s = rmax2*rmax2;

	// find the first item that satisfies d(i, q1) >= rmin1
	int i0 = FindRadius(rmin1s);

	for (int i=i0; i<(int) m_bk.size(); ++i)
	{
		NODE& n = m_bk[i];
		if (n.d1 <= rmax1s)
		{
			if ((n.d2 >= rmin2s) && (n.d2 <= rmax2s))
			{
				r = n.r;
				d = (r - x)*(r - x);
				if (d < dmin)
				{
					dmin = d;
					d = sqrt(dmin);
					imin = n.i;

//					if (d1 - d > rmin1) rmin1 = d1 - d;
					if (d1 + d < rmax1) rmax1 = d1 + d;
//					rmin1s = rmin1*rmin1;
					rmax1s = rmax1*rmax1;

					if (d2 - d > rmin2) rmin2 = d2 - d;
					if (d2 + d < rmax2) rmax2 = d2 + d;
					rmin2s = rmin2*rmin2;
					rmax2s = rmax2*rmax2;
				}
			}
		}
		else break;
	}

/*
	// do it the hard way
	int imin = 0;
	r = m_ps->Node(imin).m_rt;
	double d0 = (r - x)*(r - x);
	for (i=0; i<m_ps->Nodes(); ++i)
	{
		r = m_ps->Node(i).m_rt;
		d = (r - x)*(r - x);
		if (d < d0)
		{
			d0 = d;
			imin = i;
		}
	}
	assert(imin == m_imin);
*/

	#pragma omp critical
	m_imin = imin;

	return imin;
}

//-----------------------------------------------------------------------------

int FENNQuery::FindReference(vec3d x)
{
	int m_i0 = -1;
	double rmin1, rmin2, rmax1, rmax2;
	double rmin1s, rmin2s, rmax1s, rmax2s;
	double d, d1, d2, dmin;
	vec3d r;

	// set the initial search radii
	d1 = sqrt((m_q1 - x)*(m_q1 - x)); 
	rmin1 = 0;
	rmax1 = 2*d1;

	d2 = sqrt((m_q2 - x)*(m_q2 - x));
	rmin2 = 0;
	rmax2 = 2*d2;

	// check the last found item
	r = m_ps->Node(m_imin).m_r0;
	dmin = (r - x)*(r - x);
	d = sqrt(dmin);
	
	// adjust search radii
	if (d1 - d > rmin1) rmin1 = d1 - d;
	if (d1 + d < rmax1) rmax1 = d1 + d;
	rmin1s = rmin1*rmin1;
	rmax1s = rmax1*rmax1;

	if (d2 - d > rmin2) rmin2 = d2 - d;
	if (d2 + d < rmax2) rmax2 = d2 + d;
	rmin2s = rmin2*rmin2;
	rmax2s = rmax2*rmax2;

	// find the first item that satisfies d(i, q1) >= rmin1
	int i0 = FindRadius(rmin1s);

	for (int i=i0; i<(int) m_bk.size(); ++i)
	{
		NODE& n = m_bk[i];
		if (n.d1 <= rmax1s)
		{
			if ((n.d2 >= rmin2s) && (n.d2 <= rmax2s))
			{
				r = n.r;
				d = (r - x)*(r - x);
				if (d < dmin)
				{
					dmin = d;
					d = sqrt(dmin);
					m_imin = n.i;

//					if (d1 - d > rmin1) rmin1 = d1 - d;
					if (d1 + d < rmax1) rmax1 = d1 + d;
//					rmin1s = rmin1*rmin1;
					rmax1s = rmax1*rmax1;

					if (d2 - d > rmin2) rmin2 = d2 - d;
					if (d2 + d < rmax2) rmax2 = d2 + d;
					rmin2s = rmin2*rmin2;
					rmax2s = rmax2*rmax2;
				}
			}
		}
		else break;
	}

/*
	// do it the hard way
	int imin = 0;
	r = m_ps->Node(imin).m_r0;
	double d0 = (r - x)*(r - x);
	for (int i=0; i<m_ps->Nodes(); ++i)
	{
		r = m_ps->Node(i).m_r0;
		d = (r - x)*(r - x);
		if (d < d0)
		{
			d0 = d;
			imin = i;
		}
	}
//	assert(imin == m_imin);
*/
	return m_imin;
}

//-----------------------------------------------------------------------------

int FENNQuery::FindRadius(double r)
{
	int N = (int)m_bk.size();
	int L = N - 1;
	int i0 = 0;
	int i1 = L;
	if (m_bk[i1].d1 < r) return N;
	int i = i1 / 2;
	do
	{
		if (m_bk[i].d1 < r)
		{
			i0 = i;
			if (m_bk[i+1].d1 >= r) { ++i; break; }
		}
		else
		{
			i1 = i;
			if ((i==0) || (m_bk[i-1].d1 < r)) break;
		}
		i = (i1 + i0) / 2;
	}
	while(i0 != i1);

	return i;
}


//-----------------------------------------------------------------------------
int findNeirestNeighbors(const std::vector<vec3d>& point, const vec3d& x, int k, std::vector<int>& closestNodes)
{
	int N0 = (int) point.size();
	if (N0 < k) k = N0;

	vector<double> dist(k, 0.0);
	closestNodes.resize(k);
	int n = 0;
	for (int i = 0; i < N0; ++i)
	{
		vec3d ri = point[i] - x;
		double L2 = ri*ri;

		if (n == 0)
		{
			closestNodes[0] = i;
			dist[0] = L2;
			n++;
		}
		else if (L2 <= dist[n - 1])
		{
			int m;
			for (m = 0; m < n; ++m)
			{
				if (L2 <= dist[m])
				{
					break;
				}
			}

			if (n < k) n++;
			for (int l = n - 1; l > m; l--)
			{
				closestNodes[l] = closestNodes[l - 1];
				dist[l] = dist[l - 1];
			}

			closestNodes[m] = i;
			dist[m] = L2;
		}
		else if (n < k)
		{
			closestNodes[n] = i;
			dist[n] = L2;
			n++;
		}
	}

	return n;
}
