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
#include "FEBox.h"
#include "FEMesh.h"
#include "FEMeshPartition.h"

FEBox::FEBox()
{
}

FEBox::FEBox(const vec3d& r0, const vec3d& r1) : m_r0(r0), m_r1(r1)
{
}

FEBox::FEBox(const FEMesh& mesh)
{
	m_r0 = m_r1 = vec3d(0,0,0);
	int N = mesh.Nodes();
	if (N > 0)
	{
		m_r0 = m_r1 = mesh.Node(0).m_rt;
		for (int i=1; i<N; ++i)
		{
			const FENode& ni = mesh.Node(i);
			add(ni.m_rt);
		}
	}
}

FEBox::FEBox(const FEMeshPartition& dom)
{
	m_r0 = m_r1 = vec3d(0,0,0);
	int N = dom.Nodes();
	if (N > 0)
	{
		m_r0 = m_r1 = dom.Node(0).m_rt;
		for (int i=1; i<N; ++i)
		{
			const FENode& ni = dom.Node(i);
			add(ni.m_rt);
		}
	}
}

double FEBox::maxsize()
{
	double dx = fabs(width());
	double dy = fabs(height());
	double dz = fabs(depth());

	double D = dx;
	if (dy > D) D = dy;
	if (dz > D) D = dz;

	return D;
}

void FEBox::add(const vec3d& r)
{
	if (r.x < m_r0.x) m_r0.x = r.x; if (r.x > m_r1.x) m_r1.x = r.x;
	if (r.y < m_r0.y) m_r0.y = r.y; if (r.y > m_r1.y) m_r1.y = r.y;
	if (r.z < m_r0.z) m_r0.z = r.z; if (r.z > m_r1.z) m_r1.z = r.z;
}
