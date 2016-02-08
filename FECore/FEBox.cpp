#include "stdafx.h"
#include "FEBox.h"
#include "FEMesh.h"
#include "FEDomain.h"

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

FEBox::FEBox(const FEDomain& dom)
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
