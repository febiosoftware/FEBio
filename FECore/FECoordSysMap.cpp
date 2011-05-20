// FECoordSysMap.cpp: implementation of the FECoordSysMap class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FECoordSysMap.h"
#include "FEMesh.h"

//////////////////////////////////////////////////////////////////////
// FELocalMap
//////////////////////////////////////////////////////////////////////

FELocalMap::FELocalMap(FEMesh& mesh) : FECoordSysMap(FE_MAP_LOCAL), m_mesh(mesh)
{
	m_n[0] = 0;
	m_n[1] = 1;
	m_n[2] = 3;
}

void FELocalMap::SetLocalNodes(int n1, int n2, int n3)
{
	m_n[0] = n1;
	m_n[1] = n2;
	m_n[2] = n3;
}

mat3d FELocalMap::LocalElementCoord(FEElement& el, int n)
{
	vec3d r0[8];
	for (int i=0; i<el.Nodes(); ++i) r0[i] = m_mesh.Node(el.m_node[i]).m_r0;

	vec3d a, b, c, d;
	mat3d Q;

	a = r0[m_n[1]] - r0[m_n[0]];
	a.unit();

	if (m_n[2] != m_n[1])
	{
		d = r0[m_n[2]] - r0[m_n[0]];
	}
	else
	{
		d = vec3d(0,1,0);
		if (fabs(d*a) > 0.999) d = vec3d(1,0,0);
	}

	c = a^d;
	b = c^a;

	b.unit();
	c.unit();

	Q[0][0] = a.x; Q[0][1] = b.x; Q[0][2] = c.x;
	Q[1][0] = a.y; Q[1][1] = b.y; Q[1][2] = c.y;
	Q[2][0] = a.z; Q[2][1] = b.z; Q[2][2] = c.z;

	return Q;
}

void FELocalMap::Serialize(DumpFile& ar)
{
	if (ar.IsSaving())
	{
		ar << m_n[0] << m_n[1] << m_n[2];
	}
	else
	{
		ar >> m_n[0] >> m_n[1] >> m_n[2];
	}
}

//////////////////////////////////////////////////////////////////////
// FESphericalMap
//////////////////////////////////////////////////////////////////////

mat3d FESphericalMap::LocalElementCoord(FEElement& el, int n)
{
	vec3d r0[8];
	for (int i=0; i<el.Nodes(); ++i) r0[i] = m_mesh.Node(el.m_node[i]).m_r0;

	double* H = el.H(n);
	vec3d a;
	for (int i=0; i<el.Nodes(); ++i) a += r0[i]*H[i];
	a -= m_c;
	a.unit();

	vec3d d = r0[1] - r0[0];
	d.unit();
	if (fabs(a*d) > .99) 
	{
		d = r0[2] - r0[1];
		d.unit();
	}

	vec3d c = a^d;
	vec3d b = c^a;

	a.unit();
	b.unit();
	c.unit();

	mat3d Q;
	Q[0][0] = a.x; Q[0][1] = b.x; Q[0][2] = c.x;
	Q[1][0] = a.y; Q[1][1] = b.y; Q[1][2] = c.y;
	Q[2][0] = a.z; Q[2][1] = b.z; Q[2][2] = c.z;

	return Q;
}

void FESphericalMap::Serialize(DumpFile& ar)
{
	if (ar.IsSaving())
	{
		ar << m_c;
	}
	else
	{
		ar >> m_c;
	}
}

//////////////////////////////////////////////////////////////////////
// FEVectorMap
//////////////////////////////////////////////////////////////////////

mat3d FEVectorMap::LocalElementCoord(FEElement& el, int n)
{
	vec3d a = m_a;
	vec3d d = m_d;

	vec3d c = a^d;
	vec3d b = c^a;

	a.unit();
	b.unit();
	c.unit();

	mat3d Q;
	Q[0][0] = a.x; Q[0][1] = b.x; Q[0][2] = c.x;
	Q[1][0] = a.y; Q[1][1] = b.y; Q[1][2] = c.y;
	Q[2][0] = a.z; Q[2][1] = b.z; Q[2][2] = c.z;

	return Q;
}

void FEVectorMap::Serialize(DumpFile &ar)
{
	if (ar.IsSaving())
	{
		ar << m_a << m_d;
	}
	else
	{
		ar >> m_a >> m_d;
	}
}
