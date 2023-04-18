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
#include "FEPointBodyForce.h"
#include <FECore/FEMesh.h>

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEPointBodyForce, FEBodyForce);
	ADD_PARAMETER(m_a, "a");
	ADD_PARAMETER(m_b, "b");
	ADD_PARAMETER(m_rc, "rc");
	ADD_PARAMETER(m_inode, "node");
	ADD_PARAMETER(m_brigid, "rigid");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEPointBodyForce::FEPointBodyForce(FEModel* pfem) : FEBodyForce(pfem)
{
	m_pel = 0; 
	m_brigid = true; 
	m_inode = -1; 
}

//-----------------------------------------------------------------------------
vec3d FEPointBodyForce::force(FEMaterialPoint& mp)
{
	vec3d x = mp.m_rt;
	vec3d n = x - m_rc;
	double l = n.unit();

	double g = m_a*exp(-m_b*l);
	return n*g;
}

//-----------------------------------------------------------------------------
mat3ds FEPointBodyForce::stiffness(FEMaterialPoint &mp)
{
	vec3d x = mp.m_rt;
	vec3d n = x - m_rc;
	double l = n.unit();

	mat3ds k;
	if (l == 0.0)
	{
		k.zero();
	}
	else
	{
		double g = m_a*exp(-m_b*l);

		mat3ds nxn = dyad(n);
		mat3ds I = mat3dd(1.0);

		k = (nxn*m_b - (I - nxn)/l)*g;
	}

	return k;
}

//-----------------------------------------------------------------------------
void FEPointBodyForce::Serialize(DumpStream &ar)
{
	FEBodyForce::Serialize(ar);
	ar & m_a & m_b & m_rc;
	ar & m_inode & m_brigid;
}

//-----------------------------------------------------------------------------
bool FEPointBodyForce::Init()
{
	FEMesh& m = GetMesh();
	if (m_inode == -1)
	{
		if (!m_brigid)
		{
			// find the element in which point r0 lies
			m_pel = m.FindSolidElement(m_rc, m_rs);
		}
		else m_pel = 0;
	}
	else 
	{
		m_rc = m.Node(m_inode).m_r0;
	}

	return true;
}

//-----------------------------------------------------------------------------
// Update the position of the body force
void FEPointBodyForce::Update()
{
	if (m_inode == -1)
	{
		if (m_pel)
		{
			FEMesh& m = GetMesh();
			vec3d x[FEElement::MAX_NODES];
			for (int i=0; i<8; ++i) x[i] = m.Node(m_pel->m_node[i]).m_rt;

			double* r = m_rs;
			double H[FEElement::MAX_NODES];
			H[0] = 0.125*(1 - r[0])*(1 - r[1])*(1 - r[2]);
			H[1] = 0.125*(1 + r[0])*(1 - r[1])*(1 - r[2]);
			H[2] = 0.125*(1 + r[0])*(1 + r[1])*(1 - r[2]);
			H[3] = 0.125*(1 - r[0])*(1 + r[1])*(1 - r[2]);
			H[4] = 0.125*(1 - r[0])*(1 - r[1])*(1 + r[2]);
			H[5] = 0.125*(1 + r[0])*(1 - r[1])*(1 + r[2]);
			H[6] = 0.125*(1 + r[0])*(1 + r[1])*(1 + r[2]);
			H[7] = 0.125*(1 - r[0])*(1 + r[1])*(1 + r[2]);

			m_rc = vec3d(0,0,0);
			for (int i=0; i<8; ++i) m_rc += x[i]*H[i];
		}
	}
	else
	{
		FEMesh& m = GetMesh();
		m_rc = m.Node(m_inode).m_rt;
	}
}
