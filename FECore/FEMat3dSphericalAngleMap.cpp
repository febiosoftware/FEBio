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
#include "FEMat3dSphericalAngleMap.h"

BEGIN_FECORE_CLASS(FEMat3dSphericalAngleMap, FEMat3dValuator)
	ADD_PARAMETER(m_theta, "theta");
	ADD_PARAMETER(m_phi, "phi");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEMat3dSphericalAngleMap::FEMat3dSphericalAngleMap(FEModel* pfem) : FEMat3dValuator(pfem)
{
	m_theta = 0.0;
	m_phi = 90.0;
}

//-----------------------------------------------------------------------------
bool FEMat3dSphericalAngleMap::Init()
{
	return true;
}

//-----------------------------------------------------------------------------
void FEMat3dSphericalAngleMap::SetAngles(double theta, double phi)
{
	m_theta = theta;
	m_phi = phi;
}

//-----------------------------------------------------------------------------
mat3d FEMat3dSphericalAngleMap::operator () (const FEMaterialPoint& mp)
{
	// convert from degress to radians
	const double the = m_theta(mp)*PI / 180.;
	const double phi = m_phi(mp)*PI / 180.;

	// define the first axis (i.e. the fiber vector)
	vec3d a;
	a.x = cos(the)*sin(phi);
	a.y = sin(the)*sin(phi);
	a.z = cos(phi);
    
    vec3d b;
    b.x = -sin(the);
    b.y = cos(the);
    b.z = 0;

	vec3d c;
    c.x = -cos(the)*cos(phi);
    c.y = -sin(the)*cos(phi);
    c.z = sin(phi);

	// setup the rotation matrix
	mat3d Q;
	Q[0][0] = a.x; Q[0][1] = b.x; Q[0][2] = c.x;
	Q[1][0] = a.y; Q[1][1] = b.y; Q[1][2] = c.y;
	Q[2][0] = a.z; Q[2][1] = b.z; Q[2][2] = c.z;

    return Q;
}

//-----------------------------------------------------------------------------
FEMat3dValuator* FEMat3dSphericalAngleMap::copy()
{
	FEMat3dSphericalAngleMap* map = new FEMat3dSphericalAngleMap(GetFEModel());
	map->m_theta = m_theta;
	map->m_phi = m_phi;
	return map;
}

//-----------------------------------------------------------------------------
void FEMat3dSphericalAngleMap::Serialize(DumpStream &ar)
{
	if (ar.IsShallow()) return;
	ar & m_theta & m_phi;
}
