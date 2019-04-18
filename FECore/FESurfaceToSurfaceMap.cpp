/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, Columbia University, and others.

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
#include "FESurfaceToSurfaceMap.h"
#include "FEModel.h"
#include "FEMesh.h"
#include "FESurface.h"

BEGIN_FECORE_CLASS(FESurfaceToSurfaceMap, FEDataGenerator)
	ADD_PROPERTY(m_func, "function");
	ADD_PROPERTY(m_surf1, "bottom_surface", FEProperty::Reference);
	ADD_PROPERTY(m_surf2, "top_surface"   , FEProperty::Reference);
END_FECORE_CLASS();

FESurfaceToSurfaceMap::FESurfaceToSurfaceMap(FEModel* fem) : FEDataGenerator(fem)
{
	m_ccp = 0;
	m_npr = 0;
	m_func = 0;

	m_surf1 = nullptr;
	m_surf2 = nullptr;

	m_binverted = false;
}

FESurfaceToSurfaceMap::~FESurfaceToSurfaceMap()
{
	if (m_ccp) delete m_ccp;
	if (m_npr) delete m_npr;
}

bool FESurfaceToSurfaceMap::Init()
{
	FEModel* fem = GetFEModel();
	if (fem == 0) return false;
	if (m_func == 0) return false;
	if ((m_surf1 == 0) || (m_surf2 == 0)) return false;

	// we need to invert the second surface, otherwise the normal projections won't work
	if (m_binverted == false)
	{
		m_surf2->Invert();
		m_binverted = true;
	}

	// initialize projections
	if (m_ccp == nullptr)
	{
		m_ccp = new FEClosestPointProjection(*m_surf1);
		if (m_ccp->Init() == false) return false;
	}

	if (m_npr == nullptr)
	{
		FEMesh& mesh = fem->GetMesh();
		double R = mesh.GetBoundingBox().radius();
		m_npr = new FENormalProjection(*m_surf2);
		m_npr->SetSearchRadius(R);
		m_npr->SetTolerance(0.001);
		m_npr->Init();
	}

	return true;
}

void FESurfaceToSurfaceMap::value(const vec3d& x, double& data)
{
	vec3d r(x);

	// project x onto surface 1 using CCP
	vec3d q1(0,0,0);
	vec2d r1;
	m_ccp->Project(r, q1, r1);

	// project x onto surface 2 using ray-intersection
	vec3d N = r - q1; N.unit();
	vec3d q2 = m_npr->Project(q1, N);
	double L2 = (q2 - q1)*(q2 - q1);
	if (L2 == 0.0) L2 = 1.0;

	// find the fractional distance
	double w = ((x - q1)*(q2 - q1))/L2;

	// evaluate the function
	data = m_func->value(w);
}
