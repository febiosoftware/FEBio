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
#include "FEMeshShapeInterpolator.h"
#include <FECore/FEOctreeSearch.h>
#include <FECore/FEMesh.h>

//=============================================================================================
FEMeshShapeInterpolator::FEMeshShapeInterpolator(FEMesh* mesh) : m_mesh(mesh)
{
	m_os = nullptr;
}

FEMeshShapeInterpolator::~FEMeshShapeInterpolator()
{
	delete m_os;
}

bool FEMeshShapeInterpolator::Init()
{
	if (m_mesh == nullptr) return false;

	if (m_os == nullptr)
	{
		m_os = new FEOctreeSearch(m_mesh);
		if (m_os->Init() == false) return false;
	}

	int nodes = m_trgPoints.size();
	m_data.resize(nodes);
	for (int i = 0; i < nodes; ++i)
	{
		Data& di = m_data[i];
		vec3d ri = m_trgPoints[i];

		// find the element
		di.r[0] = di.r[1] = di.r[2] = 0.0;
		di.el = (FESolidElement*)m_os->FindElement(ri, di.r);
		if (di.el == nullptr)
		{
			assert(false);
			return false;
		}
	}

	return true;
}

void FEMeshShapeInterpolator::SetTargetPoints(const vector<vec3d>& trgPoints)
{
	m_trgPoints = trgPoints;
}

bool FEMeshShapeInterpolator::SetTargetPoint(const vec3d& r)
{
	m_trgPoints.clear();
	m_trgPoints.push_back(r);
	return Init();
}

bool FEMeshShapeInterpolator::Map(std::vector<double>& tval, function<double(int sourceNode)> src)
{
	// update solution
	int nodes = m_trgPoints.size();
	for (int i = 0; i < nodes; ++i)
	{
		Data& di = m_data[i];

		// update values
		double v[FEElement::MAX_NODES] = { 0 };
		for (int j = 0; j < di.el->Nodes(); ++j) v[j] = src(di.el->m_node[j]);
		double vl = di.el->evaluate(v, di.r[0], di.r[1], di.r[2]);
		tval[i] = vl;
	}

	return true;
}

double FEMeshShapeInterpolator::Map(int inode, function<double(int sourceNode)> f)
{
	Data& di = m_data[inode];
	double v[FEElement::MAX_NODES];
	for (int j = 0; j < di.el->Nodes(); ++j) v[j] = f(di.el->m_node[j]);
	double vl = di.el->evaluate(v, di.r[0], di.r[1], di.r[2]);
	return vl;
}

vec3d FEMeshShapeInterpolator::MapVec3d(int inode, function<vec3d(int sourceNode)> f)
{
	Data& di = m_data[inode];
	vec3d v[FEElement::MAX_NODES];
	for (int j = 0; j < di.el->Nodes(); ++j) v[j] = f(di.el->m_node[j]);
	vec3d vl = di.el->evaluate(v, di.r[0], di.r[1], di.r[2]);
	return vl;
}
