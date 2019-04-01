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
#include "FEInitialCondition.h"
#include "FEModel.h"
#include "FEMesh.h"
#include "DumpStream.h"

REGISTER_SUPER_CLASS(FEInitialCondition, FEIC_ID);

BEGIN_FECORE_CLASS(FEInitialBC, FEInitialCondition)
	ADD_PARAMETER(m_data, "value");
END_FECORE_CLASS();

FEInitialCondition::FEInitialCondition(FEModel* pfem) : FEModelComponent(pfem)
{
}

//-----------------------------------------------------------------------------
FEInitialBC::FEInitialBC(FEModel* pfem) : FEInitialCondition(pfem), m_data(FE_DOUBLE)
{
	m_dof = -1;
}

//-----------------------------------------------------------------------------
void FEInitialBC::Serialize(DumpStream& ar)
{
	FEInitialCondition::Serialize(ar);
	if (ar.IsShallow()) return;
	ar & m_dof;
	ar & m_item;
}

//-----------------------------------------------------------------------------
void FEInitialBC::SetNodes(const FENodeSet& set)
{
	int N = set.size();
	m_item.resize(N);
	for (int i=0; i<N; ++i) m_item[i] = set[i];
	m_data.Create(N, 0.0);
}

//-----------------------------------------------------------------------------
void FEInitialBC::Add(int node, double value)
{
	m_item.push_back(node);
	m_data.Add(value);
}

//-----------------------------------------------------------------------------
void FEInitialBC::Activate()
{
	FEInitialCondition::Activate();
	assert(m_dof >= 0);
	if (m_dof == -1) return;
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	int N = (int)m_item.size();
	for (size_t i=0; i<N; ++i)
	{
		FENode& node = mesh.Node(m_item[i]);
		node.set(m_dof, m_data.getValue((int)i));
	}
}

//-----------------------------------------------------------------------------
void FEInitialBCVec3D::ITEM::Serialize(DumpStream& ar)
{
	ar & nid & v0;
}

//-----------------------------------------------------------------------------
void FEInitialBCVec3D::Serialize(DumpStream& ar)
{
	FEInitialCondition::Serialize(ar);
	if (ar.IsShallow()) return;
	ar & m_dof;
	ar & m_item;
}

//-----------------------------------------------------------------------------
void FEInitialBCVec3D::Activate()
{
	assert((m_dof[0]>=0)&&(m_dof[1]>=0)&&(m_dof[2]>=0));
	FEInitialCondition::Activate();
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	for (size_t i=0; i<m_item.size(); ++i)
	{
		FENode& node = mesh.Node(m_item[i].nid);
		node.set_vec3d(m_dof[0], m_dof[1], m_dof[2], m_item[i].v0);
	}
}
