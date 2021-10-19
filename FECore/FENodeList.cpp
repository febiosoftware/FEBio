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
#include "FENodeList.h"
#include "FEMesh.h"
#include "DumpStream.h"

FENodeList::FENodeList(FEMesh* mesh) : m_mesh(mesh)
{

}

FENodeList::FENodeList(const FENodeList& nodeList)
{
	m_mesh = nodeList.m_mesh;
	m_nodes = nodeList.m_nodes;
}

FENodeList& FENodeList::operator = (const FENodeList& nodeList)
{
	m_mesh = nodeList.m_mesh;
	m_nodes = nodeList.m_nodes;
	return *this;
}

void FENodeList::Add(int n)
{
	assert(m_mesh);
	assert((n >= 0) && (n < m_mesh->Nodes()));
	m_nodes.push_back(n);
}

void FENodeList::Add(const std::vector<int>& nodeList)
{
	assert(m_mesh);
	m_nodes.insert(m_nodes.end(), nodeList.begin(), nodeList.end());
}

void FENodeList::Add(const FENodeList& nodeList)
{
	assert(m_mesh == nodeList.m_mesh);
	Add(nodeList.m_nodes);
}

void FENodeList::Clear()
{
	m_nodes.clear();
}

FENode* FENodeList::Node(int i)
{
	assert(m_mesh);
	return &m_mesh->Node(m_nodes[i]);
}

const FENode* FENodeList::Node(int i) const
{
	assert(m_mesh);
	return &m_mesh->Node(m_nodes[i]);
}

int FENodeList::Size() const
{
	return (int)m_nodes.size();
}

void FENodeList::Serialize(DumpStream& ar)
{
	if (ar.IsShallow() == false) ar & m_mesh;
	ar & m_nodes;
}

int FENodeList::GlobalToLocalID(int globalId) const
{
	for (int i = 0; i < m_nodes.size(); ++i)
	{
		if (m_nodes[i] == globalId) return i;
	}
	return -1;
}
