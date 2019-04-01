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
#include "FENodeSet.h"
#include "FEMesh.h"
#include "DumpStream.h"

//=============================================================================
// FENodeSet
//-----------------------------------------------------------------------------
FENodeSet::FENodeSet() : m_pmesh(0), m_nID(-1)
{
}

//-----------------------------------------------------------------------------
FENodeSet::FENodeSet(FEMesh* pm) : m_pmesh(pm), m_nID(-1)
{
}

//-----------------------------------------------------------------------------
FENodeSet::FENodeSet(const FENodeSet& n)
{
	m_pmesh = n.m_pmesh;
	m_Node = n.m_Node;
	SetName(n.GetName());
}

//-----------------------------------------------------------------------------
void FENodeSet::operator = (const FENodeSet& n)
{
	m_pmesh = n.m_pmesh;
	m_Node = n.m_Node;
	SetName(n.GetName());
}

//-----------------------------------------------------------------------------
void FENodeSet::create(int n)
{
	assert(n);
	m_Node.resize(n);
}

//-----------------------------------------------------------------------------
void FENodeSet::add(int id)
{
	m_Node.push_back(id);
}

//-----------------------------------------------------------------------------
void FENodeSet::add(const std::vector<int>& ns)
{
	int n0 = (int)m_Node.size();
	int n1 = (int)ns.size();
	int N = n0 + n1;
	m_Node.resize(N);
	for (int i = 0; i<n1; ++i) m_Node[n0 + i] = ns[i];
}

//-----------------------------------------------------------------------------
void FENodeSet::add(const FENodeSet& ns)
{
	int n0 = (int)m_Node.size();
	int n1 = ns.size();
	int N = n0 + n1;
	m_Node.resize(N);
	for (int i = 0; i<n1; ++i) m_Node[n0 + i] = ns[i];
}

//-----------------------------------------------------------------------------
FENode* FENodeSet::Node(int i)
{
	return &m_pmesh->Node(m_Node[i]);
}

//-----------------------------------------------------------------------------
const FENode* FENodeSet::Node(int i) const
{
	return &m_pmesh->Node(m_Node[i]);
}

//-----------------------------------------------------------------------------
void FENodeSet::Serialize(DumpStream& ar)
{
	FEItemList::Serialize(ar);
	if (ar.IsShallow()) return;
	ar & m_nID;
	ar & m_Node;
}

//-----------------------------------------------------------------------------
void FENodeSet::SaveClass(DumpStream& ar, FENodeSet* p)
{
	FEMesh* m = p->GetMesh();
	ar << m;
}

//-----------------------------------------------------------------------------
FENodeSet* FENodeSet::LoadClass(DumpStream& ar, FENodeSet* p)
{
	FEMesh* m = nullptr;
	ar >> m; assert(m);
	return new FENodeSet(m);
}
