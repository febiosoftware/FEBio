#include "stdafx.h"
#include "FENodeSet.h"
#include "FEMesh.h"

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
	m_name = n.m_name;
}

//-----------------------------------------------------------------------------
void FENodeSet::operator = (const FENodeSet& n)
{
	m_pmesh = n.m_pmesh;
	m_Node = n.m_Node;
	m_name = n.m_name;
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
void FENodeSet::SetName(const std::string& name)
{
	m_name = name;
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
	if (ar.IsSaving())
	{
		ar << m_nID;
		ar << m_name;
		ar << m_Node;
	}
	else
	{
		ar >> m_nID;
		ar >> m_name;
		ar >> m_Node;
	}
}
