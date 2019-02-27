#include "stdafx.h"
#include "FEInitialCondition.h"
#include "FEModel.h"
#include "FEMesh.h"

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
