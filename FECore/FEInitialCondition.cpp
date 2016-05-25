#include "stdafx.h"
#include "FEInitialCondition.h"
#include "FEModel.h"
#include "FEMesh.h"

BEGIN_PARAMETER_LIST(FEInitialBC, FEInitialCondition)
	ADD_PARAMETER(m_data, FE_PARAM_DATA_ARRAY, "value");
END_PARAMETER_LIST();

FEInitialCondition::FEInitialCondition(FEModel* pfem) : FEModelComponent(FEIC_ID, pfem)
{
}

//-----------------------------------------------------------------------------
FEInitialBC::FEInitialBC(FEModel* pfem) : FEInitialCondition(pfem), m_data(FE_DOUBLE)
{
	m_dof = -1;
	m_data.set(0.0);
}

//-----------------------------------------------------------------------------
void FEInitialBC::Serialize(DumpStream& ar)
{
	FEInitialCondition::Serialize(ar);
	if (ar.IsSaving())
	{
		ar << m_dof;
		int nsize = m_item.size();
		ar << nsize;
		for (size_t i=0; i<nsize; ++i) ar << m_item[i];
	}
	else
	{
		ar >> m_dof;
		int nsize = 0;
		ar >> nsize;
		m_item.resize(nsize);
		for (size_t i=0; i<nsize; ++i) ar >> m_item[i];
	}
}

//-----------------------------------------------------------------------------
void FEInitialBC::SetNodes(const FENodeSet& set)
{
	int N = set.size();
	m_item.resize(N);
	for (int i=0; i<N; ++i) m_item[i] = set[i];
	m_data.resize(N);
}

//-----------------------------------------------------------------------------
void FEInitialBC::Add(int node, double value)
{
	m_item.push_back(node);
	m_data.push_back(value);
}

//-----------------------------------------------------------------------------
void FEInitialBC::Activate()
{
	FEInitialCondition::Activate();
	assert(m_dof >= 0);
	if (m_dof == -1) return;
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	int N = m_item.size();
	for (size_t i=0; i<N; ++i)
	{
		FENode& node = mesh.Node(m_item[i]);
		node.set(m_dof, m_data.get<double>(i));
	}
}

//-----------------------------------------------------------------------------
void FEInitialBCVec3D::Serialize(DumpStream& ar)
{
	FEInitialCondition::Serialize(ar);
	if (ar.IsSaving())
	{
		ar << m_dof[0] << m_dof[1] << m_dof[2];
		int nsize = m_item.size();
		ar << nsize;
		for (size_t i=0; i<nsize; ++i)
		{
			ar << m_item[i].nid << m_item[i].v0;
		}
	}
	else
	{
		ar >> m_dof[0] >> m_dof[1] >> m_dof[2];
		int nsize = 0;
		ar >> nsize;
		m_item.resize(nsize);
		for (size_t i=0; i<nsize; ++i)
		{
			ar >> m_item[i].nid >> m_item[i].v0;
		}
	}
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
