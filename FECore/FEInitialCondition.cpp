#include "stdafx.h"
#include "FEInitialCondition.h"
#include "FEModel.h"
#include "FEMesh.h"

FEInitialCondition::FEInitialCondition(FEModel* pfem) : FEModelComponent(FEIC_ID, pfem)
{
}

//-----------------------------------------------------------------------------
FEInitialBC::FEInitialBC(FEModel* pfem) : FEInitialCondition(pfem)
{
	m_dof = -1;
}

//-----------------------------------------------------------------------------
void FEInitialBC::Serialize(DumpFile& ar)
{
	FEInitialCondition::Serialize(ar);
	if (ar.IsSaving())
	{
		ar << m_dof;
		int nsize = m_item.size();
		ar << nsize;
		for (size_t i=0; i<nsize; ++i)
		{
			ar << m_item[i].nid << m_item[i].v;
		}
	}
	else
	{
		ar >> m_dof;
		int nsize = 0;
		ar >> nsize;
		m_item.resize(nsize);
		for (size_t i=0; i<nsize; ++i)
		{
			ar >> m_item[i].nid >> m_item[i].v;
		}
	}
}

//-----------------------------------------------------------------------------
void FEInitialBC::Activate()
{
	FEInitialCondition::Activate();
	assert(m_dof >= 0);
	if (m_dof == -1) return;
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	for (size_t i=0; i<m_item.size(); ++i)
	{
		FENode& node = mesh.Node(m_item[i].nid);
		node.set(m_dof, m_item[i].v);
	}
}

//-----------------------------------------------------------------------------
void FEInitialBCVec3D::Serialize(DumpFile& ar)
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
