#include "stdafx.h"
#include "FEInitialCondition.h"
#include "FEModel.h"
#include "FEMesh.h"

FEInitialCondition::FEInitialCondition(FEModel* pfem) : FEModelComponent(FEIC_ID, pfem)
{
}

//-----------------------------------------------------------------------------
void FEInitialVelocity::Serialize(DumpFile& ar)
{
	if (ar.IsSaving())
	{
		for (size_t i=0; i<m_item.size(); ++i)
		{
			ar << m_item[i].nid << m_item[i].v0;
		}
	}
	else
	{
		for (size_t i=0; i<m_item.size(); ++i)
		{
			ar >> m_item[i].nid >> m_item[i].v0;
		}
	}
}

//-----------------------------------------------------------------------------
void FEInitialVelocity::Activate()
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	for (size_t i=0; i<m_item.size(); ++i)
	{
		FENode& node = mesh.Node(m_item[i].nid);
		node.m_vt = m_item[i].v0;
	}
}

//-----------------------------------------------------------------------------
void FEInitialPressure::Serialize(DumpFile& ar)
{
	if (ar.IsSaving())
	{
		for (size_t i=0; i<m_item.size(); ++i)
		{
			ar << m_item[i].nid << m_item[i].p0;
		}
	}
	else
	{
		for (size_t i=0; i<m_item.size(); ++i)
		{
			ar >> m_item[i].nid >> m_item[i].p0;
		}
	}
}

//-----------------------------------------------------------------------------
void FEInitialPressure::Activate()
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	for (size_t i=0; i<m_item.size(); ++i)
	{
		FENode& node = mesh.Node(m_item[i].nid);
		node.m_pt = m_item[i].p0;
	}
}

//-----------------------------------------------------------------------------
void FEInitialTemperature::Serialize(DumpFile& ar)
{
	if (ar.IsSaving())
	{
		for (size_t i=0; i<m_item.size(); ++i)
		{
			ar << m_item[i].nid << m_item[i].T0;
		}
	}
	else
	{
		for (size_t i=0; i<m_item.size(); ++i)
		{
			ar >> m_item[i].nid >> m_item[i].T0;
		}
	}
}

//-----------------------------------------------------------------------------
void FEInitialTemperature::Activate()
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	for (size_t i=0; i<m_item.size(); ++i)
	{
		FENode& node = mesh.Node(m_item[i].nid);
		node.m_T = node.m_T0 = m_item[i].T0;
	}
}

//-----------------------------------------------------------------------------
void FEInitialConcentration::Serialize(DumpFile& ar)
{
	if (ar.IsSaving())
	{
		ar << m_nsol;
		for (size_t i=0; i<m_item.size(); ++i)
		{
			ar << m_item[i].nid << m_item[i].c0;
		}
	}
	else
	{
		ar >> m_nsol;
		for (size_t i=0; i<m_item.size(); ++i)
		{
			ar >> m_item[i].nid >> m_item[i].c0;
		}
	}
}

//-----------------------------------------------------------------------------
void FEInitialConcentration::Activate()
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	for (size_t i=0; i<m_item.size(); ++i)
	{
		FENode& node = mesh.Node(m_item[i].nid);
		node.m_ct[m_nsol] = m_item[i].c0;
	}
}

//-----------------------------------------------------------------------------
void FEInitialDilatation::Serialize(DumpFile& ar)
{
    if (ar.IsSaving())
    {
        for (size_t i=0; i<m_item.size(); ++i)
        {
            ar << m_item[i].nid << m_item[i].e0;
        }
    }
    else
    {
        for (size_t i=0; i<m_item.size(); ++i)
        {
            ar >> m_item[i].nid >> m_item[i].e0;
        }
    }
}

//-----------------------------------------------------------------------------
void FEInitialDilatation::Activate()
{
    FEModel& fem = *GetFEModel();
    FEMesh& mesh = fem.GetMesh();
    for (size_t i=0; i<m_item.size(); ++i)
    {
        FENode& node = mesh.Node(m_item[i].nid);
        node.m_et = m_item[i].e0;
    }
}
