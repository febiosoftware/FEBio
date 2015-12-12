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
	FEInitialCondition::Serialize(ar);
	if (ar.IsSaving())
	{
		int nsize = m_item.size();
		ar << nsize;
		for (size_t i=0; i<nsize; ++i)
		{
			ar << m_item[i].nid << m_item[i].v0;
		}
	}
	else
	{
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
	FEInitialCondition::Serialize(ar);
	if (ar.IsSaving())
	{
		int nsize = m_item.size();
		ar << nsize;
		for (size_t i=0; i<nsize; ++i)
		{
			ar << m_item[i].nid << m_item[i].p0;
		}
	}
	else
	{
		int nsize = 0;
		ar >> nsize;
		m_item.resize(nsize);
		for (size_t i=0; i<nsize; ++i)
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
		node.set(DOF_P, m_item[i].p0);
	}
}

//-----------------------------------------------------------------------------
void FEInitialTemperature::Serialize(DumpFile& ar)
{
	FEInitialCondition::Serialize(ar);
	if (ar.IsSaving())
	{
		int nsize = m_item.size();
		ar << nsize;
		for (size_t i=0; i<nsize; ++i)
		{
			ar << m_item[i].nid << m_item[i].T0;
		}
	}
	else
	{
		int nsize = 0;
		ar >> nsize;
		m_item.resize(nsize);
		for (size_t i=0; i<nsize; ++i)
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
		node.set(DOF_T, m_item[i].T0);
	}
}

//-----------------------------------------------------------------------------
void FEInitialConcentration::Serialize(DumpFile& ar)
{
	FEInitialCondition::Serialize(ar);
	if (ar.IsSaving())
	{
		int nsize = m_item.size();
		ar << m_nsol;
		ar << nsize;
		for (size_t i=0; i<nsize; ++i)
		{
			ar << m_item[i].nid << m_item[i].c0;
		}
	}
	else
	{
		int nsize = 0;
		ar >> m_nsol;
		ar >> nsize;
		m_item.resize(nsize);
		for (size_t i=0; i<nsize; ++i)
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
		node.set(DOF_C + m_nsol, m_item[i].c0);
	}
}

//-----------------------------------------------------------------------------
void FEInitialDilatation::Serialize(DumpFile& ar)
{
	FEInitialCondition::Serialize(ar);
    if (ar.IsSaving())
    {
		int nsize = m_item.size();
		ar << nsize;
        for (size_t i=0; i<nsize; ++i)
        {
            ar << m_item[i].nid << m_item[i].e0;
        }
    }
    else
    {
		int nsize = 0;
		ar >> nsize;
		m_item.resize(nsize);
        for (size_t i=0; i<nsize; ++i)
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
        node.set(DOF_E, m_item[i].e0);
    }
}
