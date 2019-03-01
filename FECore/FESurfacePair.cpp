#include "stdafx.h"
#include "FESurfacePair.h"
#include "FEFacetSet.h"
#include "FEMesh.h"
#include "DumpStream.h"

//--------------------------------------------------------
FESurfacePair::FESurfacePair(FEMesh* pm) : m_mesh(pm)
{
	m_master = 0;
	m_slave = 0;
}

void FESurfacePair::SetName(const std::string& name)
{
	m_name = name;
}

const std::string& FESurfacePair::GetName() const
{
	return m_name;
}

FEFacetSet* FESurfacePair::GetMasterSurface()
{
	return m_master;
}

void FESurfacePair::SetMasterSurface(FEFacetSet* pf)
{
	m_master = pf;
}

FEFacetSet* FESurfacePair::GetSlaveSurface()
{
	return m_slave;
}

void FESurfacePair::SetSlaveSurface(FEFacetSet* pf)
{
	m_slave = pf;
}

void FESurfacePair::Serialize(DumpStream& ar)
{
	if (ar.IsSaving())
	{
		ar << m_name;

		if (m_master)
		{
			ar << (int)1;
			ar << m_master->GetName();
		}
		else ar << (int)0;

		if (m_slave)
		{
			ar << (int)1;
			ar << m_slave->GetName();
		}
		else ar << (int)0;
	}
	else
	{
		ar >> m_name;

		// NOTE: This assumes that facet sets have already been serialized!!
		int flag = 0;
		ar >> flag;
		if (flag == 1)
		{
			string name;
			ar >> name;
			m_master = m_mesh->FindFacetSet(name); assert(m_master);
		}
		else m_master = 0;

		ar >> flag;
		if (flag == 1)
		{
			string name;
			ar >> name;
			m_slave = m_mesh->FindFacetSet(name); assert(m_slave);
		}
		else m_slave = 0;
	}
}
