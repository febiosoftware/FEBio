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
