#include "stdafx.h"
#include "FEDomainList.h"
#include "DumpStream.h"
#include "FEDomain.h"
#include "FEModel.h"
#include "FEMesh.h"
#include <assert.h>
using namespace std;

FEDomainList::FEDomainList()
{

}

FEDomainList::FEDomainList(FEDomainList& domList)
{
	m_dom = domList.m_dom;
}

//! Clear the domain list
void FEDomainList::Clear()
{
	m_dom.clear();
}

void FEDomainList::AddDomain(FEDomain* dom)
{
	// see if this domain is already a member of this list
	if (IsMember(dom))
	{
//		assert(false);
		return;
	}

	// it's not, so let's add it
	m_dom.push_back(dom);
}

bool FEDomainList::IsMember(const FEDomain* dom) const
{
	// loop over all the domains
	for (size_t i = 0; i < m_dom.size(); ++i)
	{
		if (m_dom[i] == dom)
		{
			// found it!
			return true;
		}
	}

	// better luck next time!
	return false;
}

//! serialization
void FEDomainList::Serialize(DumpStream& ar)
{
	if (ar.IsShallow()) return;

	if (ar.IsSaving())
	{
		int N = Domains();
		ar << N;
		for (int i = 0; i < N; ++i)
		{
			int domainId = m_dom[i]->GetID();
			ar << domainId;
		}
	}
	else
	{
		FEMesh& mesh = ar.GetFEModel().GetMesh();
		int N = 0;
		ar >> N;
		m_dom.resize(N);
		for (int i = 0; i < N; ++i)
		{
			int domId;
			ar >> domId;
			FEDomain* dom = mesh.FindDomain(domId);
			if (dom == nullptr) throw DumpStream::ReadError();
			m_dom[i] = dom;
		}
	}
}
