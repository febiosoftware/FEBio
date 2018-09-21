#include "stdafx.h"
#include "FEDomainList.h"
#include <assert.h>
using namespace std;

FEDomainList::FEDomainList()
{

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
		assert(false);
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

