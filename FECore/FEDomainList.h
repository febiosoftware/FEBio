#pragma once
#include <vector>

class FEDomain;

//-----------------------------------------------------------------------------
// Class that represents a list of domains. 
class FEDomainList
{
public:
	//! constructor
	FEDomainList();

	//! Clear the domain list
	void Clear();

	//! Add a domain to the list
	void AddDomain(FEDomain* dom);

	//! See if a domain is a member of this list
	bool IsMember(const FEDomain* dom) const;

	//! See if the list is empty
	bool IsEmpty() const { return m_dom.empty(); }

	//! Return number of domains in list
	int Domains() const { return (int) m_dom.size();  }

	//! Return a domain
	FEDomain* GetDomain(int i) { return m_dom[i]; }

private:
	std::vector<FEDomain*>	m_dom;		// the actual list of domains
};
