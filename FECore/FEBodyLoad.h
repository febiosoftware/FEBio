#pragma once
#include "FEModelComponent.h"
#include "FEDomain.h"
#include "FEDomainList.h"

//-----------------------------------------------------------------------------
// forward declaration of FEModel class
class FEModel;

//-----------------------------------------------------------------------------
//! Base class for body-loads
class FECORE_API FEBodyLoad : public FEModelComponent
{
public:
	FEBodyLoad(FEModel* pfem);
	virtual ~FEBodyLoad();

	//! return number of domains this load is applied to
	int Domains() const;

	//! return a domain 
	FEDomain* Domain(int i);

	//! add a domain to which to apply this load
	void SetDomainList(FEElementSet* elset);

	// get the domain list
	FEDomainList& GetDomaintList() { return m_dom; }

public:
	//! initialization
	virtual bool Init();

	//! update
	virtual void Update();

private:
	FEDomainList	m_dom;	//!< list of domains to which to apply the body load
};
