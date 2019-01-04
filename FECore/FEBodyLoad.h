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
	DECLARE_SUPER_CLASS(FEBODYLOAD_ID);

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
	FEDomainList& GetDomainList() { return m_dom; }

	// Evaluate residual
	virtual void Residual(const FETimeInfo& tp, FEGlobalVector& R);

public:
	//! initialization
	virtual bool Init();

private:
	FEDomainList	m_dom;	//!< list of domains to which to apply the body load
};
