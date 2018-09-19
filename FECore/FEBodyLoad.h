#pragma once
#include "FEModelComponent.h"
#include "FEDomain.h"

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
	virtual void AddDomain(FEDomain* dom);

public:
	//! initialization
	virtual bool Init();

	//! update
	virtual void Update();

private:
	vector<FEDomain*>	m_dom;	//!< list of domains to which to apply the body load
};
