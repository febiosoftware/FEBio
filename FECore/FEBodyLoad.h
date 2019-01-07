#pragma once
#include "FEModelComponent.h"
#include "FEDomain.h"
#include "FEDomainList.h"

//-----------------------------------------------------------------------------
// forward declaration of FEModel class
class FEModel;
class FELinearSystem;

//-----------------------------------------------------------------------------
//! Base class for body-loads
class FECORE_API FEBodyLoad : public FEModelComponent
{
	FECORE_SUPER_CLASS

public:
	FEBodyLoad(FEModel* pfem);
	virtual ~FEBodyLoad();

	//! initialization
	bool Init() override;

public:
	//! return number of domains this load is applied to
	int Domains() const;

	//! return a domain 
	FEDomain* Domain(int i);

	//! add a domain to which to apply this load
	void SetDomainList(FEElementSet* elset);

	//! get the domain list
	FEDomainList& GetDomainList();

public: // This should be overridden by derived classes

	//! Evaluate force vector
	virtual void ForceVector(FEGlobalVector& R);

	//! evaluate stiffness matrix
	virtual void StiffnessMatrix(FELinearSystem& S);

private:
	FEDomainList	m_dom;	//!< list of domains to which to apply the body load
};
