#pragma once
#include "FEElasticShellDomain.h"
#include "FEElasticShellDomainOld.h"

//-----------------------------------------------------------------------------
//! domain class for 3D rigid shell elements
//!
class FERigidShellDomainOld : public FEElasticShellDomainOld
{
public:
	//! constructor
	FERigidShellDomainOld(FEModel* pfem);

	//! Initialize
	bool Init() override;

	//! reset data
	void Reset() override;

public:
	//! calculates the global stiffness matrix for this domain
	void StiffnessMatrix(FESolver* psolver) override;

	//! calculates the internal forces (nothing to do)
	void InternalForces(FEGlobalVector& R) override;

	// update domain data
	void Update(const FETimeInfo& tp) override;
};

//-----------------------------------------------------------------------------
//! domain class for 3D rigid shell elements
//!
class FERigidShellDomain : public FEElasticShellDomain
{
public:
	//! constructor
	FERigidShellDomain(FEModel* pfem);

	//! Initialize
	bool Init() override;

	//! reset data
	void Reset() override;

public:
	//! calculates the global stiffness matrix for this domain
	void StiffnessMatrix(FESolver* psolver) override;

	//! calculates the internal forces (nothing to do)
	void InternalForces(FEGlobalVector& R) override;

	// update domain data
	void Update(const FETimeInfo& tp) override;
};
