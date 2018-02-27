#include "stdafx.h"
#include "FERigidShellDomain.h"

//-----------------------------------------------------------------------------
FERigidShellDomainOld::FERigidShellDomainOld(FEModel* pfem) : FEElasticShellDomainOld(pfem) {}

//-----------------------------------------------------------------------------
// NOTE: Although this function doesn't do anything, we need it since 
//       for rigid shell domains we don't want to call the FEElasticShellDomain::Initialize member.
bool FERigidShellDomainOld::Initialize()
{
	// just call the base class
	return FEShellDomainOld::Initialize();
}

//-----------------------------------------------------------------------------
// We need to override it since the base class version will not work for rigid domains.
void FERigidShellDomainOld::Reset()
{
	// nothing here
}

//-----------------------------------------------------------------------------
//! Calculate stiffness contributions for rigid shells.
//! Since rigid elements don't generate stress, we don't need to do
//! anything here.
void FERigidShellDomainOld::StiffnessMatrix(FESolver* psolver)
{
	// Caught you looking!
}


//-----------------------------------------------------------------------------
//! calculate residual forces for rigid shells
//!
void FERigidShellDomainOld::InternalForces(FEGlobalVector& R)
{
	// Nothing to do.
}

//-----------------------------------------------------------------------------
//! update stresses for rigid shells.
//!
void FERigidShellDomainOld::Update(const FETimeInfo& tp)
{
	// Nothing to see here. Please move on.
}

//===========================================================================================================

//-----------------------------------------------------------------------------
FERigidShellDomain::FERigidShellDomain(FEModel* pfem) : FEElasticShellDomain(pfem) {}

//-----------------------------------------------------------------------------
// NOTE: Although this function doesn't do anything, we need it since 
//       for rigid shell domains we don't want to call the FEElasticShellDomain::Initialize member.
bool FERigidShellDomain::Initialize()
{
	// just call the base class
	return FEShellDomain::Initialize();
}

//-----------------------------------------------------------------------------
// We need to override it since the base class version will not work for rigid domains.
void FERigidShellDomain::Reset()
{
	// nothing here
}

//-----------------------------------------------------------------------------
//! Calculate stiffness contributions for rigid shells.
//! Since rigid elements don't generate stress, we don't need to do
//! anything here.
void FERigidShellDomain::StiffnessMatrix(FESolver* psolver)
{
	// Caught you looking!
}


//-----------------------------------------------------------------------------
//! calculate residual forces for rigid shells
//!
void FERigidShellDomain::InternalForces(FEGlobalVector& R)
{
	// Nothing to do.
}

//-----------------------------------------------------------------------------
//! update stresses for rigid shells.
//!
void FERigidShellDomain::Update(const FETimeInfo& tp)
{
	// Nothing to see here. Please move on.
}
