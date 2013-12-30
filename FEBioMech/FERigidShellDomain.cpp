#include "stdafx.h"
#include "FERigidShellDomain.h"

//-----------------------------------------------------------------------------
bool FERigidShellDomain::Initialize(FEModel& fem)
{
	return FEShellDomain::Initialize(fem);
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
void FERigidShellDomain::UpdateStresses(FEModel &fem)
{
	// Nothing to see here. Please move on.
}
