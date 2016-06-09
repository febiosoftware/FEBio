#include "stdafx.h"
#include "FERigidSolidDomain.h"
#include "FERigidMaterial.h"
#include <FECore/FERigidBody.h>

//-----------------------------------------------------------------------------
FERigidSolidDomain::FERigidSolidDomain(FEModel* pfem) : FEElasticSolidDomain(pfem) {}

//-----------------------------------------------------------------------------
// NOTE: Although this function doesn't do anything we need it because we don't
//       want to call the FEElasticSolidDomain::Initialize function.
bool FERigidSolidDomain::Initialize(FEModel& fem)
{
	return FESolidDomain::Initialize(fem);
}

//-----------------------------------------------------------------------------
// We need to override it since the base class version will not work for rigid domains.
void FERigidSolidDomain::Reset()
{
	// nothing to reset here
}

//-----------------------------------------------------------------------------
//! Calculates the stiffness matrix for 3D rigid elements.
//! Rigid elements don't generate stress, so there is nothing to do here
void FERigidSolidDomain::StiffnessMatrix(FESolver* psolver)
{
	// Caught you looking!
}

//-----------------------------------------------------------------------------
// Rigid bodies do not generate stress so there is nothing to do here
void FERigidSolidDomain::InternalForces(FESolver* psolver, vector<double>& R)
{
	// what you looking at ?!
}

//-----------------------------------------------------------------------------
//! We don't need to update the stresses for rigid elements
//!
void FERigidSolidDomain::Update(const FETimeInfo& tp)
{
	// Nothing to see here. Please move on.
}
