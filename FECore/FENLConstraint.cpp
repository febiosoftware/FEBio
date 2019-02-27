#include "stdafx.h"
#include "FENLConstraint.h"

REGISTER_SUPER_CLASS(FENLConstraint, FENLCONSTRAINT_ID);

//-----------------------------------------------------------------------------
FENLConstraint::FENLConstraint(FEModel* pfem) : FEModelComponent(pfem)
{
	static int ncount = 1;
	SetID(ncount++);
}

//-----------------------------------------------------------------------------
FENLConstraint::~FENLConstraint()
{
}
